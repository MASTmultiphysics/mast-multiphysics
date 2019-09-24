/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// C++ includes
#include <iomanip>

// MAST includes
#include "examples/base/input_wrapper.h"
#include "examples/fluid/meshing/cylinder.h"
#include "level_set/level_set_discipline.h"
#include "level_set/level_set_system_initialization.h"
#include "level_set/level_set_eigenproblem_assembly.h"
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_perimeter_output.h"
#include "level_set/level_set_boundary_velocity.h"
#include "level_set/indicator_function_constrain_dofs.h"
#include "level_set/level_set_constrain_dofs.h"
#include "level_set/level_set_intersection.h"
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/ks_stress_output.h"
#include "elasticity/smooth_ramp_stress_output.h"
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/compliance_output.h"
#include "elasticity/structural_system_initialization.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/constant_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/slepc_eigen_solver.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/npsol_optimization_interface.h"
#include "optimization/function_evaluation.h"


// libMesh includes
#include "libmesh/fe_type.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"


void
_optim_obj(int*    mode,
           int*    n,
           double* x,
           double* f,
           double* g,
           int*    nstate);
void
_optim_con(int*    mode,
           int*    ncnln,
           int*    n,
           int*    ldJ,
           int*    needc,
           double* x,
           double* c,
           double* cJac,
           int*    nstate);

//
// BEGIN_TRANSLATE 2D Level-set topology optimization
//
//   \tableofcontents
//
//  This example computes the optimal topology of a structure subject to
//  specified boundary conditions (Dirichlet and Neumann). A level-set function
//  is used to implicitly define the geometry inside a mesh using the
//  immersed boundary approach.
//
//  Level Set Mesh Function
class PhiMeshFunction:
public MAST::FieldFunction<Real> {
public:
    PhiMeshFunction():
    MAST::FieldFunction<Real>("phi"), _phi(nullptr) { }
    virtual ~PhiMeshFunction(){ if (_phi) delete _phi;}
    
    void init(MAST::SystemInitialization& sys, const libMesh::NumericVector<Real>& sol) {
        if (!_phi) _phi = new MAST::MeshFieldFunction(sys, "phi");
        else _phi->clear();
        _phi->init(sol);
    }
    
    MAST::MeshFieldFunction& get_mesh_function() {return *_phi;}
    
    virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
        libmesh_assert(_phi);
        RealVectorX v1;
        (*_phi)(p, t, v1);
        v = v1(0);
    }
    
protected:
    MAST::MeshFieldFunction *_phi;
};


class ElementParameterDependence:
public MAST::AssemblyBase::ElemParameterDependence {
public:
    ElementParameterDependence(const MAST::FilterBase& filter):
    MAST::AssemblyBase::ElemParameterDependence(true), _filter(filter) {}
    virtual ~ElementParameterDependence() {}
    virtual bool if_elem_depends_on_parameter(const libMesh::Elem& e,
                                              const MAST::FunctionBase& p) const {
        const MAST::LevelSetParameter
        &p_ls = dynamic_cast<const MAST::LevelSetParameter&>(p);
        
        return _filter.if_elem_in_domain_of_influence(e, *p_ls.level_set_node());
    }
    
private:
    const MAST::FilterBase& _filter;
};


class TopologyOptimizationLevelSet2D:
public MAST::FunctionEvaluation {
    
protected:
    
    bool                                      _initialized;
    MAST::Examples::GetPotWrapper&            _input;
    
    Real                                      _length;
    Real                                      _height;
    Real                                      _obj_scaling;
    Real                                      _stress_penalty;
    Real                                      _perimeter_penalty;
    Real                                      _stress_lim;
    Real                                      _p_val, _vm_rho;
    Real                                      _ref_eig_val;
    unsigned int                              _n_eig_vals;
    
    libMesh::UnstructuredMesh*                _mesh;
    libMesh::UnstructuredMesh*                _level_set_mesh;
    
    libMesh::EquationSystems*                 _eq_sys;
    libMesh::EquationSystems*                 _level_set_eq_sys;
    
    MAST::NonlinearSystem*                    _sys;
    MAST::NonlinearSystem*                    _level_set_sys;
    MAST::NonlinearSystem*                    _level_set_sys_on_str_mesh;
    MAST::NonlinearSystem*                    _indicator_sys;
    
    MAST::StructuralSystemInitialization*     _sys_init;
    MAST::LevelSetSystemInitialization*       _level_set_sys_init_on_str_mesh;
    MAST::LevelSetSystemInitialization*       _level_set_sys_init;
    MAST::HeatConductionSystemInitialization* _indicator_sys_init;
    
    MAST::PhysicsDisciplineBase*              _discipline;
    MAST::PhysicsDisciplineBase*              _indicator_discipline;
    MAST::LevelSetDiscipline*                 _level_set_discipline;
    
    MAST::FilterBase*                         _filter;
    
    MAST::MaterialPropertyCardBase*           _m_card;
    MAST::ElementPropertyCardBase*            _p_card;
    
    PhiMeshFunction*                          _level_set_function;
    MAST::LevelSetBoundaryVelocity*           _level_set_vel;
    libMesh::ExodusII_IO*                     _output;
    
    libMesh::FEType                           _fetype;
    libMesh::FEType                           _level_set_fetype;
    
    std::vector<MAST::Parameter*>             _params_for_sensitivity;
    std::map<std::string, MAST::Parameter*>   _parameters;
    std::set<MAST::FunctionBase*>             _field_functions;
    std::set<MAST::BoundaryConditionBase*>    _boundary_conditions;
    std::set<unsigned int>                    _dv_dof_ids;
    
    std::vector<std::pair<unsigned int, MAST::Parameter*>>  _dv_params;

public:
    
    //  \section  ex_5_init_mesh Mesh Generation
    //  This creates the mesh for the specified problem type.
    //
    void _init_mesh() {
        
        // The mesh is created using classes written in MAST. The particular
        // mesh to be used can be selected using the input parameter
        // ` mesh=val `, where `val` can be one of the following:
        //   - `inplane` inplane structure with load on top and left and right boundaries constrained
        //   - `bracket` L-bracket
        //
        std::string
        s  = _input("mesh", "type of mesh to be analyzed {inplane, bracket}", "inplane");
        
        if (s == "inplane" || s == "truss")
            _init_mesh_inplane();
        else if (s == "bracket")
            _init_mesh_bracket();
        else if (s == "eye_bar")
            _init_mesh_eye_bar();
        else
            libmesh_error();
    }

    //
    //  \subsection ex_5_inplane_mesh Inplane problem
    //
    void _init_mesh_inplane()  {
        
        _mesh = new libMesh::SerialMesh(this->comm());
        
        //
        // identify the element type from the input file or from the order
        // of the element
        //
        unsigned int
        nx_divs = _input("nx_divs", "number of elements along x-axis", 20),
        ny_divs = _input("ny_divs", "number of elements along y-axis", 20);
        
        _length = _input("length", "length of domain along x-axis", 0.3),
        _height = _input("height", "length of domain along y-axis", 0.3);
        
        std::string
        t = _input("elem_type", "type of geometric element in the mesh", "quad4");
        
        libMesh::ElemType
        e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
        
        //
        // if high order FE is used, libMesh requires atleast a second order
        // geometric element.
        //
        if (_fetype.order > 1 && e_type == libMesh::QUAD4)
            e_type = libMesh::QUAD9;
        else if (_fetype.order > 1 && e_type == libMesh::TRI3)
            e_type = libMesh::TRI6;
        
        //
        // initialize the mesh with one element
        //
        libMesh::MeshTools::Generation::build_square(*_mesh,
                                                     nx_divs, ny_divs,
                                                     0, _length,
                                                     0, _height,
                                                     e_type);
        
        //
        // mesh on which the level-set function is defined
        //
        _level_set_mesh = new libMesh::SerialMesh(this->comm());
        
        nx_divs = _input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10);
        ny_divs = _input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
        e_type  = libMesh::QUAD4;
        
        // initialize the mesh with one element
        libMesh::MeshTools::Generation::build_square(*_level_set_mesh,
                                                     nx_divs, ny_divs,
                                                     0, _length,
                                                     0, _height,
                                                     e_type);
    }
    
    //
    //  \subsection ex_5_bracket_mesh Bracket
    //
    void _init_mesh_bracket() {

        {
            unsigned int
            nx_divs = _input("nx_divs", "number of elements along x-axis", 20),
            ny_divs = _input("ny_divs", "number of elements along y-axis", 20);
            
            if (nx_divs%10 != 0 || ny_divs%10 != 0) libmesh_error();
        }
        
        {
            unsigned int
            nx_divs = _input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
            ny_divs = _input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
            
            if (nx_divs%10 != 0 || ny_divs%10 != 0) libmesh_error();
            
        }

        _init_mesh_inplane();
        _delete_elems_from_bracket_mesh(*_mesh);
        _delete_elems_from_bracket_mesh(*_level_set_mesh);
    }

    void _delete_elems_from_bracket_mesh(libMesh::MeshBase &mesh) {
        
        Real
        tol     = 1.e-12,
        x       = -1.,
        y       = -1.,
        length  = _input("length", "length of domain along x-axis", 0.3),
        width   = _input( "height", "length of domain along y-axis", 0.3),
        l_frac  = 0.4,
        w_frac  = 0.4,
        x_lim   = length * l_frac,
        y_lim   =  width * (1.-w_frac);
        
        //
        // now, remove elements that are outside of the L-bracket domain
        //
        libMesh::MeshBase::element_iterator
        e_it   = mesh.elements_begin(),
        e_end  = mesh.elements_end();
        
        for ( ; e_it!=e_end; e_it++) {
            
            libMesh::Elem* elem = *e_it;
            x = length;
            y = 0.;
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                const libMesh::Node& n = elem->node_ref(i);
                if (x > n(0)) x = n(0);
                if (y < n(1)) y = n(1);
            }
          
            //
            // delete element if the lowest x,y locations are outside of the bracket
            // domain
            //
            if (x >= x_lim && y<= y_lim)
                mesh.delete_elem(elem);
        }
        
        mesh.prepare_for_use();
        
        //
        // add the two additional boundaries to the boundary info so that
        // we can apply loads on them
        //
        bool
        facing_right = false,
        facing_down  = false;
        
        e_it   = mesh.elements_begin();
        e_end  = mesh.elements_end();
        
        for ( ; e_it != e_end; e_it++) {
            
            libMesh::Elem* elem = *e_it;
            
            if (!elem->on_boundary()) continue;
            
            for (unsigned int i=0; i<elem->n_sides(); i++) {
                
                if (elem->neighbor_ptr(i)) continue;
                
                std::unique_ptr<libMesh::Elem> s(elem->side_ptr(i).release());
                
                const libMesh::Point p = s->centroid();
                
                facing_right = true;
                facing_down  = true;
                for (unsigned int j=0; j<s->n_nodes(); j++) {
                    const libMesh::Node& n = s->node_ref(j);
                    
                    if (n(0) < x_lim ||  n(1) > y_lim) {
                        facing_right = false;
                        facing_down  = false;
                    }
                    else if (std::fabs(n(0) - p(0)) > tol)
                        facing_right = false;
                    else if (std::fabs(n(1) - p(1)) > tol)
                        facing_down = false;
                }
                
                if (facing_right) mesh.boundary_info->add_side(elem, i, 4);
                if (facing_down) mesh.boundary_info->add_side(elem, i, 5);
            }
        }
        
        mesh.boundary_info->sideset_name(4) = "facing_right";
        mesh.boundary_info->sideset_name(5) = "facing_down";
    }

    
    //
    //  \subsection ex_5_eyebar_mesh Eyebar
    //
    void _init_mesh_eye_bar() {
        
        _mesh = new libMesh::SerialMesh(this->comm());

        //
        // identify the element type from the input file or from the order
        // of the element
        //
        unsigned int
        n_radial_divs  = _input("n_radial_divs", "number of elements along radial direction", 20),
        n_quarter_divs = _input("n_quarter_divs", "number of elements along height", 20);
        
        Real
        radius   = 1.5,
        h_ratio  = _input("h_ratio", "ratio of radial element size at cylinder and at edge", 2);
        _height  = 8.;
        _length  = _height*2;
        
        std::string
        t = _input("elem_type", "type of geometric element in the mesh", "quad4");
        
        libMesh::ElemType
        e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
        
        //
        // if high order FE is used, libMesh requires atleast a second order
        // geometric element.
        //
        if (_fetype.order > 1 && e_type == libMesh::QUAD4)
            e_type = libMesh::QUAD9;
        else if (_fetype.order > 1 && e_type == libMesh::TRI3)
            e_type = libMesh::TRI6;
        
        MAST::Examples::CylinderMesh2D cylinder;
        cylinder.mesh(radius, _height/2.,
                      n_radial_divs, n_quarter_divs, h_ratio,
                      *_mesh, e_type,
                      true, _height, n_quarter_divs*2);
        
        //
        // add the boundary ids for Dirichlet conditions
        //
        libMesh::MeshBase::const_element_iterator
        e_it   = _mesh->elements_begin(),
        e_end  = _mesh->elements_end();
        
        Real
        tol  = radius * 1.e-8;
        
        for (; e_it != e_end; e_it++) {
            
            libMesh::Elem* elem = *e_it;
            
            std::unique_ptr<libMesh::Elem> edge(elem->side_ptr(1));
            libMesh::Point p = edge->centroid();
            
            if (std::fabs(p(0)-_height*1.5) < tol &&
                std::fabs(p(1)) <= 1.) // on the right edge
                _mesh->boundary_info->add_side(elem, 1, 0);
            
            // check for the circumference of the circle where load will be
            // applied
            edge.reset(elem->side_ptr(3).release());
            p = edge->centroid();
            
            if ((std::fabs(p.norm()-radius) < 1.e-2) &&
                p(0) < 0.) // left semi-circle
                _mesh->boundary_info->add_side(elem, 3, 5);
        }
        
        _mesh->boundary_info->sideset_name(0) = "dirichlet";
        _mesh->boundary_info->sideset_name(5) = "load";

        // mesh on which the level-set function is defined
        _level_set_mesh = new libMesh::SerialMesh(this->comm());
        
        n_radial_divs  = _input("level_set_n_radial_divs", "number of elements along radial direction", 10),
        n_quarter_divs = _input("level_set_n_quarter_divs", "number of elements along height", 10);
        e_type  = libMesh::QUAD4;
        
        //
        // initialize the mesh with one element
        //
        cylinder.mesh(radius, _height/2,
                      n_radial_divs, n_quarter_divs, h_ratio,
                      *_level_set_mesh, e_type,
                      true, _height, n_quarter_divs*2);
    }

    //
    //  \section  ex_5_system_discipline  System and Discipline
    //
    void _init_system_and_discipline() {
        
        //
        // make sure that the mesh has been initialized
        //
        libmesh_assert(_mesh);
        
        //
        // create the equation system
        //
        _eq_sys    = new  libMesh::EquationSystems(*_mesh);
        
        //
        // create the libmesh system and set the preferences for structural
        // eigenvalue problems
        //
        _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
        _sys->set_eigenproblem_type(libMesh::GHEP);
        
        //
        // initialize the system to the right set of variables
        //
        _sys_init       = new MAST::StructuralSystemInitialization(*_sys,
                                                                   _sys->name(),
                                                                   _fetype);
        _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
        
        //
        // Initialize the system for level set function.
        // A level set function is defined on a coarser mesh than the structural
        // mesh.
        // A level set function is assumed to be a first-order Lagrange finite element
        //
        _level_set_fetype      = libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE);
        _level_set_eq_sys      = new libMesh::EquationSystems(*_level_set_mesh);
        _level_set_sys         = &(_level_set_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
        _level_set_sys->extra_quadrature_order = 2;
        _level_set_sys_init    = new MAST::LevelSetSystemInitialization(*_level_set_sys,
                                                                        _level_set_sys->name(),
                                                                        _level_set_fetype);
        _level_set_discipline  = new MAST::LevelSetDiscipline(*_eq_sys);
        
        //
        // A system with level set function is defined on the strucutral mesh
        // for the purpose of plotting.
        //
        _level_set_sys_on_str_mesh      = &(_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
        _level_set_sys_init_on_str_mesh = new MAST::LevelSetSystemInitialization(*_level_set_sys_on_str_mesh,
                                                                                 _level_set_sys->name(),
                                                                                 _level_set_fetype);
        
        //
        //  an indicator function is used to locate unconnected free-floating
        // domains of material. The indicator function solves a heat condution
        // problem. Regions with uniformly zero temperature are marked as
        // unconnected domains.
        //
        _indicator_sys                  = &(_eq_sys->add_system<MAST::NonlinearSystem>("indicator"));
        _indicator_sys_init             = new MAST::HeatConductionSystemInitialization(*_indicator_sys,
                                                                                       _indicator_sys->name(),
                                                                                       _fetype);
        _indicator_discipline           = new MAST::PhysicsDisciplineBase(*_eq_sys);
    }

    
    void _init_eq_sys() {
        
        _eq_sys->init();
        _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
        _sys->set_exchange_A_and_B(true);
        
        _level_set_eq_sys->init();
    }
    

    //
    //   variables added to the mesh
    //
    void _init_fetype() {
        
        // FEType to initialize the system. Get the order and type of element.
        std::string
        order_str   = _input("fe_order", "order of finite element shape basis functions",    "first"),
        family_str  = _input("fe_family",      "family of finite element shape functions", "lagrange");
        
        libMesh::Order
        o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
        libMesh::FEFamily
        fe = libMesh::Utility::string_to_enum<libMesh::FEFamily>(family_str);
        _fetype = libMesh::FEType(o, fe);
    }
    

    
    //
    //  \section  ex_5_dirichlet Dirichlet Constraints
    //
    void _init_dirichlet_conditions() {
        
        std::string
        s  = _input("mesh", "type of mesh to be analyzed {inplane, truss, bracket, eye_bar}", "inplane");
        
        if (s == "inplane")
            _init_dirichlet_conditions_inplane();
        else if (s == "truss")
            _init_dirichlet_conditions_truss();
        else if (s == "bracket")
            _init_dirichlet_conditions_bracket();
        else if (s == "eye_bar")
            _init_dirichlet_conditions_eye_bar();
        else
            libmesh_error();
    }
    
    //
    //  \subsection ex_5_inplane_dirichlet Inplane
    //
    void _init_dirichlet_conditions_inplane() {
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for structural system
        ///////////////////////////////////////////////////////////////////////
        MAST::DirichletBoundaryCondition
        *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(1, _sys_init->vars());
        _discipline->add_dirichlet_bc(1,  *dirichlet);

        dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(3, _sys_init->vars());
        _discipline->add_dirichlet_bc(3,  *dirichlet);

        _discipline->init_system_dirichlet_bc(*_sys);
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for indicator system
        ///////////////////////////////////////////////////////////////////////
        dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(1, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(1,  *dirichlet);
        _boundary_conditions.insert(dirichlet);
        
        dirichlet   = new MAST::DirichletBoundaryCondition;   // left boundary
        dirichlet->init(3, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(3,  *dirichlet);
        _boundary_conditions.insert(dirichlet);

        _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
    }
    
    //
    //  \subsection ex_5_truss_dirichlet Truss
    //
    void _init_dirichlet_conditions_truss() {
        
        Real
        dirichlet_length_fraction = _input("truss_dirichlet_length_fraction", "length fraction of the truss boundary where dirichlet condition is applied", 0.05);
        
        // identify the boundaries for dirichlet condition
        libMesh::MeshBase::const_element_iterator
        e_it   = _mesh->elements_begin(),
        e_end  = _mesh->elements_end();
        
        for ( ; e_it != e_end; e_it++) {
            
            const libMesh::Elem* e = *e_it;
            
            if ((*e->node_ptr(0))(1) < 1.e-8 &&
                e->centroid()(0) <= _length*dirichlet_length_fraction)
                _mesh->boundary_info->add_side(e, 0, 6);
            else if ((*e->node_ptr(1))(1) < 1.e-8 &&
                     e->centroid()(0) >= _length*(1.-dirichlet_length_fraction))
                _mesh->boundary_info->add_side(e, 0, 7);
            
            if ((*e->node_ptr(0))(0) < 1.e-8 &&
                (*e->node_ptr(0))(1) < 1.e-8 &&
                e->centroid()(0) <= _length*dirichlet_length_fraction)
                _mesh->boundary_info->add_side(e, 0, 8);
        }
        
        _mesh->boundary_info->sideset_name(6) = "left_dirichlet";
        _mesh->boundary_info->sideset_name(7) = "right_dirichlet";
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for structural system
        ///////////////////////////////////////////////////////////////////////
        std::vector<unsigned int> vars = {1, 2, 3, 4, 5};
        MAST::DirichletBoundaryCondition
        *dirichlet  = new MAST::DirichletBoundaryCondition;   // left support
        dirichlet->init(6, vars);
        _discipline->add_dirichlet_bc(6,  *dirichlet);
        
        dirichlet  = new MAST::DirichletBoundaryCondition;   // right support
        dirichlet->init(7, vars);
        _discipline->add_dirichlet_bc(7,  *dirichlet);

        vars = {0};
        dirichlet  = new MAST::DirichletBoundaryCondition;   // left support
        dirichlet->init(8, vars);
        _discipline->add_dirichlet_bc(8,  *dirichlet);

        _discipline->init_system_dirichlet_bc(*_sys);
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for indicator system
        ///////////////////////////////////////////////////////////////////////
        dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(6, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(6,  *dirichlet);
        _boundary_conditions.insert(dirichlet);
        
        dirichlet   = new MAST::DirichletBoundaryCondition;   // left boundary
        dirichlet->init(7, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(7,  *dirichlet);
        _boundary_conditions.insert(dirichlet);
        
        _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
    }

    
    //
    //  \subsection ex_5_bracket_dirichlet Bracket
    //
    void _init_dirichlet_conditions_bracket() {
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for structural system
        ///////////////////////////////////////////////////////////////////////
        MAST::DirichletBoundaryCondition
        *dirichlet  = new MAST::DirichletBoundaryCondition;   // bottom boundary
        dirichlet->init(0, _sys_init->vars());
        _discipline->add_dirichlet_bc(0,  *dirichlet);
        
        _discipline->init_system_dirichlet_bc(*_sys);

        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for indicator system
        ///////////////////////////////////////////////////////////////////////
        dirichlet  = new MAST::DirichletBoundaryCondition;   // bottom boundary
        dirichlet->init(0, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(0,  *dirichlet);
        _boundary_conditions.insert(dirichlet);
        
        _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
    }
    
    
    //
    //  \subsection ex_5_eyebar_dirichlet Eyebar
    //
    void _init_dirichlet_conditions_eye_bar() {
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for structural system
        ///////////////////////////////////////////////////////////////////////
        MAST::DirichletBoundaryCondition
        *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(0, _sys_init->vars());
        _discipline->add_dirichlet_bc(0,  *dirichlet);
        
        _discipline->init_system_dirichlet_bc(*_sys);
        
        ///////////////////////////////////////////////////////////////////////
        // initialize Dirichlet conditions for indicator system
        ///////////////////////////////////////////////////////////////////////
        dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
        dirichlet->init(0, _indicator_sys_init->vars());
        _indicator_discipline->add_dirichlet_bc(0,  *dirichlet);
        _boundary_conditions.insert(dirichlet);
        
        _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
    }

    

    //
    //  \section  ex_5_loading Loading
    //
    //
    void _init_loads() {
        
        std::string
        s  = _input("mesh", "type of mesh to be analyzed {inplane, truss, bracket, eye_bar}", "inplane");
        
        if (s == "inplane" || s == "truss")
            _init_loads_inplane();
        else if (s == "bracket")
            _init_loads_bracket();
        else if (s == "eye_bar")
            _init_loads_eye_bar();
        else
            libmesh_error();
    }
    
    
    //  \subsection ex_5_inplane_loading Inplane
    //
    class FluxLoad:
    public MAST::FieldFunction<Real> {
    public:
        FluxLoad(const std::string& nm, Real p, Real l1, Real fraction):
        MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _frac(fraction) { }
        virtual ~FluxLoad() {}
        virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            if (fabs(p(0)-_l1*0.5) <= 0.5*_frac*_l1) v = _p;
            else v = 0.;
        }
        virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
            v = 0.;
        }
    protected:
        Real _p, _l1, _frac;
    };
    

    void _init_loads_inplane() {
        
        Real
        frac    = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        p_val   =  _input("pressure", "pressure on side of domain",   2.e4);
        
        FluxLoad
        *press_f         = new FluxLoad( "pressure", p_val, _length, frac),
        *flux_f          = new FluxLoad("heat_flux", -2.e6, _length, frac);
        
        // initialize the load
        MAST::BoundaryConditionBase
        *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE),
        *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
        
        p_load->add(*press_f);
        _discipline->add_side_load(2, *p_load);
        
        f_load->add(*flux_f);
        _indicator_discipline->add_side_load(2, *f_load);
        
        _field_functions.insert(press_f);
        _field_functions.insert(flux_f);
    }
    
    
    //
    //  \subsection ex_5_bracket_loading Bracket
    //
    class BracketLoad:
    public MAST::FieldFunction<Real> {
    public:
        BracketLoad(const std::string& nm, Real p, Real l1, Real fraction):
        MAST::FieldFunction<Real>(nm), _p(p), _l1(l1), _frac(fraction) { }
        virtual ~BracketLoad() {}
        virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            if (fabs(p(0) >= _l1*(1.-_frac))) v = _p;
            else v = 0.;
        }
        virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
            v = 0.;
        }
    protected:
        Real _p, _l1, _frac;
    };
    
    
    
    void _init_loads_bracket() {
        
        Real
        length  = _input("length", "length of domain along x-axis", 0.3),
        frac    = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.125),
        p_val   = _input("pressure", "pressure on side of domain",   5.e7);
        
        BracketLoad
        *press_f         = new BracketLoad( "pressure", p_val, length, frac),
        *flux_f          = new BracketLoad("heat_flux", -2.e6, length, frac);
        
        //
        // initialize the load
        //
        MAST::BoundaryConditionBase
        *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE),
        *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
        
        p_load->add(*press_f);
        _discipline->add_side_load(5, *p_load);
        
        f_load->add(*flux_f);
        _indicator_discipline->add_side_load(5, *f_load);
        
        _field_functions.insert(press_f);
        _field_functions.insert(flux_f);
    }

    
    //
    //  \subsection ex_5_eyebar_loading Eyebar
    //
    class EyebarLoad:
    public MAST::FieldFunction<Real> {
    public:
        EyebarLoad():
        MAST::FieldFunction<Real>("pressure") { }
        virtual ~EyebarLoad() {}
        virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
            if (p(0) <= 0.) v = (-std::pow(p(1), 2) + std::pow(1.5, 2))*1.e6;
            else v = 0.;
        }
        virtual void derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, Real& v) const {
            v = 0.;
        }
    };

    
    
    void _init_loads_eye_bar() {
        
        EyebarLoad
        *press_f         = new EyebarLoad();
        
        // initialize the load
        MAST::BoundaryConditionBase
        *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
        
        p_load->add(*press_f);
        _discipline->add_side_load(5, *p_load);
        
        _field_functions.insert(press_f);
    }

    
    //
    //   \section  ex_5_properties Properties
    //
    //
    //
    //   \subsection ex_5_material_properties Material Properties
    //

    void _init_material() {
        
        Real
        Eval      = _input("E", "modulus of elasticity", 72.e9),
        rhoval    = _input("rho", "material density", 2700.),
        nu_val    = _input("nu", "Poisson's ratio",  0.33),
        kappa_val = _input("kappa", "shear correction factor",  5./6.),
        kval      = _input("k", "thermal conductivity",  1.e-2),
        cpval     = _input("cp", "thermal capacitance",  864.);
        
        
        MAST::Parameter
        *E         = new MAST::Parameter("E",          Eval),
        *rho       = new MAST::Parameter("rho",      rhoval),
        *nu        = new MAST::Parameter("nu",       nu_val),
        *kappa     = new MAST::Parameter("kappa", kappa_val),
        *k         = new MAST::Parameter("k",          kval),
        *cp        = new MAST::Parameter("cp",        cpval);
        
        MAST::ConstantFieldFunction
        *E_f     = new MAST::ConstantFieldFunction(    "E",      *E),
        *rho_f   = new MAST::ConstantFieldFunction(  "rho",    *rho),
        *nu_f    = new MAST::ConstantFieldFunction(   "nu",     *nu),
        *kappa_f = new MAST::ConstantFieldFunction("kappa",  *kappa),
        *k_f     = new MAST::ConstantFieldFunction( "k_th",      *k),
        *cp_f    = new MAST::ConstantFieldFunction(   "cp",     *cp);
        
        _parameters[    E->name()]     = E;
        _parameters[  rho->name()]     = rho;
        _parameters[   nu->name()]     = nu;
        _parameters[kappa->name()]     = kappa;
        _parameters[    k->name()]     = k;
        _parameters[   cp->name()]     = cp;
        _field_functions.insert(E_f);
        _field_functions.insert(rho_f);
        _field_functions.insert(nu_f);
        _field_functions.insert(kappa_f);
        _field_functions.insert(k_f);
        _field_functions.insert(cp_f);

        _m_card = new MAST::IsotropicMaterialPropertyCard;
        _m_card->add(*E_f);
        _m_card->add(*rho_f);
        _m_card->add(*nu_f);
        _m_card->add(*kappa_f);
        _m_card->add(*k_f);
        _m_card->add(*cp_f);
    }

    
    //
    //   \subsection ex_5_section_properties Section Properties
    //

    void _init_section_property(){
        
        
        
        Real
        th_v      =  _input("th", "thickness of 2D element",  0.001);
        
        MAST::Parameter
        *th       = new MAST::Parameter("th", th_v),
        *zero     = new MAST::Parameter("zero", 0.);
        
        MAST::ConstantFieldFunction
        *th_f     = new MAST::ConstantFieldFunction("h",       *th),
        *hoff_f   = new MAST::ConstantFieldFunction("off",   *zero);
        
        
        _parameters[th->name()]    = th;
        _parameters[zero->name()]  = zero;
        _field_functions.insert(th_f);
        _field_functions.insert(hoff_f);
        
        MAST::Solid2DSectionElementPropertyCard
        *p_card   = new MAST::Solid2DSectionElementPropertyCard;
        
        _p_card   = p_card;
        
        // set nonlinear strain if requested
        bool
        nonlinear = _input("if_nonlinear", "flag to turn on/off nonlinear strain", false);
        if (nonlinear) p_card->set_strain(MAST::NONLINEAR_STRAIN);
        
        p_card->add(*th_f);
        p_card->add(*hoff_f);
        p_card->set_material(*_m_card);
        _discipline->set_property_for_subdomain(0, *p_card);
        _indicator_discipline->set_property_for_subdomain(0, *p_card);
    }
    

    //
    //   \section  ex_5_initial_solution Initial Level Set
    //
    //
    class Phi:
    public MAST::FieldFunction<RealVectorX> {
        
    public:
        Phi(Real x0,
            Real y0,
            Real l1,
            Real l2,
            Real nx_mesh,
            Real ny_mesh,
            Real nx_holes,
            Real ny_holes):
        MAST::FieldFunction<RealVectorX>("Phi"),
        _x0  (x0),
        _y0  (y0),
        _l1  (l1),
        _l2  (l2),
        _nx_mesh  (nx_mesh),
        _ny_mesh  (ny_mesh),
        _nx_holes (nx_holes),
        _ny_holes (ny_holes),
        _pi  (acos(-1.)) {

            Real
            dx = _l1/(1.*_nx_holes);
            
            for (unsigned int i=0; i<_nx_holes; i++)
                _x_axis_hole_locations.insert(_x0+(i+.5)*dx);
            
            //
            // now, along the y-axis
            //
            dx = _l2/(1.*_ny_holes);
            for (unsigned int i=0; i<_ny_holes; i++)
                _y_axis_hole_locations.insert(_y0+(i+0.5)*dx);
        }
        virtual ~Phi() {}
        virtual void operator()(const libMesh::Point& p,
                                const Real t,
                                RealVectorX& v) const {
            
            libmesh_assert_less_equal(t, 1);
            libmesh_assert_equal_to(v.size(), 1);
            
            //
            // the libMesh solution projection routine for Lagrange elements
            // will query the function value at the nodes. So, we figure
            // out which nodes should have zero values set to them.
            // if there is one hole in any direction, it will be in the
            // center of the domain. If there are more than 1, then two of
            // the holes will be on the boundary and others will fill the
            // interior evenly.
            //
            const Real
            dx_mesh = _l1/(1.*_nx_holes),
            dy_mesh = _l2/(1.*_ny_holes);
            
            std::set<Real>::const_iterator
            x_it_low = _x_axis_hole_locations.lower_bound(p(0)-dx_mesh),
            y_it_low = _y_axis_hole_locations.lower_bound(p(1)-dy_mesh);
            
            unsigned int
            n = 0;
            //
            // see if the x-location needs a hole
            //
            for ( ; x_it_low != _x_axis_hole_locations.end(); x_it_low++) {
                if (std::fabs(*x_it_low - p(0)) <= dx_mesh*0.25) {
                    n++;
                    break;
                }
            }
            
            //
            // now check the y-location
            //
            for ( ; y_it_low != _y_axis_hole_locations.end(); y_it_low++) {
                if (std::fabs(*y_it_low - p(1)) <= dy_mesh*0.25) {
                    n++;
                    break;
                }
            }
            
            if (n == 2)
                v(0) = -1.e0;
            else
                v(0) = 1.e0;
        }
    protected:
        Real
        _x0,
        _y0,
        _l1,
        _l2,
        _nx_mesh,
        _ny_mesh,
        _nx_holes,
        _ny_holes,
        _pi;
        std::set<Real> _x_axis_hole_locations;
        std::set<Real> _y_axis_hole_locations;
    };
    
    
    void initialize_solution() {
        
        //
        // initialize solution of the level set problem
        //
        unsigned int
        nx_h    = _input("initial_level_set_n_holes_in_x",
                            "number of holes along x-direction for initial level-set field", 6),
        ny_h    = _input("initial_level_set_n_holes_in_y",
                            "number of holes along y-direction for initial level-set field", 6),
        nx_m    = _input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
        ny_m    = _input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);

        std::string
        s  = _input("mesh", "type of mesh to be analyzed {inplane, truss, bracket, eyebar}", "inplane");

        std::unique_ptr<Phi> phi;
        
        if (s == "inplane" || s == "truss" || s == "bracket")
            phi.reset(new Phi(0., 0., _length, _height, nx_m, ny_m, nx_h, ny_h));
        else if (s == "eye_bar")
            phi.reset(new Phi(-0.5*_height, -0.5*_height, _length, _height, nx_m, ny_m, nx_h, ny_h));
        else
            libmesh_error();
        
        _level_set_sys_init->initialize_solution(*phi);
    }
    
    
    void _init_phi_dvs() {
        
        std::string
        s  = _input("mesh", "type of mesh to be analyzed {inplane, truss, bracket, eye_bar}", "inplane");
        
        if (s == "inplane" || s == "truss")
            _init_phi_dvs_inplane();
        else if (s == "bracket")
            _init_phi_dvs_bracket();
        else if (s == "eye_bar")
            _init_phi_dvs_eye_bar();
        else
            libmesh_error();
        
        Real
        filter_radius          = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        _filter                = new MAST::FilterBase(*_level_set_sys, filter_radius, _dv_dof_ids);
        libMesh::NumericVector<Real>& vec = _level_set_sys->add_vector("base_values");
        vec = *_level_set_sys->solution;
        vec.close();
    }
    
    //
    //  \subsection ex_5_inplane_initial_level_set Inplane
    //
    void _init_phi_dvs_inplane() {
        
        //
        // this assumes that level set is defined using lagrange shape functions
        //
        libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
        
        Real
        frac          = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_level_set_mesh->is_replicated());
        
        std::vector<Real> local_phi(_level_set_sys->solution->size());
        _level_set_sys->solution->localize(local_phi);
        
        // iterate over all the node values
        libMesh::MeshBase::const_node_iterator
        it  = _level_set_mesh->nodes_begin(),
        end = _level_set_mesh->nodes_end();
        
        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_level_set_mesh->n_nodes());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(0, 0, 0);
            
            // only if node is not on the upper edge
            if ((n(1)+filter_radius >= _height) &&
                (n(0)-filter_radius <= _length*.5*(1.+frac))   &&
                (n(0)+filter_radius >= _length*.5*(1.-frac))) {
                
                // set value at the material points to a small positive number
                if (dof_id >= _level_set_sys->solution->first_local_index() &&
                    dof_id <  _level_set_sys->solution->last_local_index())
                    _level_set_sys->solution->set(dof_id, 1.e0);
            }
            else {

                std::ostringstream oss;
                oss << "dv_" << _n_vars;
                val  = local_phi[dof_id];
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);

                _n_vars++;
            }
        }
        
        _level_set_sys->solution->close();
    }
    
    //
    //  \subsection ex_5_truss_initial_level_set Truss
    //
    void _init_phi_dvs_truss() {
        
        //
        // this assumes that level set is defined using lagrange shape functions
        //
        libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
        
        Real
        frac          = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_level_set_mesh->is_replicated());
        
        std::vector<Real> local_phi(_level_set_sys->solution->size());
        _level_set_sys->solution->localize(local_phi);
        
        // iterate over all the node values
        libMesh::MeshBase::const_node_iterator
        it  = _level_set_mesh->nodes_begin(),
        end = _level_set_mesh->nodes_end();
        
        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_level_set_mesh->n_nodes());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(0, 0, 0);
            
            // only if node is not on the upper edge
            if ((n(1)-filter_radius <= 0.) &&
                (n(0)-filter_radius <= _length*.5*(1.+frac))   &&
                (n(0)+filter_radius >= _length*.5*(1.-frac))) {
                
                // set value at the material points to a small positive number
                if (dof_id >= _level_set_sys->solution->first_local_index() &&
                    dof_id <  _level_set_sys->solution->last_local_index())
                    _level_set_sys->solution->set(dof_id, 1.e0);
            }
            else {
                
                std::ostringstream oss;
                oss << "dv_" << _n_vars;
                val  = local_phi[dof_id];
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);

                _n_vars++;
            }
        }
        
        _level_set_sys->solution->close();
    }
    
    //
    //  \subsection ex_5_bracket_initial_level_set Bracket
    //
    void _init_phi_dvs_bracket() {
        
        libmesh_assert(_initialized);
        //
        // this assumes that level set is defined using lagrange shape functions
        //
        libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
        
        Real
        tol           = 1.e-12,
        length        = _input("length", "length of domain along x-axis", 0.3),
        height        = _input("height", "length of domain along y-axis", 0.3),
        l_frac        = 0.4,//_input("length_fraction", "fraction of length along x-axis that is in the bracket", 0.4),
        h_frac        = 0.4,//_input( "height_fraction", "fraction of length along y-axis that is in the bracket", 0.4),
        x_lim         = length * l_frac,
        y_lim         =  height * (1.-h_frac),
        frac          = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.125),
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_level_set_mesh->is_replicated());
        
        std::vector<Real> local_phi(_level_set_sys->solution->size());
        _level_set_sys->solution->localize(local_phi);
        
        //
        // iterate over all the node values
        //
        libMesh::MeshBase::const_node_iterator
        it  = _level_set_mesh->nodes_begin(),
        end = _level_set_mesh->nodes_end();
        
        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_level_set_mesh->n_nodes());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(0, 0, 0);
            
            if ((n(1)-filter_radius) <= y_lim && (n(0)+filter_radius) >= length*(1.-frac)) {
          
                //
                // set value at the constrained points to a small positive number
                // material here
                //
                if (dof_id >= _level_set_sys->solution->first_local_index() &&
                    dof_id <  _level_set_sys->solution->last_local_index())
                    _level_set_sys->solution->set(dof_id, 1.e0);
            }
            else {
                
                std::ostringstream oss;
                oss << "dv_" << _n_vars;
                val = local_phi[dof_id];
                
                //
                // on the boundary, set everything to be zero, so that there
                // is always a boundary there that the optimizer can move
                //
                if (n(0) < tol                     ||  // left boundary
                    std::fabs(n(0) - length) < tol ||  // right boundary
                    std::fabs(n(1) - height) < tol ||  // top boundary
                    (n(0) >= x_lim && n(1) <= y_lim)) {
                    
                    if (dof_id >= _level_set_sys->solution->first_local_index() &&
                        dof_id <  _level_set_sys->solution->last_local_index())
                        _level_set_sys->solution->set(dof_id, -1.0);
                    val = -1.0;
                }
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);
                
                _n_vars++;
            }
        }
        
        _level_set_sys->solution->close();
    }
    
    
    //
    //  \subsection ex_5_eyebar_initial_level_set Eyebar
    //
    void _init_phi_dvs_eye_bar() {
        
        libmesh_assert(_initialized);
        //
        // this assumes that level set is defined using lagrange shape functions
        //
        libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
        
        Real
        tol           = 1.e-6,
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_level_set_mesh->is_replicated());
        
        std::vector<Real> local_phi(_level_set_sys->solution->size());
        _level_set_sys->solution->localize(local_phi);
        
        //
        // iterate over all the node values
        //
        libMesh::MeshBase::const_node_iterator
        it  = _level_set_mesh->nodes_begin(),
        end = _level_set_mesh->nodes_end();
        
        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_level_set_mesh->n_nodes());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(0, 0, 0);
            
            if (((n.norm() <= 1.5+filter_radius) && n(0) <= 0.) ||  // circle
                (std::fabs(n(0)-_height*1.5) < filter_radius &&  // right edge
                 std::fabs(n(1)) <= 1.+filter_radius)) { // dirichlet constraint
                    
                    //
                    // set value at the constrained points to a small positive number
                    // material here
                    //
                    if (dof_id >= _level_set_sys->solution->first_local_index() &&
                        dof_id <  _level_set_sys->solution->last_local_index())
                        _level_set_sys->solution->set(dof_id, 1.e0);
                }
            else {
                
                std::ostringstream oss;
                oss << "dv_" << _n_vars;
                val = local_phi[dof_id];
                
                //
                // on the boundary, set everything to be zero, so that there
                // is always a boundary there that the optimizer can move
                //
                if (std::fabs(n(0)+_height*0.5) < tol    ||  // left boundary
                    std::fabs(n(1)-_height*0.5) < tol    ||  // top boundary
                    std::fabs(n(1)+_height*0.5) < tol    ||  // bottom boundary
                    std::fabs(n(0)-_height*1.5) < tol) {     // right boundary
                    
                    if (dof_id >= _level_set_sys->solution->first_local_index() &&
                        dof_id <  _level_set_sys->solution->last_local_index())
                        _level_set_sys->solution->set(dof_id, -1.);
                    val = -1.;
                }
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);
                
                _n_vars++;
            }
        }
        
        _level_set_sys->solution->close();
    }

    
    //
    //   \subsection ex_5_design_variable_init   Design Variables
    //
    //   initializes the design variable vector, called by the
    //   optimization interface.
    //
    void init_dvar(std::vector<Real>& x,
                   std::vector<Real>& xmin,
                   std::vector<Real>& xmax) {
        
        x.resize(_n_vars);
        xmin.resize(_n_vars);
        xmax.resize(_n_vars);
        
        std::fill(xmin.begin(), xmin.end(),   -1.e0);
        std::fill(xmax.begin(), xmax.end(),    1.e0);

        //
        // now, check if the user asked to initialize dvs from a previous file
        //
        std::string
        nm    =  _input("restart_optimization_file", "filename with optimization history for restart", "");
        
        if (nm.length()) {
            
            unsigned int
            iter = _input("restart_optimization_iter", "restart iteration number from file", 0);
            this->initialize_dv_from_output_file(nm, iter, x);
        }
        else {
            
            for (unsigned int i=0; i<_n_vars; i++)
                x[i] = (*_dv_params[i].second)();
        }
    }

    //
    //  \section  ex_5_analysis Function Evaluation and Sensitivity
    //
    //
    //   \subsection ex_5_element_error_metric Element Error Metric
    //
    void
    _compute_element_errors(libMesh::ErrorVector& error) {
        
        MAST::LevelSetIntersection intersection;
        
        libMesh::MeshBase::const_element_iterator
        it  = _mesh->active_elements_begin(),
        end = _mesh->active_elements_end();
        
        for ( ; it != end; it++) {
            
            const libMesh::Elem* elem = *it;
            intersection.init( *_level_set_function, *elem, _sys->time,
                              _mesh->max_elem_id(),
                              _mesh->max_node_id());
            if (intersection.if_intersection_through_elem())
                error[elem->id()] = 1.-intersection.get_positive_phi_volume_fraction();
            intersection.clear();
        }
    }
    
    class ElemFlag: public libMesh::MeshRefinement::ElementFlagging {
    public:
        ElemFlag(libMesh::MeshBase& mesh, MAST::FieldFunction<Real>& phi, unsigned int max_h):
        _mesh(mesh), _phi(phi), _max_h(max_h) {}
        virtual ~ElemFlag() {}
        virtual void flag_elements () {
            
            MAST::LevelSetIntersection intersection;
            
            libMesh::MeshBase::element_iterator
            it  = _mesh.active_elements_begin(),
            end = _mesh.active_elements_end();
            
            for ( ; it != end; it++) {
                
                libMesh::Elem* elem = *it;
                intersection.init( _phi, *elem, 0.,
                                  _mesh.max_elem_id(),
                                  _mesh.max_node_id());
                if (intersection.if_intersection_through_elem()) {
                    
                    Real vol_frac = intersection.get_positive_phi_volume_fraction();
                    if (vol_frac < 0.5 && elem->level() < _max_h)
                        elem->set_refinement_flag(libMesh::Elem::REFINE);
                    else if (vol_frac > 0.90)
                        elem->set_refinement_flag(libMesh::Elem::COARSEN);
                }
                else
                    elem->set_refinement_flag(libMesh::Elem::COARSEN);
                intersection.clear();
            }
        }
        
    protected:
        libMesh::MeshBase& _mesh;
        MAST::FieldFunction<Real>& _phi;
        unsigned int _max_h;
    };
    

    //
    //  \subsection ex_5_function_evaluation Function Evaluation
    //
    void evaluate(const std::vector<Real>& dvars,
                  Real& obj,
                  bool eval_obj_grad,
                  std::vector<Real>& obj_grad,
                  std::vector<Real>& fvals,
                  std::vector<bool>& eval_grads,
                  std::vector<Real>& grads) {
        
        libMesh::out << "New Evaluation" << std::endl;
        
        // copy DVs to level set function
        libMesh::NumericVector<Real>
        &base_phi = _level_set_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, dvars[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_level_set_sys->solution);
        _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
        _sys->solution->zero();
        
        //*********************************************************************
        // DO NOT zero out the gradient vector, since GCMMA needs it for the  *
        // subproblem solution                                                *
        //*********************************************************************
        MAST::LevelSetNonlinearImplicitAssembly                  nonlinear_assembly(true);
        MAST::LevelSetNonlinearImplicitAssembly                  level_set_assembly(false);
        MAST::LevelSetEigenproblemAssembly                       eigen_assembly;
        MAST::LevelSetStressAssembly                             stress_assembly;
        MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
        MAST::HeatConductionNonlinearAssemblyElemOperations      conduction_elem_ops;
        MAST::StructuralModalEigenproblemAssemblyElemOperations  modal_elem_ops;
        
        //
        // reinitialize the dof constraints before solution of the linear system
        // FIXME: we should be able to clear the constraint object from the
        // system before it goes out of scope, but libMesh::System does not
        // have a clear method. So, we are going to leave it as is, hoping
        // that libMesh::System will not attempt to use it (most likely, we
        // shoudl be ok).
        //
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        SNESConvergedReason r;
        /*{
            libMesh::out << "Indicator Function" << std::endl;
            nonlinear_assembly.set_discipline_and_system(*_indicator_discipline, *_indicator_sys_init);
            conduction_elem_ops.set_discipline_and_system(*_indicator_discipline, *_indicator_sys_init);
            nonlinear_assembly.set_level_set_function(*_level_set_function);
            
            MAST::LevelSetConstrainDofs constrain(*_indicator_sys_init, *_level_set_function);
            constrain.constrain_all_negative_indices(true);
            _indicator_sys->attach_constraint_object(constrain);
            _indicator_sys->reinit_constraints();
            _indicator_sys->solve(conduction_elem_ops, nonlinear_assembly);
            r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
            (*_indicator_sys->nonlinear_solver).get_converged_reason();
            nonlinear_assembly.clear_level_set_function();
            nonlinear_assembly.clear_discipline_and_system();
        }
        // if the solver diverged due to linear solve, then there is a problem with
        // this geometry and we need to return with a high value set for the
        // constraints
        if (r == SNES_DIVERGED_LINEAR_SOLVE) {
            
            obj = 1.e10;
            for (unsigned int i=0; i<_n_ineq; i++)
                fvals[i] = 1.e10;
            return;
        }
        
        
        /////////////////////////////////////////////////////////////////////
        // now, use the indicator function to constrain dofs in the structural
        // system
        /////////////////////////////////////////////////////////////////////
        MAST::MeshFieldFunction indicator(*_indicator_sys_init, "indicator");
        indicator.init(*_indicator_sys->solution);
        MAST::IndicatorFunctionConstrainDofs constrain(*_sys_init, *_level_set_function, indicator);
        MAST::LevelSetConstrainDofs constrain(*_sys_init, *_level_set_function);
        _sys->attach_constraint_object(constrain);
        _sys->reinit_constraints();
        _sys->initialize_condensed_dofs(*_discipline);*/
        
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        nonlinear_assembly.set_level_set_function(*_level_set_function, *_filter);
        nonlinear_assembly.set_level_set_velocity_function(*_level_set_vel);
        //nonlinear_assembly.set_indicator_function(indicator);
        eigen_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        eigen_assembly.set_level_set_function(*_level_set_function);
        eigen_assembly.set_level_set_velocity_function(*_level_set_vel);
        stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        stress_assembly.init(*_level_set_function, nonlinear_assembly.if_use_dof_handler()?&nonlinear_assembly.get_dof_handler():nullptr);
        level_set_assembly.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
        level_set_assembly.set_level_set_function(*_level_set_function, *_filter);
        level_set_assembly.set_level_set_velocity_function(*_level_set_vel);
        nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        modal_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        //nonlinear_assembly.plot_sub_elems(true, false, true);
        
        
        libMesh::MeshRefinement refine(*_mesh);
        
        libMesh::out << "before refinement" << std::endl;
        _mesh->print_info();

        bool
        continue_refining    = true;
        Real
        threshold            = 0.05;
        unsigned int
        n_refinements        = 0,
        max_refinements      = _input("max_refinements","maximum refinements", 3);
        
        while (n_refinements < max_refinements && continue_refining) {
            
            // The ErrorVector is a particular StatisticsVector
            // for computing error information on a finite element mesh.
            libMesh::ErrorVector error(_mesh->max_elem_id(), _mesh);
            _compute_element_errors(error);
            libMesh::out
            << "After refinement: " << n_refinements << std::endl
            << "max error:    " << error.maximum()
            << ",  mean error: " << error.mean() << std::endl;

            if (error.maximum() > threshold) {
                
                ElemFlag flag(*_mesh, *_level_set_function, max_refinements);
                refine.max_h_level()      = max_refinements;
                refine.refine_fraction()  = 1.;
                refine.coarsen_fraction() = 0.5;
                refine.flag_elements_by (flag);
                if (refine.refine_and_coarsen_elements())
                    _eq_sys->reinit ();

                _mesh->print_info();
                
                n_refinements++;
            }
            else
                continue_refining = false;
        }

        
        MAST::LevelSetVolume                            volume(level_set_assembly.get_intersection());
        MAST::LevelSetPerimeter                         perimeter(level_set_assembly.get_intersection());
        MAST::StressStrainOutputBase                    stress;
        MAST::ComplianceOutput                          compliance;
        volume.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
        perimeter.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
        stress.set_discipline_and_system(*_discipline, *_sys_init);
        volume.set_participating_elements_to_all();
        perimeter.set_participating_elements_to_all();
        stress.set_participating_elements_to_all();
        stress.set_aggregation_coefficients(_p_val, 1., _vm_rho, _stress_lim) ;
        compliance.set_participating_elements_to_all();
        compliance.set_discipline_and_system(*_discipline, *_sys_init);

        //////////////////////////////////////////////////////////////////////
        // evaluate the stress constraint
        //////////////////////////////////////////////////////////////////////
        // tell the thermal jacobian scaling object about the assembly object
        
        libMesh::out << "Static Solve" << std::endl;
        _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
        r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
        (*_sys->nonlinear_solver).get_converged_reason();
        
        // if the solver diverged due to linear solve, then there is a problem with
        // this geometry and we need to return with a high value set for the
        // constraints
        if (r == SNES_DIVERGED_LINEAR_SOLVE ||
            _sys->final_nonlinear_residual() > 1.e-1) {
            
            obj = 1.e10;
            for (unsigned int i=0; i<_n_ineq; i++)
                fvals[i] = 1.e10;
            return;
        }
        
        nonlinear_assembly.calculate_output(*_sys->solution, stress);
        //nonlinear_assembly.calculate_output(*_sys->solution, compliance);
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the objective
        //////////////////////////////////////////////////////////////////////
        level_set_assembly.set_evaluate_output_on_negative_phi(false);
        level_set_assembly.calculate_output(*_level_set_sys->solution, volume);
        level_set_assembly.set_evaluate_output_on_negative_phi(true);
        level_set_assembly.calculate_output(*_level_set_sys->solution, perimeter);
        level_set_assembly.set_evaluate_output_on_negative_phi(false);

        Real
        max_vm = stress.get_maximum_von_mises_stress(),
        vm_agg = stress.output_total(),
        vf     = _input("volume_fraction", "volume fraction", 0.3),
        vol    = volume.output_total(),
        per    = perimeter.output_total(),
        comp   = compliance.output_total();
        
        obj       = _obj_scaling * (vol + _perimeter_penalty * per);
        //_obj_scaling    * (vol+ _perimeter_penalty * per) +
        //_stress_penalty * (vm_agg);///_stress_lim - 1.);
        
        fvals[0]  =  stress.output_total()/_stress_lim - 1.;  // g = sigma/sigma0-1 <= 0
        //fvals[0]  =  stress.output_total()/_length/_height;  // g <= 0 for the smooth ramp function
        //fvals[0]  = vol/_length/_height - vf; // vol/vol0 - a <=
        libMesh::out << "volume: " << vol    << "  perim: "  << per    << std::endl;
        libMesh::out << "max: "    << max_vm << "  constr: " << vm_agg/_stress_lim - 1.
        << std::endl;
        libMesh::out << "compliance: " << comp << std::endl;

        if (_n_eig_vals) {
            
            //////////////////////////////////////////////////////////////////////
            // evaluate the eigenvalue constraint
            //////////////////////////////////////////////////////////////////////
            libMesh::out << "Eigen Solve" << std::endl;
            _sys->eigenproblem_solve(modal_elem_ops, eigen_assembly);
            Real eig_imag = 0.;
            //
            // hopefully, the solver found the requested number of eigenvalues.
            // if not, then we will set zero values for the ones it did not.
            //
            unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());
            std::vector<Real> eig(_n_eig_vals, 0.);
            
            // get the converged eigenvalues
            for (unsigned int i=0; i<n_conv; i++)      _sys->get_eigenvalue(0, eig[i], eig_imag);
            //
            //  eig > eig0
            //  -eig < -eig0
            //  -eig/eig0 < -1
            // -eig/eig0 + 1 < 0
            //
            for (unsigned int i=0; i<_n_eig_vals; i++)
                fvals[i+1] = -eig[i]/_ref_eig_val + 1.;
        }
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the objective sensitivities, if requested
        //////////////////////////////////////////////////////////////////////
        if (eval_obj_grad) {
            
            std::vector<Real>
            grad1(obj_grad.size(), 0.);
            
            _evaluate_volume_sensitivity(&volume, &perimeter, level_set_assembly, obj_grad);
            
            /*_evaluate_compliance_sensitivity(compliance,
                                             nonlinear_elem_ops,
                                             nonlinear_assembly,
                                             grad1);
            
            for (unsigned int i=0; i<obj_grad.size(); i++)
                obj_grad[i] += _obj_scaling * grad1[i];*/
        }
        
        //////////////////////////////////////////////////////////////////////
        // check to see if the sensitivity of constraint is requested
        //////////////////////////////////////////////////////////////////////
        bool if_grad_sens = false;
        for (unsigned int i=0; i<eval_grads.size(); i++)
            if_grad_sens = (if_grad_sens || eval_grads[i]);
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the sensitivities for constraints
        //////////////////////////////////////////////////////////////////////
        if (if_grad_sens)
            _evaluate_stress_sensitivity(stress,
                                         nonlinear_elem_ops,
                                         nonlinear_assembly,
                                         modal_elem_ops,
                                         eigen_assembly,
                                         grads);
            
            //_evaluate_volume_sensitivity(&volume, nullptr, level_set_assembly, grads);

        
        //
        // also the stress data for plotting
        //
        stress_assembly.update_stress_strain_data(stress, *_sys->solution);
    }

    //
    //  \subsection ex_5_volume_sensitivity Sensitivity of Material Volume
    //
    void _evaluate_volume_sensitivity(MAST::LevelSetVolume*    volume,
                                      MAST::LevelSetPerimeter* perimeter,
                                      MAST::LevelSetNonlinearImplicitAssembly& assembly,
                                      std::vector<Real>& grad) {
        
        std::fill(grad.begin(), grad.end(), 0.);
        
        //
        // iterate over each DV, create a sensitivity vector and calculate the
        // volume sensitivity explicitly
        //
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi_base(_level_set_sys->solution->zero_clone().release()),
        dphi_filtered(_level_set_sys->solution->zero_clone().release());
        
        ElementParameterDependence dep(*_filter);
        assembly.attach_elem_parameter_dependence_object(dep);
        
        for (unsigned int i=0; i<_n_vars; i++) {
            
            dphi_base->zero();
            dphi_filtered->zero();
            //
            // set the value only if the dof corresponds to a local node
            //
            if (_dv_params[i].first >=  dphi_base->first_local_index() &&
                _dv_params[i].first <   dphi_base->last_local_index())
                dphi_base->set(_dv_params[i].first, 1.);
            dphi_base->close();
            _filter->compute_filtered_values(*dphi_base, *dphi_filtered);
            
            _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi_filtered);

            // if the volume output was specified then compute the sensitivity
            // and add to the grad vector
            if (volume) {
                
                assembly.set_evaluate_output_on_negative_phi(false);
                assembly.calculate_output_direct_sensitivity(*_level_set_sys->solution,
                                                             dphi_filtered.get(),
                                                             *_dv_params[i].second,
                                                             *volume);
                
                grad[i] = _obj_scaling * volume->output_sensitivity_total(*_dv_params[i].second);
                //grad[i] = volume->output_sensitivity_total(*_dv_params[i].second)/_length/_height;
            }
            
            // if the perimeter output was specified then compute the sensitivity
            // and add to the grad vector
            if (perimeter) {
                assembly.set_evaluate_output_on_negative_phi(true);
                assembly.calculate_output_direct_sensitivity(*_level_set_sys->solution,
                                                             dphi_filtered.get(),
                                                             *_dv_params[i].second,
                                                             *perimeter);
                assembly.set_evaluate_output_on_negative_phi(false);
                
                grad[i] += _obj_scaling * _perimeter_penalty *
                perimeter->output_sensitivity_total(*_dv_params[i].second);
            }
        }
        
        assembly.clear_elem_parameter_dependence_object();
    }
    
    
    
    //
    //  \subsection ex_5_stress_sensitivity Sensitivity of Stress and Eigenvalues
    //
    void
    _evaluate_stress_sensitivity
    (MAST::StressStrainOutputBase& stress,
     MAST::AssemblyElemOperations& nonlinear_elem_ops,
     MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
     MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
     MAST::LevelSetEigenproblemAssembly& eigen_assembly,
     std::vector<Real>& grads) {
        
        unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());
        
        _sys->adjoint_solve(nonlinear_elem_ops, stress, nonlinear_assembly, false);
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi_base(_level_set_sys->solution->zero_clone().release()),
        dphi_filtered(_level_set_sys->solution->zero_clone().release());

        ElementParameterDependence dep(*_filter);
        nonlinear_assembly.attach_elem_parameter_dependence_object(dep);

        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<_n_vars; i++) {
            
            dphi_base->zero();
            dphi_filtered->zero();
            //
            // set the value only if the dof corresponds to a local node
            //
            if (_dv_params[i].first >=  dphi_base->first_local_index() &&
                _dv_params[i].first <   dphi_base->last_local_index())
                dphi_base->set(_dv_params[i].first, 1.);
            dphi_base->close();
            _filter->compute_filtered_values(*dphi_base, *dphi_filtered);
            
            //
            // initialize the level set perturbation function to create a velocity
            // field
            _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi_filtered);
            
            //////////////////////////////////////////////////////////////////////
            // stress sensitivity
            //////////////////////////////////////////////////////////////////////
            grads[1*i+0] = 1./_stress_lim*
            nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                    _sys->get_adjoint_solution(),
                                                                    *_dv_params[i].second,
                                                                    nonlinear_elem_ops,
                                                                    stress);
            stress.clear_sensitivity_data();
            
            //////////////////////////////////////////////////////////////////////
            // eigenvalue sensitivity, only if the values were requested
            //////////////////////////////////////////////////////////////////////
            if (_n_eig_vals) {
                
                std::vector<Real> sens;
                _sys->eigenproblem_sensitivity_solve(eigen_elem_ops,
                                                     eigen_assembly,
                                                     *_dv_params[i].second,
                                                     sens);
                for (unsigned int j=0; j<n_conv; j++)
                    grads[_n_ineq*i+j+1] = -sens[j]/_ref_eig_val;
            }
        }
        
        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    
    void
    _evaluate_compliance_sensitivity
    (MAST::ComplianceOutput& compliance,
     MAST::AssemblyElemOperations& nonlinear_elem_ops,
     MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
     std::vector<Real>& grads) {
        
        // Adjoint solution for compliance = - X
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi_base(_level_set_sys->solution->zero_clone().release()),
        dphi_filtered(_level_set_sys->solution->zero_clone().release());
        
        ElementParameterDependence dep(*_filter);
        nonlinear_assembly.attach_elem_parameter_dependence_object(dep);

        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<_n_vars; i++) {
            
            dphi_base->zero();
            dphi_filtered->zero();
            //
            // set the value only if the dof corresponds to a local node
            //
            if (_dv_params[i].first >=  dphi_base->first_local_index() &&
                _dv_params[i].first <   dphi_base->last_local_index())
                dphi_base->set(_dv_params[i].first, 1.);
            dphi_base->close();
            _filter->compute_filtered_values(*dphi_base, *dphi_filtered);
            
            //
            // initialize the level set perturbation function to create a velocity
            // field
            _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi_filtered);
            
            //////////////////////////////////////////////////////////////////////
            // compliance sensitivity
            //////////////////////////////////////////////////////////////////////
            grads[i] = -1. *
            nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                    *_sys->solution,
                                                                    *_dv_params[i].second,
                                                                    nonlinear_elem_ops,
                                                                    compliance);
        }
        
        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    //
    //  \subsection ex_5_design_output  Output of Design Iterate
    //
    void output(unsigned int iter,
                const std::vector<Real>& x,
                Real obj,
                const std::vector<Real>& fval,
                bool if_write_to_optim_file) {
        
        libmesh_assert_equal_to(x.size(), _n_vars);
        
        Real
        sys_time     = _sys->time;
        
        std::string
        output_name  = _input("output_file_root", "prefix of output file names", "output"),
        modes_name   = output_name + "modes.exo";
        
        std::ostringstream oss;
        oss << "output_optim.e-s." << std::setfill('0') << std::setw(5) << iter ;
        
        //
        // copy DVs to level set function
        //
        libMesh::NumericVector<Real>
        &base_phi = _level_set_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, x[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_level_set_sys->solution);
        _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
        _level_set_sys_init_on_str_mesh->initialize_solution(_level_set_function->get_mesh_function());
        
        std::vector<bool> eval_grads(this->n_ineq(), false);
        std::vector<Real> f(this->n_ineq(), 0.), grads;
        this->evaluate(x, obj, false, grads, f, eval_grads, grads);
        
        _sys->time = iter;
        _sys_init->get_stress_sys().time = iter;
        // "1" is the number of time-steps in the file, as opposed to the time-step number.
        libMesh::ExodusII_IO(*_mesh).write_timestep(oss.str(), *_eq_sys, 1, (1.*iter));
        
        if (_n_eig_vals) {
            
            //////////////////////////////////////////////////////////////////////////
            // eigenvalue analysis: write modes to file
            //////////////////////////////////////////////////////////////////////////
            libMesh::ExodusII_IO writer(*_mesh);
            Real eig_r, eig_i;
            for (unsigned int i=0; i<_sys->get_n_converged_eigenvalues(); i++) {
                _sys->get_eigenpair(i, eig_r, eig_i, *_sys->solution);
                writer.write_timestep(modes_name, *_eq_sys, i+1, i);
            }
            _sys->solution->zero();
        }
        
        //
        // set the value of time back to its original value
        //
        _sys->time = sys_time;
        
        //
        // increment the parameter values
        //
        unsigned int
        update_freq = _input("update_freq_optim_params", "number of iterations after which the optimization parameters are updated", 50),
        factor = iter/update_freq ;
        if (factor > 0 && iter%update_freq == 0) {
            
            Real
            p_val           = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0),
            vm_rho          = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0),
            constr_penalty  = _input("constraint_penalty", "constraint penalty in GCMMA",      50.),
            max_penalty     = _input("max_constraint_penalty", "maximum constraint penalty in GCMMA",      1.e7),
            initial_step    = _input("initial_rel_step", "initial relative step length in GCMMA",      0.5),
            min_step        = _input("minimum_rel_step", "minimum relative step length in GCMMA",      0.001);
            
            constr_penalty = std::min(constr_penalty*pow(10, factor), max_penalty);
            initial_step   = std::max(initial_step-0.01*factor, min_step);
            _p_val         = std::min(p_val+2*factor, 10.);
            _vm_rho        = std::min(vm_rho+factor*0.5, 2.);
            libMesh::out
            << "Updated values: c = " << constr_penalty
            << "  step = " << initial_step
            << "  p = " << _p_val
            << "  rho = " << _vm_rho << std::endl;
            
            _optimization_interface->set_real_parameter   ( "constraint_penalty",   constr_penalty);
            _optimization_interface->set_real_parameter   ("initial_rel_step",        initial_step);
        }

        MAST::FunctionEvaluation::output(iter, x, obj/_obj_scaling, f, if_write_to_optim_file);
    }

#if MAST_ENABLE_SNOPT == 1
    MAST::FunctionEvaluation::funobj
    get_objective_evaluation_function() {
    
        return _optim_obj;
    }

    MAST::FunctionEvaluation::funcon
    get_constraint_evaluation_function() {
    
        return _optim_con;
    }
#endif
    
    
    //
    // \section  ex_5_initialization  Initialization
    //
    //   \subsection ex_5_constructor  Constructor
    //
    
    TopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in,
                                   MAST::Examples::GetPotWrapper& input):
    MAST::FunctionEvaluation             (comm_in),
    _initialized                         (false),
    _input                               (input),
    _length                              (0.),
    _height                              (0.),
    _obj_scaling                         (0.),
    _stress_penalty                      (0.),
    _perimeter_penalty                   (0.),
    _stress_lim                          (0.),
    _p_val                               (0.),
    _vm_rho                              (0.),
    _ref_eig_val                         (0.),
    _n_eig_vals                          (0),
    _mesh                                (nullptr),
    _level_set_mesh                      (nullptr),
    _eq_sys                              (nullptr),
    _level_set_eq_sys                    (nullptr),
    _sys                                 (nullptr),
    _level_set_sys                       (nullptr),
    _level_set_sys_on_str_mesh           (nullptr),
    _indicator_sys                       (nullptr),
    _sys_init                            (nullptr),
    _level_set_sys_init_on_str_mesh      (nullptr),
    _level_set_sys_init                  (nullptr),
    _indicator_sys_init                  (nullptr),
    _discipline                          (nullptr),
    _indicator_discipline                (nullptr),
    _level_set_discipline                (nullptr),
    _filter                              (nullptr),
    _m_card                              (nullptr),
    _p_card                              (nullptr),
    _level_set_function                  (nullptr),
    _level_set_vel                       (nullptr),
    _output                              (nullptr) {
        
        libmesh_assert(!_initialized);
        
        //
        // call the initialization routines for each component
        //
        _init_fetype();
        _init_mesh();
        _init_system_and_discipline();
        _init_dirichlet_conditions();
        _init_eq_sys();
        _init_material();
        _init_loads();
        _init_section_property();
        _initialized = true;
        
        //
        // ask structure to use Mindlin bending operator
        //
        dynamic_cast<MAST::ElementPropertyCard2D&>(*_p_card).set_bending_model(MAST::MINDLIN);
        
        /////////////////////////////////////////////////
        // now initialize the design data.
        /////////////////////////////////////////////////
        
        //
        // first, initialize the level set functions over the domain
        //
        this->initialize_solution();
        
        //
        // next, define a new parameter to define design variable for nodal level-set
        // function value
        //
        this->_init_phi_dvs();
        
        _obj_scaling           = 1./_length/_height;
        _stress_penalty        = _input("stress_penalty", "penalty value for stress_constraint", 0.);
        _perimeter_penalty     = _input("perimeter_penalty", "penalty value for perimeter in the objective function", 0.);
        _stress_lim            = _input("vm_stress_limit", "limit von-mises stress value", 2.e8);
        _p_val                 = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0);
        _vm_rho                = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0);
        _level_set_vel         = new MAST::LevelSetBoundaryVelocity(2);
        _level_set_function    = new PhiMeshFunction;
        _output                = new libMesh::ExodusII_IO(*_mesh);
        
        _n_eig_vals            = _input("n_eig", "number of eigenvalues to constrain", 0);
        if (_n_eig_vals) {
            //
            // set only if the user requested eigenvalue constraints
            //
            _ref_eig_val           = _input("eigenvalue_low_bound", "lower bound enforced on eigenvalue constraints", 1.e3);
            _sys->set_n_requested_eigenvalues(_n_eig_vals);
        }
        
        //
        // two inequality constraints: stress and eigenvalue.
        //
        _n_ineq = 1+_n_eig_vals;
        
        std::string
        output_name = _input("output_file_root", "prefix of output file names", "output");
        output_name += "_optim_history.txt";
        this->set_output_file(output_name);
        
    }
    
    //
    //   \subsection ex_5_destructor  Destructor
    //
    ~TopologyOptimizationLevelSet2D() {
        
        {
            std::set<MAST::BoundaryConditionBase*>::iterator
            it   = _boundary_conditions.begin(),
            end  = _boundary_conditions.end();
            for ( ; it!=end; it++)
                delete *it;
        }
        
        {
            std::set<MAST::FunctionBase*>::iterator
            it   = _field_functions.begin(),
            end  = _field_functions.end();
            for ( ; it!=end; it++)
                delete *it;
        }
        
        {
            std::map<std::string, MAST::Parameter*>::iterator
            it   = _parameters.begin(),
            end  = _parameters.end();
            for ( ; it!=end; it++)
                delete it->second;
        }
        
        if (!_initialized)
            return;
        
        delete _m_card;
        delete _p_card;
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _sys_init;
        
        delete _level_set_function;
        delete _level_set_vel;
        delete _level_set_sys_init;
        delete _indicator_sys_init;
        delete _indicator_discipline;
        delete _level_set_discipline;
        delete _filter;
        delete _level_set_eq_sys;
        delete _level_set_mesh;
        delete _output;
        delete _level_set_sys_init_on_str_mesh;
        
        for (unsigned int i=0; i<_dv_params.size(); i++)
            delete _dv_params[i].second;
    }
    

};


//
//   \subsection ex_5_wrappers_snopt  Wrappers for SNOPT
//

TopologyOptimizationLevelSet2D* _my_func_eval = nullptr;

#if MAST_ENABLE_SNOPT == 1

unsigned int
it_num = 0;

void
_optim_obj(int*    mode,
           int*    n,
           double* x,
           double* f,
           double* g,
           int*    nstate) {

    //
    // make sure that the global variable has been setup
    //
    libmesh_assert(_my_func_eval);

    //
    // initialize the local variables
    //
    Real
    obj = 0.;

    unsigned int
    n_vars  =  _my_func_eval->n_vars(),
    n_con   =  _my_func_eval->n_eq()+_my_func_eval->n_ineq();

    libmesh_assert_equal_to(*n, n_vars);

    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (0);

    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), false);
    
    //
    // copy the dvars
    //
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];


    _my_func_eval->_evaluate_wrapper(dvars,
                                     obj,
                                     *mode>0,       // request the derivatives of obj
                                     obj_grad,
                                     fvals,
                                     eval_grads,
                                     grads);

    //
    // now copy them back as necessary
    //
    *f  = obj;
    if (*mode > 0) {
        
        // output data to the file
        _my_func_eval->_output_wrapper(it_num, dvars, obj, fvals, true);
        it_num++;

        for (unsigned int i=0; i<n_vars; i++)
            g[i] = obj_grad[i];
    }

    if (obj > 1.e5) *mode = -1;
}






void
_optim_con(int*    mode,
           int*    ncnln,
           int*    n,
           int*    ldJ,
           int*    needc,
           double* x,
           double* c,
           double* cJac,
           int*    nstate) {

    //
    // make sure that the global variable has been setup
    //
    libmesh_assert(_my_func_eval);

    //
    // initialize the local variables
    //
    Real
    obj = 0.;

    unsigned int
    n_vars  =  _my_func_eval->n_vars(),
    n_con   =  _my_func_eval->n_eq()+_my_func_eval->n_ineq();

    libmesh_assert_equal_to(    *n, n_vars);
    libmesh_assert_equal_to(*ncnln, n_con);

    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (n_vars*n_con, 0.);

    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), *mode>0);

    //
    // copy the dvars
    //
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];


    _my_func_eval->_evaluate_wrapper(dvars,
                                     obj,
                                     false,       // request the derivatives of obj
                                     obj_grad,
                                     fvals,
                                     eval_grads,
                                     grads);

    //
    // now copy them back as necessary
    //
    // first the constraint functions
    //
    for (unsigned int i=0; i<n_con; i++)
        c[i] = fvals[i];

    if (*mode > 0) {
        //
        // next, the constraint gradients
        //
        for (unsigned int i=0; i<n_con*n_vars; i++)
            cJac[i] = grads[i];
    }
    
    if (obj > 1.e5) *mode = -1;
}
#endif

//
//   \subsection ex_5_main Main function
//

int main(int argc, char* argv[]) {

    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    TopologyOptimizationLevelSet2D top_opt(init.comm(), input);
    _my_func_eval = &top_opt;
    
    //MAST::NLOptOptimizationInterface optimizer(NLOPT_LD_SLSQP);
    std::unique_ptr<MAST::OptimizationInterface>
    optimizer;
    
    std::string
    s          = input("optimizer", "optimizer to use in the example", "gcmma");

    if (s == "gcmma") {

        optimizer.reset(new MAST::GCMMAOptimizationInterface);
        
        unsigned int
        max_inner_iters        = input("max_inner_iters", "maximum inner iterations in GCMMA", 15);
        
        Real
        constr_penalty         = input("constraint_penalty", "constraint penalty in GCMMA", 50.),
        initial_rel_step       = input("initial_rel_step", "initial step size in GCMMA", 1.e-2),
        asymptote_reduction    = input("asymptote_reduction", "reduction of aymptote in GCMMA", 0.7),
        asymptote_expansion    = input("asymptote_expansion", "expansion of asymptote in GCMMA", 1.2);
        
        optimizer->set_real_parameter   ("constraint_penalty",  constr_penalty);
        optimizer->set_real_parameter   ("initial_rel_step",  initial_rel_step);
        optimizer->set_real_parameter   ("asymptote_reduction",  asymptote_reduction);
        optimizer->set_real_parameter   ("asymptote_expansion",  asymptote_expansion);
        optimizer->set_integer_parameter(   "max_inner_iters", max_inner_iters);
    }
    else if (s == "snopt") {
        
        optimizer.reset(new MAST::NPSOLOptimizationInterface);
    }
    else {
        
        libMesh::out
        << "Unrecognized optimizer specified: " << s << std::endl;
        libmesh_error();
    }
    
    if (optimizer.get()) {
        
        optimizer->attach_function_evaluation_object(top_opt);

        //std::vector<Real> xx1(top_opt.n_vars()), xx2(top_opt.n_vars());
        //top_opt.init_dvar(xx1, xx2, xx2);
        //top_opt.initialize_dv_from_output_file("output1.txt", 24, xx1);
        //top_opt.verify_gradients(xx1);
        optimizer->optimize();
        //top_opt.parametric_line_study("output1.txt", 0, 450, 500);
    }
    
    // END_TRANSLATE
    return 0;
}
