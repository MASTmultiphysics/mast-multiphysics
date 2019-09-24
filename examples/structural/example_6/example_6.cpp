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
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/ks_stress_output.h"
#include "elasticity/smooth_ramp_stress_output.h"
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/compliance_output.h"
#include "elasticity/structural_system_initialization.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/nonlinear_implicit_assembly.h"
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
// BEGIN_TRANSLATE 2D SIMP topology optimization
//
//   \tableofcontents
//
//  This example computes the optimal topology of a structure subject to
//  specified boundary conditions (Dirichlet and Neumann). An element-wise density
//  is used to parameterize the topology.
//
//  Elasticity function with the penalty term
class ElasticityFunction:
public MAST::FieldFunction<Real> {
public:
    ElasticityFunction(Real E0, Real rho_min, Real penalty,
                       MAST::MeshFieldFunction& rho,
                       MAST::MeshFieldFunction& drho):
    MAST::FieldFunction<Real>("E"),
    _E0(E0),
    _rho_min(rho_min),
    _penalty(penalty),
    _rho(rho),
    _drho(drho) { }
    virtual ~ElasticityFunction(){}
    void set_penalty_val(Real penalty) {_penalty = penalty;}
    
    virtual bool depends_on(const MAST::FunctionBase& f) const { return true;}
    virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {

        RealVectorX v1;
        _rho(p, t, v1);
        
        v = _E0 * (_rho_min + (1.-_rho_min) * pow(v1(0), _penalty));
    }

    virtual void derivative(const MAST::FunctionBase& f,
                            const libMesh::Point& p, const Real t, Real& v) const {
        
        RealVectorX v1, dv1;
        _rho(p, t, v1);
        _drho(p, t, dv1);
        
        v = _E0 * (1.-_rho_min) * _penalty * pow(v1(0), _penalty-1.) * dv1(0);
    }

protected:
    Real                    _E0; // value of the material Young's modulus
    Real                    _rho_min; // lower limit on density
    Real                    _penalty; // value of penalty term
    MAST::MeshFieldFunction &_rho;
    MAST::MeshFieldFunction &_drho;
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



class TopologyOptimizationSIMP2D:
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
    Real                                      _vf;      // volume fraction
    Real                                      _rho_min; // lower limit on density

    ElasticityFunction*                       _Ef;
    libMesh::UnstructuredMesh*                _mesh;
    
    libMesh::EquationSystems*                 _eq_sys;
    
    MAST::NonlinearSystem*                    _sys;
    libMesh::System*                          _density_sys;
    
    MAST::StructuralSystemInitialization*     _sys_init;
    
    MAST::PhysicsDisciplineBase*              _discipline;
    
    MAST::FilterBase*                         _filter;
    
    MAST::MaterialPropertyCardBase*           _m_card;
    MAST::ElementPropertyCardBase*            _p_card;
    
    MAST::MeshFieldFunction*                  _density_function;
    MAST::MeshFieldFunction*                  _density_sens_function;
    libMesh::ExodusII_IO*                     _output;
    
    libMesh::FEType                           _fetype;
    libMesh::FEType                           _density_fetype;
    
    std::vector<MAST::Parameter*>             _params_for_sensitivity;
    std::map<std::string, MAST::Parameter*>   _parameters;
    std::set<MAST::FunctionBase*>             _field_functions;
    std::set<MAST::BoundaryConditionBase*>    _boundary_conditions;
    std::set<unsigned int>                    _dv_dof_ids;
    
    std::vector<std::pair<unsigned int, MAST::Parameter*>>  _dv_params;

public:
    
    //  \section  ex_6_init_mesh Mesh Generation
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
    //  \subsection  ex_6_inplane_mesh Inplane problem
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
    }
    
    //
    //  \subsection  ex_6_bracket_mesh Bracket
    //
    void _init_mesh_bracket() {

        {
            unsigned int
            nx_divs = _input("nx_divs", "number of elements along x-axis", 20),
            ny_divs = _input("ny_divs", "number of elements along y-axis", 20);
            
            if (nx_divs%10 != 0 || ny_divs%10 != 0) libmesh_error();
        }
        
        _init_mesh_inplane();
        _delete_elems_from_bracket_mesh(*_mesh);
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
    //  \subsection  ex_6_eyebar_mesh Eyebar
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
    }

    //
    //  \section  ex_6_system_discipline  System and Discipline
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
        _density_fetype      = libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE);
        _density_sys         = &(_eq_sys->add_system<libMesh::ExplicitSystem>("density"));
        _density_sys->add_variable("rho", _density_fetype);
    }

    
    void _init_eq_sys() {
        
        _eq_sys->init();
        _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
        _sys->set_exchange_A_and_B(true);
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
    //  \section  ex_6_dirichlet Dirichlet Constraints
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
    //  \subsection  ex_6_inplane_dirichlet Inplane
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
    }
    
    //
    //  \subsection  ex_6_truss_dirichlet Truss
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
    }

    
    //
    //  \subsection  ex_6_bracket_dirichlet Bracket
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
    }
    
    
    //
    //  \subsection  ex_6_eyebar_dirichlet Eyebar
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
    }

    

    //
    //  \section  ex_6_loading Loading
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
    
    
    //  \subsection  ex_6_inplane_loading Inplane
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
        *press_f         = new FluxLoad( "pressure", p_val, _length, frac);
        
        // initialize the load
        MAST::BoundaryConditionBase
        *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
        
        p_load->add(*press_f);
        _discipline->add_side_load(2, *p_load);
        
        _field_functions.insert(press_f);
    }
    
    
    //
    //  \subsection  ex_6_bracket_loading Bracket
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
        *press_f         = new BracketLoad( "pressure", p_val, length, frac);
        
        //
        // initialize the load
        //
        MAST::BoundaryConditionBase
        *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
        
        p_load->add(*press_f);
        _discipline->add_side_load(5, *p_load);
        
        _field_functions.insert(press_f);
    }

    
    //
    //  \subsection  ex_6_eyebar_loading Eyebar
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
    //   \section  ex_6_properties Properties
    //
    //
    //
    //   \subsection  ex_6_material_properties Material Properties
    //

    void _init_material() {
        
        _rho_min  = _input("rho_min", "lower limit on density variable", 1.e-8);

        Real
        Eval      = _input("E", "modulus of elasticity", 72.e9),
        penalty   = _input("rho_penalty", "SIMP modulus of elasticity penalty", 4.),
        rhoval    = _input("rho", "material density", 2700.),
        nu_val    = _input("nu", "Poisson's ratio",  0.33),
        kappa_val = _input("kappa", "shear correction factor",  5./6.),
        kval      = _input("k", "thermal conductivity",  1.e-2),
        cpval     = _input("cp", "thermal capacitance",  864.);
        
        
        MAST::Parameter
        *rho       = new MAST::Parameter("rho",      rhoval),
        *nu        = new MAST::Parameter("nu",       nu_val),
        *kappa     = new MAST::Parameter("kappa", kappa_val),
        *k         = new MAST::Parameter("k",          kval),
        *cp        = new MAST::Parameter("cp",        cpval);
        
        MAST::ConstantFieldFunction
        *rho_f   = new MAST::ConstantFieldFunction(  "rho",    *rho),
        *nu_f    = new MAST::ConstantFieldFunction(   "nu",     *nu),
        *kappa_f = new MAST::ConstantFieldFunction("kappa",  *kappa),
        *k_f     = new MAST::ConstantFieldFunction( "k_th",      *k),
        *cp_f    = new MAST::ConstantFieldFunction(   "cp",     *cp);

        _Ef      = new ElasticityFunction(Eval, _rho_min, penalty,
                                          *_density_function,
                                          *_density_sens_function);
        
        _parameters[  rho->name()]     = rho;
        _parameters[   nu->name()]     = nu;
        _parameters[kappa->name()]     = kappa;
        _parameters[    k->name()]     = k;
        _parameters[   cp->name()]     = cp;
        _field_functions.insert(_Ef);
        _field_functions.insert(rho_f);
        _field_functions.insert(nu_f);
        _field_functions.insert(kappa_f);
        _field_functions.insert(k_f);
        _field_functions.insert(cp_f);

        _m_card = new MAST::IsotropicMaterialPropertyCard;
        _m_card->add(*_Ef);
        _m_card->add(*rho_f);
        _m_card->add(*nu_f);
        _m_card->add(*kappa_f);
        _m_card->add(*k_f);
        _m_card->add(*cp_f);
    }

    
    //
    //   \subsection  ex_6_section_properties Section Properties
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
    }
    

    //
    //   \section  ex_6_initial_solution Initial Density field
    //
    //
    
    
    void initialize_solution() {
        
        //
        // initialize density field to a constant value of the specified
        // volume fraction
        //
        _vf    = _input("volume_fraction", "upper limit for the voluem fraction", 0.5);

        _density_sys->solution->zero();
        _density_sys->solution->add(_vf);
        _density_sys->solution->close();
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
        _filter                = new MAST::FilterBase(*_density_sys, filter_radius, _dv_dof_ids);
        libMesh::NumericVector<Real>& vec = _density_sys->add_vector("base_values");
        vec = *_density_sys->solution;
        vec.close();
    }
    
    //
    //  \subsection  ex_6_inplane_initial_level_set Inplane
    //
    void _init_phi_dvs_inplane() {
        
        //
        // this assumes that density variable has a constant value per element
        //
        libmesh_assert_equal_to(_density_fetype.family, libMesh::LAGRANGE);

        Real
        frac          = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        sys_num = _density_sys->number(),
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_mesh->is_replicated());
        
        std::vector<Real> local_phi(_density_sys->solution->size());
        _density_sys->solution->localize(local_phi);
        
        // iterate over all the element values
        libMesh::MeshBase::const_node_iterator
        it  = _mesh->nodes_begin(),
        end = _mesh->nodes_end();
        
        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_mesh->n_elem());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;

            dof_id                     = n.dof_number(sys_num, 0, 0);
            
            // only if node is not on the upper edge
            if ((n(1)+filter_radius >= _height) &&
                (n(0)-filter_radius <= _length*.5*(1.+frac))   &&
                (n(0)+filter_radius >= _length*.5*(1.-frac))) {
                
                // set value at the material points to a small positive number
                if (dof_id >= _density_sys->solution->first_local_index() &&
                    dof_id <  _density_sys->solution->last_local_index())
                    _density_sys->solution->set(dof_id, 1.e0);
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
        
        _density_sys->solution->close();
    }
    
    //
    //  \subsection  ex_6_truss_initial_level_set Truss
    //
    void _init_phi_dvs_truss() {
        
        //
        // this assumes that density variable has a constant value per element
        //
        libmesh_assert_equal_to(_density_fetype.family, libMesh::LAGRANGE);

        Real
        frac          = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        sys_num = _density_sys->number(),
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_mesh->is_replicated());
        
        std::vector<Real> local_phi(_density_sys->solution->size());
        _density_sys->solution->localize(local_phi);
        
        // iterate over all the element values
        libMesh::MeshBase::const_node_iterator
        it  = _mesh->nodes_begin(),
        end = _mesh->nodes_end();

        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_mesh->n_elem());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;

            dof_id                     = n.dof_number(sys_num, 0, 0);
            
            // only if node is not on the upper edge
            if ((n(1)-filter_radius <= 0.) &&
                (n(0)-filter_radius <= _length*.5*(1.+frac))   &&
                (n(0)+filter_radius >= _length*.5*(1.-frac))) {
                
                // set value at the material points to a small positive number
                if (dof_id >= _density_sys->solution->first_local_index() &&
                    dof_id <  _density_sys->solution->last_local_index())
                    _density_sys->solution->set(dof_id, 1.e0);
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
        
        _density_sys->solution->close();
    }
    
    //
    //  \subsection  ex_6_bracket_initial_level_set Bracket
    //
    void _init_phi_dvs_bracket() {
        
        libmesh_assert(_initialized);

        //
        // this assumes that density variable has a constant value per element
        //
        libmesh_assert_equal_to(_density_fetype.family, libMesh::LAGRANGE);

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
        sys_num = _density_sys->number(),
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_mesh->is_replicated());
        
        std::vector<Real> local_phi(_density_sys->solution->size());
        _density_sys->solution->localize(local_phi);
        
        // iterate over all the element values
        libMesh::MeshBase::const_node_iterator
        it  = _mesh->nodes_begin(),
        end = _mesh->nodes_end();

        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_mesh->n_elem());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;

            dof_id                     = n.dof_number(sys_num, 0, 0);
            
            if ((n(1)-filter_radius) <= y_lim && (n(0)+filter_radius) >= length*(1.-frac)) {
          
                //
                // set value at the constrained points to a small positive number
                // material here
                //
                if (dof_id >= _density_sys->solution->first_local_index() &&
                    dof_id <  _density_sys->solution->last_local_index())
                    _density_sys->solution->set(dof_id, 1.e0);
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
                    
                    if (dof_id >= _density_sys->solution->first_local_index() &&
                        dof_id <  _density_sys->solution->last_local_index())
                        _density_sys->solution->set(dof_id, _rho_min);
                    val = _rho_min;
                }
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);
                
                _n_vars++;
            }
        }
        
        _density_sys->solution->close();
    }
    
    
    //
    //  \subsection  ex_6_eyebar_initial_level_set Eyebar
    //
    void _init_phi_dvs_eye_bar() {
        
        libmesh_assert(_initialized);

        //
        // this assumes that density variable has a constant value per element
        //
        libmesh_assert_equal_to(_density_fetype.family, libMesh::LAGRANGE);

        Real
        tol           = 1.e-6,
        filter_radius = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        
        unsigned int
        sys_num = _density_sys->number(),
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        //
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        //
        libmesh_assert(_mesh->is_replicated());
        
        std::vector<Real> local_phi(_density_sys->solution->size());
        _density_sys->solution->localize(local_phi);
        
        // iterate over all the element values
        // iterate over all the element values
        libMesh::MeshBase::const_node_iterator
        it  = _mesh->nodes_begin(),
        end = _mesh->nodes_end();

        //
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        //
        _dv_params.reserve(_mesh->n_elem());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(sys_num, 0, 0);
            
            
            
            if (((n.norm() <= 1.5+filter_radius) && n(0) <= 0.) ||  // circle
                (std::fabs(n(0)-_height*1.5) < filter_radius &&  // right edge
                 std::fabs(n(1)) <= 1.+filter_radius)) { // dirichlet constraint
                    
                    //
                    // set value at the constrained points to material
                    //
                    if (dof_id >= _density_sys->solution->first_local_index() &&
                        dof_id <  _density_sys->solution->last_local_index())
                        _density_sys->solution->set(dof_id, 1.e0);
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
                    
                    if (dof_id >= _density_sys->solution->first_local_index() &&
                        dof_id <  _density_sys->solution->last_local_index())
                        _density_sys->solution->set(dof_id, _rho_min);
                    val = _rho_min;
                }
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::LevelSetParameter(oss.str(), val, &n);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                _dv_dof_ids.insert(dof_id);
                
                _n_vars++;
            }
        }
        
        _density_sys->solution->close();
    }

    
    //
    //   \subsection  ex_6_design_variable_init   Design Variables
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
        
        std::fill(xmin.begin(), xmin.end(),    _rho_min);
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
    //  \subsection  ex_6_function_evaluation Function Evaluation
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
        &base_phi = _density_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, dvars[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_density_sys->solution);
        _density_function->clear();
        _density_function->init(*_density_sys->solution);
        _sys->solution->zero();
        
        //*********************************************************************
        // DO NOT zero out the gradient vector, since GCMMA needs it for the  *
        // subproblem solution                                                *
        //*********************************************************************
        MAST::NonlinearImplicitAssembly                          nonlinear_assembly;
        MAST::StressAssembly                                     stress_assembly;
        MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
        
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        
        MAST::StressStrainOutputBase                    stress;
        MAST::ComplianceOutput                          compliance;
        stress.set_discipline_and_system(*_discipline, *_sys_init);
        stress.set_participating_elements_to_all();
        stress.set_aggregation_coefficients(_p_val, 1., _vm_rho, _stress_lim) ;
        compliance.set_participating_elements_to_all();
        compliance.set_discipline_and_system(*_discipline, *_sys_init);

        //////////////////////////////////////////////////////////////////////
        // evaluate the stress constraint
        //////////////////////////////////////////////////////////////////////
        
        libMesh::out << "Static Solve" << std::endl;

        Real
        penalty          = _input("rho_penalty", "SIMP modulus of elasticity penalty", 4.),
        stress_penalty   = _input("stress_rho_penalty", "SIMP modulus of elasticity penalty for stress evaluation", 0.5);
        // set the elasticity penalty for solution
        _Ef->set_penalty_val(penalty);
        
        _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
        SNESConvergedReason
        r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
        (*_sys->nonlinear_solver).get_converged_reason();
        
        // if the solver diverged due to linear solve, then there is a problem with
        // this geometry and we need to return with a high value set for the
        // constraints
        if (r == SNES_DIVERGED_LINEAR_SOLVE ||
            _sys->final_nonlinear_residual() > 1.e-1) {
            
            obj = 1.e11;
            for (unsigned int i=0; i<_n_ineq; i++)
                fvals[i] = 1.e11;
            return;
        }
        
        // evaluate compliance
        //nonlinear_assembly.calculate_output(*_sys->solution, compliance);

        // set the elasticity penalty for stress evaluation
        _Ef->set_penalty_val(stress_penalty);
        nonlinear_assembly.calculate_output(*_sys->solution, stress);
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the objective
        //////////////////////////////////////////////////////////////////////
        Real
        max_vm = stress.get_maximum_von_mises_stress(),
        vm_agg = stress.output_total(),
        vol    = 0.,
        comp   = compliance.output_total();
        
        _evaluate_volume(&vol, nullptr);
        
        //obj       = _obj_scaling * comp;
        obj         = _obj_scaling * vol;
        //_obj_scaling    * (vol+ _perimeter_penalty * per) +
        //_stress_penalty * (vm_agg);///_stress_lim - 1.);
        
        fvals[0]  =  stress.output_total()/_stress_lim - 1.;  // g = sigma/sigma0-1 <= 0
        //fvals[0]  =  stress.output_total();  // g = sigma/sigma0-1 <= 0
        //fvals[0]  = vol/_length/_height - _vf; // vol/vol0 - a <=
        libMesh::out << "volume: " << vol << std::endl;
        libMesh::out << "max: "    << max_vm << "  constr: " << vm_agg///_stress_lim - 1.
        << std::endl;
        libMesh::out << "compliance: " << comp << std::endl;

        //////////////////////////////////////////////////////////////////////
        // evaluate the objective sensitivities, if requested
        //////////////////////////////////////////////////////////////////////
        if (eval_obj_grad) {
            
            _evaluate_volume(nullptr, &obj_grad);
            for (unsigned int i=0; i<grads.size(); i++) obj_grad[i] /= (_length*_height);
//            std::vector<Real>
//            grad1(obj_grad.size(), 0.);
//
//            _evaluate_compliance_sensitivity(compliance,
//                                             nonlinear_elem_ops,
//                                             nonlinear_assembly,
//                                             grad1);
//
//            for (unsigned int i=0; i<obj_grad.size(); i++)
//                obj_grad[i] += _obj_scaling * grad1[i];
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
        if (if_grad_sens) {
            
            _evaluate_stress_sensitivity(penalty,
                                         stress_penalty,
                                         stress,
                                         nonlinear_elem_ops,
                                         nonlinear_assembly,
                                         grads);

            //_evaluate_volume(nullptr, &grads);
            //for (unsigned int i=0; i<grads.size(); i++) grads[i] /= (_length*_height);
        }
        
        //
        // also the stress data for plotting
        //
        stress_assembly.update_stress_strain_data(stress, *_sys->solution);
        
        _density_function->clear();
    }

    //
    //  \subsection  ex_6_volume_sensitivity Sensitivity of Material Volume
    //
    void _evaluate_volume(Real               *volume,
                          std::vector<Real>  *grad) {

        libMesh::DofMap
        &dof_map  =  _density_sys->get_dof_map();
        
        const std::vector<libMesh::dof_id_type>
        &send_list = dof_map.get_send_list();
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        local_sol(libMesh::NumericVector<Real>::build(_density_sys->comm()).release()),
        local_dsol(libMesh::NumericVector<Real>::build(_density_sys->comm()).release());
        
        local_sol->init(_density_sys->n_dofs(),
                        _density_sys->n_local_dofs(),
                        send_list,
                        false,
                        libMesh::GHOSTED);

        local_dsol->init(_density_sys->n_dofs(),
                         _density_sys->n_local_dofs(),
                         send_list,
                         false,
                         libMesh::GHOSTED);

        
        if (volume) {

            *volume = 0.;
            
            unsigned int
            sys_num = _density_sys->number();
            
            _density_sys->solution->localize(*local_sol, send_list);
            
            libMesh::MeshBase::element_iterator
            it    =  _mesh->active_local_elements_begin(),
            end   =  _mesh->active_local_elements_end();
            
            Real
            rho = 0.;
            
            for ( ; it != end; it++) {
                
                const libMesh::Elem& e = **it;
                
                // compute the average element density value
                rho = 0.;
                for (unsigned int i=0; i<e.n_nodes(); i++) {
                    const libMesh::Node& n = *e.node_ptr(i);
                    rho += local_sol->el(n.dof_number(sys_num, 0, 0));
                }
                rho /= (1. * e.n_nodes());

                // use this density value to compute the volume
                *volume  +=  e.volume() * rho;
            }
            
            this->comm().sum(*volume);
        }
        
        
        if (grad) {
            
            std::fill(grad->begin(), grad->end(), 0.);
            
            //
            // iterate over each DV, create a sensitivity vector and calculate the
            // volume sensitivity explicitly
            //
            std::unique_ptr<libMesh::NumericVector<Real>>
            dphi_base(_density_sys->solution->zero_clone().release()),
            dphi_filtered(_density_sys->solution->zero_clone().release());
            
            ElementParameterDependence dep(*_filter);
            
            for (unsigned int i=0; i<_n_vars; i++) {
                
                dphi_base->zero();
                dphi_filtered->zero();
                local_dsol->zero();
                //
                // set the value only if the dof corresponds to a local node
                //
                if (_dv_params[i].first >=  dphi_base->first_local_index() &&
                    _dv_params[i].first <   dphi_base->last_local_index())
                    dphi_base->set(_dv_params[i].first, 1.);
                dphi_base->close();
                _filter->compute_filtered_values(*dphi_base, *dphi_filtered);
                
                dphi_filtered->localize(*local_dsol, send_list);
                
                unsigned int
                sys_num = _density_sys->number();
                
                libMesh::MeshBase::element_iterator
                it    =  _mesh->active_local_elements_begin(),
                end   =  _mesh->active_local_elements_end();
                
                Real
                drho = 0.;
                
                for ( ; it != end; it++) {
                    
                    const libMesh::Elem& e = **it;
                    
                    // do not compute if the element is not in the domain
                    // of influence of the parameter
                    if (!dep.if_elem_depends_on_parameter(e, *_dv_params[i].second))
                        continue;
                    
                    // compute the average element density value
                    drho = 0.;
                    for (unsigned int i=0; i<e.n_nodes(); i++) {
                        const libMesh::Node& n = *e.node_ptr(i);
                        drho += local_dsol->el(n.dof_number(sys_num, 0, 0));
                    }
                    drho /= (1. * e.n_nodes());
                    
                    // use this density value to compute the volume
                    (*grad)[i]  +=  e.volume() * drho;
                }
                
                this->comm().sum((*grad)[i]);
            }
        }
    }
    
    
    
    //
    //  \subsection  ex_6_stress_sensitivity Sensitivity of Stress and Eigenvalues
    //
    void
    _evaluate_stress_sensitivity
    (const Real                    penalty,
     const Real                    stress_penalty,
     MAST::StressStrainOutputBase& stress,
     MAST::AssemblyElemOperations& nonlinear_elem_ops,
     MAST::NonlinearImplicitAssembly& nonlinear_assembly,
     std::vector<Real>& grads) {
        
        _sys->adjoint_solve(nonlinear_elem_ops, stress, nonlinear_assembly, false);
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi_base(_density_sys->solution->zero_clone().release()),
        dphi_filtered(_density_sys->solution->zero_clone().release());
        
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
            
            _density_sens_function->init(*dphi_filtered);
            
            //////////////////////////////////////////////////////////////////////
            // stress sensitivity
            //////////////////////////////////////////////////////////////////////
            // set the elasticity penalty for solution, which is needed for
            // computation of the residual sensitivity
            _Ef->set_penalty_val(penalty);

            grads[1*i+0] = 1./_stress_lim*
            nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                    _sys->get_adjoint_solution(),
                                                                    *_dv_params[i].second,
                                                                    nonlinear_elem_ops,
                                                                    stress,
                                                                    false);
            
            _Ef->set_penalty_val(stress_penalty);
            nonlinear_assembly.calculate_output_direct_sensitivity(*_sys->solution,
                                                                   nullptr,
                                                                   *_dv_params[i].second,
                                                                   stress);
            grads[1*i+0] += 1./_stress_lim* stress.output_sensitivity_total(*_dv_params[i].second);
            
            stress.clear_sensitivity_data();
            _density_sens_function->clear();
        }
        
        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    
    void
    _evaluate_compliance_sensitivity
    (MAST::ComplianceOutput&                  compliance,
     MAST::AssemblyElemOperations&            nonlinear_elem_ops,
     MAST::NonlinearImplicitAssembly&         nonlinear_assembly,
     std::vector<Real>& grads) {
        
        // Adjoint solution for compliance = - X
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi_base(_density_sys->solution->zero_clone().release()),
        dphi_filtered(_density_sys->solution->zero_clone().release());
        
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
            
            _density_sens_function->init(*dphi_filtered);
            
            //////////////////////////////////////////////////////////////////////
            // compliance sensitivity
            //////////////////////////////////////////////////////////////////////
            grads[i] = -1. *
            nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                    *_sys->solution,
                                                                    *_dv_params[i].second,
                                                                    nonlinear_elem_ops,
                                                                    compliance);
            _density_sens_function->clear();
        }
        
        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    //
    //  \subsection  ex_6_design_output  Output of Design Iterate
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
        &base_phi = _density_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, x[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_density_sys->solution);
        _density_function->init(*_density_sys->solution);
        
        std::vector<bool> eval_grads(this->n_ineq(), false);
        std::vector<Real> f(this->n_ineq(), 0.), grads;
        this->evaluate(x, obj, false, grads, f, eval_grads, grads);
        
        _sys->time = iter;
        _sys_init->get_stress_sys().time = iter;
        // "1" is the number of time-steps in the file, as opposed to the time-step number.
        libMesh::ExodusII_IO(*_mesh).write_timestep(oss.str(), *_eq_sys, 1, (1.*iter));
        
        _density_function->clear();
        
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

        MAST::FunctionEvaluation::output(iter, x, obj/_obj_scaling, fval, if_write_to_optim_file);
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
    // \section  ex_6_initialization  Initialization
    //
    //   \subsection  ex_6_constructor  Constructor
    //
    
    TopologyOptimizationSIMP2D(const libMesh::Parallel::Communicator& comm_in,
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
    _vf                                  (0.),
    _rho_min                             (0.),
    _mesh                                (nullptr),
    _eq_sys                              (nullptr),
    _sys                                 (nullptr),
    _sys_init                            (nullptr),
    _discipline                          (nullptr),
    _filter                              (nullptr),
    _m_card                              (nullptr),
    _p_card                              (nullptr),
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
        
        // density function is used by elasticity modulus function. So, we
        // initialize this here
        _density_function        = new MAST::MeshFieldFunction(*_density_sys, "rho");
        _density_sens_function   = new MAST::MeshFieldFunction(*_density_sys, "rho");

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
        _output                = new libMesh::ExodusII_IO(*_mesh);
        
        //
        // two inequality constraints: stress and eigenvalue.
        //
        _n_ineq = 1;
        
        std::string
        output_name = _input("output_file_root", "prefix of output file names", "output");
        output_name += "_optim_history.txt";
        this->set_output_file(output_name);
        
    }
    
    //
    //   \subsection  ex_6_destructor  Destructor
    //
    ~TopologyOptimizationSIMP2D() {
        
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
        
        delete _filter;
        delete _output;
        delete _density_function;
        delete _density_sens_function;

        for (unsigned int i=0; i<_dv_params.size(); i++)
            delete _dv_params[i].second;
    }
    

};


//
//   \subsection  ex_6_wrappers_snopt  Wrappers for SNOPT
//

TopologyOptimizationSIMP2D* _my_func_eval = nullptr;

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

    if (obj > 1.e10) *mode = -1;
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
    
    if (obj > 1.e10) *mode = -1;
}
#endif

//
//   \subsection  ex_6_main Main function
//

int main(int argc, char* argv[]) {

    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    TopologyOptimizationSIMP2D top_opt(init.comm(), input);
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
        //top_opt.verify_gradients(xx1);
        optimizer->optimize();
    }
    
    // END_TRANSLATE
    return 0;
}
