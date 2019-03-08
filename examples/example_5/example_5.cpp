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


// MAST includes
#include "examples/base/input_wrapper.h"
#include "level_set/level_set_discipline.h"
#include "level_set/level_set_system_initialization.h"
#include "level_set/level_set_eigenproblem_assembly.h"
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_boundary_velocity.h"
#include "level_set/indicator_function_constrain_dofs.h"
#include "level_set/level_set_constrain_dofs.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/structural_system_initialization.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/slepc_eigen_solver.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "optimization/nlopt_optimization_interface.h"
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




class Phi:
public MAST::FieldFunction<RealVectorX> {
    
public:
    Phi(Real l1,
        Real l2,
        Real nx_mesh,
        Real ny_mesh,
        Real nx_holes,
        Real ny_holes):
    MAST::FieldFunction<RealVectorX>("Phi"),
    _l1  (l1),
    _l2  (l2),
    _nx_mesh  (nx_mesh),
    _ny_mesh  (ny_mesh),
    _nx_holes (nx_holes),
    _ny_holes (ny_holes),
    _pi  (acos(-1.)) {
        
        // initialize the locations at which the holes will be nucleated
        // first, along the x-axis
        if (_nx_holes == 1)
            _x_axis_hole_locations.insert(l1 * 0.5);
        else if (_nx_holes >= 2) {
            
            // add holes at the beginning and end
            _x_axis_hole_locations.insert(0.);
            _x_axis_hole_locations.insert(_l1);
            
            // now, add holes at uniformly spaced locations
            // in the domain
            Real
            dx = _l1/(1.*(_nx_holes-1));
            for (unsigned int i=2; i<_nx_holes; i++)
                _x_axis_hole_locations.insert(dx*(i-1));
        }
        
        
        // now, along the y-axis
        if (_ny_holes == 1)
            _y_axis_hole_locations.insert(l2 * 0.5);
        else if (_ny_holes >= 2) {
            
            // add holes at the beginning and end
            _y_axis_hole_locations.insert(0.);
            _y_axis_hole_locations.insert(_l2);
            
            // now, add holes at uniformly spaced locations
            // in the domain
            Real
            dx = _l2/(1.*(_ny_holes-1));
            for (unsigned int i=2; i<_ny_holes; i++)
                _y_axis_hole_locations.insert(dx*(i-1));
        }
    }
    virtual ~Phi() {}
    virtual void operator()(const libMesh::Point& p,
                            const Real t,
                            RealVectorX& v) const {
        
        libmesh_assert_less_equal(t, 1);
        libmesh_assert_equal_to(v.size(), 1);
        
        
        /*// circle
         //v(0) = -(pow(p(0)-_l1*.5, 2) + pow(p(1)-_l2*.5, 2) - pow(_r, 2));
         
         // waves
         Real
         c  = 0.5,
         pi = acos(-1.),
         x  = p(0)-.5*_l1,
         y  = p(1)-.5*_l2,
         r  = pow(pow(x,2)+pow(y,2),.5);
         //v(0) = 1.*cos(nx*r*_pi/_l1);
         v(0) = cos(2.*_nx*pi*x/_l1)+cos(2.*_ny*pi*y/_l2)+c;
         
         // linear
         //v(0) = (p(0)-_l1*0.5)*(-10.);
         //v(0) = (p(0)+p(1)-_l1)*(-10.);*/
        
        // the libMesh solution projection routine for Lagrange elements
        // will query the function value at the nodes. So, we figure
        // out which nodes should have zero values set to them.
        // if there is one hole in any direction, it will be in the
        // center of the domain. If there are more than 1, then two of
        // the holes will be on the boundary and others will fill the
        // interior evenly.
        
        const Real
        tol     = 1.e-6*std::min(_l1, _l2),
        dx_mesh = _l1/(1.*_nx_mesh),
        dy_mesh = _l2/(1.*_ny_mesh);
        
        std::set<Real>::const_iterator
        x_it_low = _x_axis_hole_locations.lower_bound(p(0)-dx_mesh),
        y_it_low = _y_axis_hole_locations.lower_bound(p(1)-dy_mesh);
        
        unsigned int
        n = 0;
        // see if the x-location needs a hole
        for ( ; x_it_low != _x_axis_hole_locations.end(); x_it_low++) {
            if (std::fabs(*x_it_low - p(0)) <= dx_mesh*0.5) {
                n++;
                break;
            }
        }
        
        // now check the y-location
        for ( ; y_it_low != _y_axis_hole_locations.end(); y_it_low++) {
            if (std::fabs(*y_it_low - p(1)) <= dy_mesh*0.5) {
                n++;
                break;
            }
        }
        
        if (n == 2)
            v(0) = -0.01;
        else
            v(0) = 0.01;
    }
protected:
    Real
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

class Vel: public MAST::FieldFunction<Real> {
public:
    Vel(): MAST::FieldFunction<Real>("vel") {}
    
    virtual void operator() (const libMesh::Point& p,
                             const Real t,
                             Real& v) const {
        
        // waves
        Real
        nt = 8.,
        th = atan2(p(0)-.15, p(1)-.15);
        v  = sin(nt*th/2.);
        
        // constant
        // v    = 1.;
    }
    
protected:
    
};


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



class TopologyOptimizationLevelSet2D:
public MAST::FunctionEvaluation {
    
protected:
    
    bool                                      _initialized;
    MAST::Examples::GetPotWrapper&            _input;
    
    Real                                      _length;
    Real                                      _height;
    Real                                      _obj_scaling;
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
    
    MAST::SystemInitialization*               _sys_init;
    MAST::LevelSetSystemInitialization*       _level_set_sys_init_on_str_mesh;
    MAST::LevelSetSystemInitialization*       _level_set_sys_init;
    MAST::HeatConductionSystemInitialization* _indicator_sys_init;
    
    MAST::PhysicsDisciplineBase*              _discipline;
    MAST::PhysicsDisciplineBase*              _indicator_discipline;
    MAST::LevelSetDiscipline*                 _level_set_discipline;
    
    MAST::MaterialPropertyCardBase*           _m_card;
    MAST::ElementPropertyCardBase*            _p_card;
    
    PhiMeshFunction*                          _level_set_function;
    MAST::LevelSetBoundaryVelocity*           _level_set_vel;
    libMesh::ExodusII_IO*                     _output;
    
    libMesh::FEType                           _fetype;
    libMesh::FEType                           _level_set_fetype;
    
    std::vector<MAST::Parameter*>             _params_for_sensitivity;
    std::map<std::string, MAST::Parameter*>   _parameters;
    std::vector<MAST::Parameter*>             _dv_parameters;
    std::set<MAST::FunctionBase*>             _field_functions;
    std::set<MAST::BoundaryConditionBase*>    _boundary_conditions;
    
    
    std::vector<std::pair<unsigned int, MAST::Parameter*>>  _dv_params;
    
    

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
    
    
    
    void _init_mesh()  {
        
        _mesh = new libMesh::SerialMesh(this->comm());
        
        // identify the element type from the input file or from the order
        // of the element
        unsigned int
        nx_divs = _input("nx_divs", "number of elements along x-axis", 10),
        ny_divs = _input("ny_divs", "number of elements along y-axis", 10);
        
        _length = _input("length", "length of domain along x-axis", 0.3),
        _height = _input("height", "length of domain along y-axis", 0.3);
        
        std::string
        t = _input("elem_type", "type of geometric element in the mesh", "quad4");
        
        libMesh::ElemType
        e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
        
        // if high order FE is used, libMesh requires atleast a second order
        // geometric element.
        if (_fetype.order > 1 && e_type == libMesh::QUAD4)
            e_type = libMesh::QUAD9;
        else if (_fetype.order > 1 && e_type == libMesh::TRI3)
            e_type = libMesh::TRI6;
        
        // initialize the mesh with one element
        libMesh::MeshTools::Generation::build_square(*_mesh,
                                                     nx_divs, ny_divs,
                                                     0, _length,
                                                     0, _height,
                                                     e_type);

        // mesh on which the level-set function is defined
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
    
    
    
    void _init_system_and_discipline() {
        
        
        // make sure that the mesh has been initialized
        libmesh_assert(_mesh);
        
        // create the equation system
        _eq_sys    = new  libMesh::EquationSystems(*_mesh);
        
        // create the libmesh system and set the preferences for structural
        // eigenvalue problems
        _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
        _sys->set_eigenproblem_type(libMesh::GHEP);
        
        // initialize the system to the right set of variables
        _sys_init       = new MAST::StructuralSystemInitialization(*_sys,
                                                                   _sys->name(),
                                                                   _fetype);
        _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
        
        
        // Initialize the system for level set function.
        // A level set function is defined on a coarser mesh than the structural
        // mesh.
        // A level set function is assumed to be a first-order Lagrange finite element
        _level_set_fetype      = libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE);
        _level_set_eq_sys      = new libMesh::EquationSystems(*_level_set_mesh);
        _level_set_sys         = &(_level_set_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
        _level_set_sys_init    = new MAST::LevelSetSystemInitialization(*_level_set_sys,
                                                                        _level_set_sys->name(),
                                                                        _level_set_fetype);
        _level_set_discipline  = new MAST::LevelSetDiscipline(*_eq_sys);
        
        // A system with level set function is defined on the strucutral mesh
        // for the purpose of plotting.
        _level_set_sys_on_str_mesh      = &(_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
        _level_set_sys_init_on_str_mesh = new MAST::LevelSetSystemInitialization(*_level_set_sys_on_str_mesh,
                                                                                 _level_set_sys->name(),
                                                                                 _level_set_fetype);
        
        //  an indicator function is used to locate unconnected free-floating
        // domains of material. The indicator function solves a heat condution
        // problem. Regions with uniformly zero temperature are marked as
        // unconnected domains.
        _indicator_sys                  = &(_eq_sys->add_system<MAST::NonlinearSystem>("indicator"));
        _indicator_sys_init             = new MAST::HeatConductionSystemInitialization(*_indicator_sys,
                                                                                       _indicator_sys->name(),
                                                                                       _fetype);
        _indicator_discipline           = new MAST::PhysicsDisciplineBase(*_eq_sys);
    }
    
    
    void _init_dirichlet_conditions() {
        
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
    
    
    
    void _init_eq_sys() {
        
        _eq_sys->init();
        _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
        _sys->set_exchange_A_and_B(true);

        _level_set_eq_sys->init();
    }
    
    
    
    void _init_loads() {
        
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
    
    
    
    void initialize_solution() {
        
        // initialize solution of the level set problem
        unsigned int
        nx_h    = _input("initial_level_set_n_holes_in_x",
                            "number of holes along x-direction for initial level-set field", 2.),
        ny_h    = _input("initial_level_set_n_holes_in_y",
                            "number of holes along y-direction for initial level-set field", 2.),
        nx_m    = _input("level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
        ny_m    = _input("level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
        
        Phi phi(_length, _height, nx_m, ny_m, nx_h, ny_h);
        _level_set_sys_init->initialize_solution(phi);
    }

    
    void _init_phi_dvs() {
        
        // this assumes that level set is defined using lagrange shape functions
        libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
        
        Real
        tol     = 1.e-6,
        frac    = _input("load_length_fraction", "fraction of boundary length on which pressure will act", 0.2);
        
        unsigned int
        dof_id  = 0;
        
        Real
        val     = 0.;
        
        // all ranks will have DVs defined for all variables. So, we should be
        // operating on a replicated mesh
        libmesh_assert(_level_set_mesh->is_replicated());
        
        std::vector<Real> local_phi(_level_set_sys->solution->size());
        _level_set_sys->solution->localize(local_phi);
        
        // iterate over all the node values
        libMesh::MeshBase::const_node_iterator
        it  = _level_set_mesh->nodes_begin(),
        end = _level_set_mesh->nodes_end();
        
        // maximum number of dvs is the number of nodes on the level set function
        // mesh. We will evaluate the actual number of dvs
        _dv_params.reserve(_level_set_mesh->n_nodes());
        _n_vars = 0;
        
        for ( ; it!=end; it++) {
            
            const libMesh::Node& n = **it;
            
            dof_id                     = n.dof_number(0, 0, 0);
            
            // only if node is not on the upper edge
            if ((std::fabs(n(1)-_height) > tol) ||
                (n(0) > _length*.5*(1.+frac))   ||
                (n(0) < _length*.5*(1.-frac))) {
                
                std::ostringstream oss;
                oss << "dv_" << _n_vars;
                val  = local_phi[dof_id];
                
                //            // on the traction free boundary, set everything to be zero, so that there
                //            // is always a boundary there that the optimizer can move
                //            if (n(1) < tol                     ||
                //                std::fabs(n(1) - height) < tol) {
                //
                //                if (dof_id >= _level_set_sys->solution->first_local_index() &&
                //                    dof_id <  _level_set_sys->solution->last_local_index())
                //                    _level_set_sys->solution->set(dof_id, 0.);
                //                val = 0.;
                //            }
                
                _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
                _dv_params[_n_vars].first  = dof_id;
                _dv_params[_n_vars].second = new MAST::Parameter(oss.str(), val);
                _dv_params[_n_vars].second->set_as_topology_parameter(true);
                
                _n_vars++;
            }
            else {
                // set value at the material points to a small positive number
                if (dof_id >= _level_set_sys->solution->first_local_index() &&
                    dof_id <  _level_set_sys->solution->last_local_index())
                    _level_set_sys->solution->set(dof_id, 0.01);
            }
        }
        
        _level_set_sys->solution->close();
    }
    
    
    
    void _evaluate_volume_sensitivity(MAST::LevelSetVolume& volume,
                                      MAST::LevelSetNonlinearImplicitAssembly& assembly,
                                      std::vector<Real>& obj_grad) {
        
        // iterate over each DV, create a sensitivity vector and calculate the
        // volume sensitivity explicitly
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi(_level_set_sys->solution->zero_clone().release());
        
        for (unsigned int i=0; i<_n_vars; i++) {
            
            dphi->zero();
            // set the value only if the dof corresponds to a local node
            if (_dv_params[i].first >=  dphi->first_local_index() &&
                _dv_params[i].first <   dphi->last_local_index())
                dphi->set(_dv_params[i].first, 1.);
            dphi->close();
            
            _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
            
            assembly.calculate_output_direct_sensitivity(*_level_set_sys->solution,
                                                         *dphi,
                                                         *_dv_params[i].second,
                                                         volume);
            obj_grad[i] = _obj_scaling * volume.output_sensitivity_total(*_dv_params[i].second);
        }
    }
    


    
    void
    _evaluate_constraint_sensitivity
    (MAST::StressStrainOutputBase& stress,
     MAST::AssemblyElemOperations& nonlinear_elem_ops,
     MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
     MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
     MAST::LevelSetEigenproblemAssembly& eigen_assembly,
     const std::vector<bool>& eval_grads,
     std::vector<Real>& grads) {
        
        unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());
        
        _sys->adjoint_solve(nonlinear_elem_ops, stress, nonlinear_assembly, false);
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        dphi(_level_set_sys->solution->zero_clone().release());
        
        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<_n_vars; i++) {
            
            dphi->zero();
            dphi->set(_dv_params[i].first, 1.);
            dphi->close();
            
            // initialize the level set perturbation function to create a velocity
            // field
            _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
            
            //////////////////////////////////////////////////////////////////////
            // stress sensitivity
            //////////////////////////////////////////////////////////////////////
            grads[_n_ineq*i+0] = 1./_stress_lim*
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
    }

public:
    
    TopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in,
                                   MAST::Examples::GetPotWrapper& input):
    MAST::FunctionEvaluation             (comm_in),
    _initialized                         (false),
    _input                               (input),
    _length                              (0.),
    _height                              (0.),
    _obj_scaling                         (0.),
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
    _m_card                              (nullptr),
    _p_card                              (nullptr),
    _level_set_function                  (nullptr),
    _level_set_vel                       (nullptr),
    _output                              (nullptr) {
        
        libmesh_assert(!_initialized);
        
        // call the initialization routines for each component
        _init_fetype();
        _init_mesh();
        _init_system_and_discipline();
        _init_dirichlet_conditions();
        _init_eq_sys();
        _init_material();
        _init_loads();
        _init_section_property();
        
        
        // ask structure to use Mindlin bending operator
        dynamic_cast<MAST::ElementPropertyCard2D&>(*_p_card).set_bending_model(MAST::MINDLIN);
        
        /////////////////////////////////////////////////
        // now initialize the design data.
        /////////////////////////////////////////////////
        
        // first, initialize the level set functions over the domain
        this->initialize_solution();
        
        // next, define a new parameter to define design variable for nodal level-set
        // function value
        this->_init_phi_dvs();
        
        _obj_scaling           = 100./_length/_height;
        _stress_lim            = _input("vm_stress_limit", "limit von-mises stress value", 2.e8);
        _p_val                 = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0);
        _vm_rho                = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0);
        _level_set_vel         = new MAST::LevelSetBoundaryVelocity(2);
        _level_set_function    = new PhiMeshFunction;
        _output                = new libMesh::ExodusII_IO(*_mesh);
        
        _n_eig_vals            = _input("n_eig", "number of eigenvalues to constrain", 5);
        if (_n_eig_vals) {
            // set only if the user requested eigenvalue constraints
            _ref_eig_val           = _input("eigenvalue_low_bound", "lower bound enforced on eigenvalue constraints", 1.e3);
            _sys->set_n_requested_eigenvalues(_n_eig_vals);
        }
        
        // two inequality constraints: stress and eigenvalue.
        _n_ineq = 1+_n_eig_vals;
        
        std::string
        output_name = _input("output_file_root", "prefix of output file names", "output");
        output_name += "_optim_history.txt";
        this->set_output_file(output_name);
        
        _initialized = true;
    }
    
    
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
        delete _level_set_eq_sys;
        delete _level_set_mesh;
        delete _output;
        delete _level_set_sys_init_on_str_mesh;
        
        for (unsigned int i=0; i<_dv_params.size(); i++)
            delete _dv_params[i].second;
    }

    
    
    
    /*!
     *   initializes the design variable vector, called by the
     *   optimization interface.
     */
    void init_dvar(std::vector<Real>& x,
                   std::vector<Real>& xmin,
                   std::vector<Real>& xmax) {
        
        // one DV for each element
        x.resize(_n_vars);
        xmin.resize(_n_vars);
        xmax.resize(_n_vars);
        
        std::fill(xmin.begin(), xmin.end(),   -1.);
        std::fill(xmax.begin(), xmax.end(),    1.);
        
        // now, check if the user asked to initialize dvs from a previous file
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
    

    /*!
     *   \p grads(k): Derivative of f_i(x) with respect
     *   to x_j, where k = (j-1)*M + i.
     */
    void evaluate(const std::vector<Real>& dvars,
                  Real& obj,
                  bool eval_obj_grad,
                  std::vector<Real>& obj_grad,
                  std::vector<Real>& fvals,
                  std::vector<bool>& eval_grads,
                  std::vector<Real>& grads) {
        
        libMesh::out << "New Evaluation" << std::endl;
        
        // copy DVs to level set function
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= _level_set_sys->solution->first_local_index() &&
                _dv_params[i].first <  _level_set_sys->solution->last_local_index())
                _level_set_sys->solution->set(_dv_params[i].first, dvars[i]);
        _level_set_sys->solution->close();
        _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
        _sys->solution->zero();
        
        /**********************************************************************
         * DO NOT zero out the gradient vector, since GCMMA needs it for the  *
         * subproblem solution                                                *
         **********************************************************************/
        MAST::LevelSetNonlinearImplicitAssembly                  nonlinear_assembly;
        MAST::LevelSetNonlinearImplicitAssembly                  level_set_assembly;
        MAST::LevelSetEigenproblemAssembly                       eigen_assembly;
        MAST::LevelSetStressAssembly                             stress_assembly;
        MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
        MAST::HeatConductionNonlinearAssemblyElemOperations      conduction_elem_ops;
        MAST::StructuralModalEigenproblemAssemblyElemOperations  modal_elem_ops;
        
        // reinitialize the dof constraints before solution of the linear system
        // FIXME: we should be able to clear the constraint object from the
        // system before it goes out of scope, but libMesh::System does not
        // have a clear method. So, we are going to leave it as is, hoping
        // that libMesh::System will not attempt to use it (most likely, we
        // shoudl be ok).
        
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        SNESConvergedReason r;
        {
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
        //MAST::IndicatorFunctionConstrainDofs constrain(*_sys_init, *_level_set_function, indicator);
        //MAST::LevelSetConstrainDofs constrain(*_sys_init, *_level_set_function);
        //_sys->attach_constraint_object(constrain);
        //_sys->reinit_constraints();
        //_sys->initialize_condensed_dofs(*_discipline);
        
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        nonlinear_assembly.set_level_set_function(*_level_set_function);
        nonlinear_assembly.set_level_set_velocity_function(*_level_set_vel);
        //nonlinear_assembly.set_indicator_function(indicator);
        eigen_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        eigen_assembly.set_level_set_function(*_level_set_function);
        eigen_assembly.set_level_set_velocity_function(*_level_set_vel);
        stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        stress_assembly.init(*_level_set_function, nonlinear_assembly.get_dof_handler());
        level_set_assembly.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
        level_set_assembly.set_level_set_function(*_level_set_function);
        level_set_assembly.set_level_set_velocity_function(*_level_set_vel);
        nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        modal_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        //nonlinear_assembly.plot_sub_elems(true, false, true);
        
        
        
        
        MAST::LevelSetVolume                            volume(level_set_assembly.get_intersection());
        MAST::StressStrainOutputBase                    stress;
        volume.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
        stress.set_discipline_and_system(*_discipline, *_sys_init);
        volume.set_participating_elements_to_all();
        stress.set_participating_elements_to_all();
        stress.set_aggregation_coefficients(_p_val, _vm_rho, _stress_lim);
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the objective
        //////////////////////////////////////////////////////////////////////
        level_set_assembly.calculate_output(*_level_set_sys->solution, volume);
        obj       = volume.output_total() * _obj_scaling;
        
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
        fvals[0]  =  stress.output_total()/_stress_lim - 1.;  // g = sigma/sigma0-1 <= 0
        
        //stress_assembly.update_stress_strain_data(stress, *_sys->solution);
        //libMesh::ExodusII_IO(*_mesh).write_equation_systems("indicator.exo", *_eq_sys);
        //libMesh::ExodusII_IO(*_level_set_mesh).write_equation_systems("phi.exo", *_level_set_eq_sys);
        
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
        if (eval_obj_grad)
            _evaluate_volume_sensitivity(volume, level_set_assembly, obj_grad);
        
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
            _evaluate_constraint_sensitivity(stress,
                                             nonlinear_elem_ops,
                                             nonlinear_assembly,
                                             modal_elem_ops,
                                             eigen_assembly,
                                             eval_grads,
                                             grads);
        
        // also the stress data for plotting
        stress_assembly.update_stress_strain_data(stress, *_sys->solution);
    }

    
    
    //void level_set_solve();
    
    
    void output(unsigned int iter,
                const std::vector<Real>& x,
                Real obj,
                const std::vector<Real>& fval,
                bool if_write_to_optim_file) {
        
        libmesh_assert_equal_to(x.size(), _n_vars);
        
        std::string
        output_name  = _input("output_file_root", "prefix of output file names", "output"),
        modes_name   = output_name + "modes.exo";
        output_name += "_optim.exo";
        
        // copy DVs to level set function
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= _level_set_sys->solution->first_local_index() &&
                _dv_params[i].first <  _level_set_sys->solution->last_local_index())
                _level_set_sys->solution->set(_dv_params[i].first, x[i]);
        _level_set_sys->solution->close();
        _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
        _level_set_sys_init_on_str_mesh->initialize_solution(_level_set_function->get_mesh_function());
        
        std::vector<bool> eval_grads(this->n_ineq(), false);
        std::vector<Real> f(this->n_ineq(), 0.), grads;
        this->evaluate(x, obj, false, grads, f, eval_grads, grads);
        
        _output->write_timestep(output_name, *_eq_sys, iter+1, (1.*iter));
        
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
        
        MAST::FunctionEvaluation::output(iter, x, obj/_obj_scaling, fval, if_write_to_optim_file);
    }

};


//
//
//void
//TopologyOptimizationLevelSet2D::level_set_solve() {
//
//    libmesh_assert(_initialized);
//
//    bool
//    output      = _input("if_output", "if write output to a file", false),
//    propagate   = _input("if_propagate", "if propagate level set, or reinitialize it", true);
//
//
//    std::string
//    output_name = _input("output_file_root", "prefix of output file names", "output");
//    output_name += "_level_set.exo";
//
//    // create the nonlinear assembly object
//    std::unique_ptr<MAST::TransientAssembly> level_set_assembly;
//    if (propagate)
//        level_set_assembly.reset(new MAST::TransientAssembly);
//    else {
//        MAST::LevelSetReinitializationTransientAssembly
//        *assembly = new MAST::LevelSetReinitializationTransientAssembly;
//
//        libMesh::NumericVector<Real>
//        &base_sol = _level_set_sys->add_vector("base_sol");
//        base_sol  = *_level_set_sys->solution;
//
//        assembly->set_reference_solution(base_sol);
//        _level_set_discipline->set_level_set_propagation_mode(false);
//        level_set_assembly.reset(assembly);
//    }
//    MAST::LevelSetTransientAssemblyElemOperations            level_set_elem_ops;
//
//    // Transient solver for time integration
//    MAST::FirstOrderNewmarkTransientSolver  level_set_solver;
//
//    // now solve the system
//    level_set_assembly->set_discipline_and_system(*_level_set_discipline,
//                                                  *_level_set_sys_init);
//
//    // file to write the solution for visualization
//    libMesh::ExodusII_IO exodus_writer(*_mesh);
//
//    // time solver parameters
//    unsigned int
//    t_step                         = 0,
//    n_steps                        = _input("level_set_n_transient_steps", "number of transient time-steps", 100);
//    level_set_solver.dt            = _input("level_set_dt", "time-step size",    1.e-3);
//    level_set_solver.beta          = 0.5;
//
//    // set the previous state to be same as the current state to account for
//    // zero velocity as the initial condition
//    level_set_solver.solution(1).zero();
//    level_set_solver.solution(1).add(1., level_set_solver.solution());
//    level_set_solver.solution(1).close();
//
//
//    // loop over time steps
//    while (t_step <= n_steps) {
//
//        libMesh::out
//        << "Time step: " << t_step
//        << " :  t = " << _level_set_sys->time
//        << std::endl;
//
//        // write the time-step
//        if (output) {
//
//            exodus_writer.write_timestep(output_name,
//                                         *_eq_sys,
//                                         t_step+1,
//                                         _level_set_sys->time);
//        }
//
//        level_set_solver.solve(*level_set_assembly);
//
//        level_set_solver.advance_time_step();
//        t_step++;
//    }
//
//    level_set_assembly->clear_discipline_and_system();
//}

int main(int argc, char* argv[]) {

    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    TopologyOptimizationLevelSet2D top_opt(init.comm(), input);
    
    MAST::NLOptOptimizationInterface optimizer(NLOPT_LD_SLSQP);

    optimizer.attach_function_evaluation_object(top_opt);
    optimizer.optimize();
    
    // END_TRANSLATE
    return 0;
}
