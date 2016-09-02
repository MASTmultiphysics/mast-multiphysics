/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include <iostream>
#include <ostream>


// MAST includes
#include "examples/fsi/beam_flutter_nonuniform_aero_base_solution/beam_euler_fsi_flutter_nonuniform_aero_base_solution.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "base/nonlinear_system.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "fluid/small_disturbance_pressure_function.h"
#include "solver/complex_solver_base.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/flexible_surface_motion.h"
#include "boundary_condition/function_surface_motion.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/ug_flutter_root.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "aeroelasticity/ug_flutter_solver.h"
#include "aeroelasticity/displacement_function_base.h"
#include "examples/base/augment_ghost_elem_send_list.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_solver.h"


extern libMesh::LibMeshInit* __init;


namespace MAST {
    class PanelShape:
    public MAST::DisplacementFunctionBase {
        
    public:
        
        PanelShape(Real                n_sin_index,
                   Real                t_by_c,
                   Real                length):
        MAST::DisplacementFunctionBase("panel_shape"),
        _n          (n_sin_index),
        _t_by_c     (t_by_c),
        _length     (length),
        _pi         (acos(-1.)) { }
        
        
        virtual ~PanelShape()   { }
        
        
        virtual void operator ()(const libMesh::Point& p,
                                 const Real t,
                                 RealVectorX& w) const {
            
            w.setZero();
            w(1) = .5*_t_by_c* sin(_n*_pi*p(0)/_length);
        }

        
        virtual void dwdx(const libMesh::Point& p,
                          const Real t,
                          RealVectorX& w) const {
            
            w.setZero();
            w(1) = .5*_t_by_c*_n*_pi/_length * cos(_n*_pi*p(0)/_length);
        }

        
        virtual void dwdy(const libMesh::Point& p,
                          const Real t,
                          RealVectorX& w) const {
            
            w.setZero();
        }

        
        virtual void dwdz(const libMesh::Point& p,
                          const Real t,
                          RealVectorX& w) const {
            
            w.setZero();
        }

        
    protected:
        
        unsigned int _n;
        Real         _length;
        Real         _t_by_c;
        Real         _pi;
    };
    
}




MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis::
BeamEulerFSIFlutterNonuniformAeroBaseAnalysis():
_structural_mesh                        (nullptr),
_fluid_mesh                             (nullptr),
_structural_eq_sys                      (nullptr),
_fluid_eq_sys                           (nullptr),
_structural_sys                         (nullptr),
_fluid_sys                              (nullptr),
_structural_sys_init                    (nullptr),
_structural_discipline                  (nullptr),
_fluid_sys_init                         (nullptr),
_fluid_discipline                       (nullptr),
_flight_cond                            (nullptr),
_far_field                              (nullptr),
_symm_wall                              (nullptr),
_slip_wall                              (nullptr),
_pressure                               (nullptr),
_panel_shape                            (nullptr),
_steady_motion_function                 (nullptr),
_unsteady_motion_function               (nullptr),
_small_dist_pressure_function           (nullptr),
_omega                                  (nullptr),
_velocity                               (nullptr),
_b_ref                                  (nullptr),
_omega_f                                (nullptr),
_velocity_f                             (nullptr),
_b_ref_f                                (nullptr),
_zero_freq_func                         (nullptr),
_freq_function                          (nullptr),
_length                                 (0.),
_k_lower                                (0.),
_k_upper                                (0.),
_n_k_divs                               (0),
_thy                                    (nullptr),
_thz                                    (nullptr),
_rho                                    (nullptr),
_E                                      (nullptr),
_nu                                     (nullptr),
_zero                                   (nullptr),
_mach                                   (nullptr),
_rho_air                                (nullptr),
_gamma_air                              (nullptr),
_thy_f                                  (nullptr),
_thz_f                                  (nullptr),
_rho_f                                  (nullptr),
_E_f                                    (nullptr),
_nu_f                                   (nullptr),
_hyoff_f                                (nullptr),
_hzoff_f                                (nullptr),
_mach_f                                 (nullptr),
_rho_air_f                              (nullptr),
_gamma_air_f                            (nullptr),
_flutter_solver                         (nullptr),
_flutter_root                           (nullptr),
_m_card                                 (nullptr),
_p_card                                 (nullptr),
_dirichlet_left                         (nullptr),
_dirichlet_right                        (nullptr),
_augment_send_list_obj                  (nullptr) {
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE FLUID DATA
    //////////////////////////////////////////////////////////////////////

    // initialize the libMesh object
    _fluid_mesh              = new libMesh::ParallelMesh(__init->comm());
    _fluid_eq_sys            = new libMesh::EquationSystems(*_fluid_mesh);
    
    // add the system to be used for analysis
    _fluid_sys = &(_fluid_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    _fluid_sys->set_init_B_matrix();
    
    
    // initialize the flow conditions
    GetPot infile("input.in");
    
    
    const unsigned int
    dim                 = 2,
    n_sin_index         = infile("n_sin_index", 1),
    nx_divs             = 3,
    ny_divs             = 1,
    panel_bc_id         = 10,
    symmetry_bc_id      = 11;
    
    Real
    t_c                 = infile("t_by_c", 0.01);
    
    
    libMesh::ElemType
    elem_type           =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "QUAD4"));
    
    libMesh::FEFamily
    fe_type             =
    libMesh::Utility::string_to_enum<libMesh::FEFamily>(infile("fe_family", "LAGRANGE"));
    
    libMesh::Order
    fe_order            =
    libMesh::Utility::string_to_enum<libMesh::Order>(infile("fe_order", "FIRST"));
    
    std::vector<Real>
    x_div_loc        (nx_divs+1),
    x_relative_dx    (nx_divs+1),
    y_div_loc        (ny_divs+1),
    y_relative_dx    (ny_divs+1);
    
    std::vector<unsigned int>
    x_divs           (nx_divs),
    y_divs           (ny_divs);
    
    std::auto_ptr<MAST::MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions);
    
    std::vector<MAST::MeshInitializer::CoordinateDivisions*>
    divs(dim);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
        
        if (i_div < nx_divs) //  this is only till nx_divs
            x_divs[i_div]       = infile( "x_div_nelem", 0, i_div);
    }
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    
    
    // now read in the values: y-coord
    for (unsigned int i_div=0; i_div<ny_divs+1; i_div++) {
        
        y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
        y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
        
        if (i_div < ny_divs) //  this is only till ny_divs
            y_divs[i_div]    = infile( "y_div_nelem",  0, i_div);
    }
    
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    
    
    
    // initialize the mesh
    MAST::PanelMesh2D().init(0.,               // t/c
                             false,            // if cos bump
                             0,                // n max bumps
                             panel_bc_id,
                             symmetry_bc_id,
                             divs,
                             *_fluid_mesh,
                             elem_type);
    
    _fluid_discipline   = new MAST::ConservativeFluidDiscipline(*_fluid_eq_sys);
    _fluid_sys_init     = new MAST::ConservativeFluidSystemInitialization(*_fluid_sys,
                                                                          _fluid_sys->name(),
                                                                          libMesh::FEType(fe_order, fe_type),
                                                                          dim);
    
    _augment_send_list_obj = new MAST::AugmentGhostElementSendListObj(*_fluid_sys);
    _fluid_sys->get_dof_map().attach_extra_send_list_object(*_augment_send_list_obj);
    
    // initialize the equation system for analysis
    _fluid_eq_sys->init();
    
    // print the information
    _fluid_eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
    _symm_wall     = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    _slip_wall     = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
    
    _flight_cond    =  new MAST::FlightCondition;
    for (unsigned int i=0; i<3; i++) {
        
        _flight_cond->body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
        _flight_cond->body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
        _flight_cond->body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
        _flight_cond->body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
        _flight_cond->body_angular_rates(i) = infile("body_angular_rates", 0., i);
    }
    
    _flight_cond->ref_chord       = infile("ref_c",    1.);
    _flight_cond->altitude        = infile( "alt",     0.);
    _flight_cond->mach            = infile("mach",     .5);
    _flight_cond->gas_property.cp = infile(  "cp",  1003.);
    _flight_cond->gas_property.cv = infile(  "cv",   716.);
    _flight_cond->gas_property.T  = infile("temp",   300.);
    _flight_cond->gas_property.rho= infile( "rho",   1.05);
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _fluid_discipline->set_flight_condition(*_flight_cond);
    
    // define parameters
    _omega             = new MAST::Parameter("omega",       0.);
    _zero_omega        = new MAST::Parameter("omega",       0.);
    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
    _b_ref             = new MAST::Parameter("b_ref",       1.);
    
    
    // now define the constant field functions based on this
    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
    _zero_omega_f      = new MAST::ConstantFieldFunction("omega",  *_zero_omega);
    
    // initialize the frequency function
    _zero_freq_func    = new MAST::FrequencyFunction("freq",
                                                     *_zero_omega_f,
                                                     *_velocity_f,
                                                     *_b_ref_f);
    
    _freq_function     = new MAST::FrequencyFunction("freq",
                                                     *_omega_f,
                                                     *_velocity_f,
                                                     *_b_ref_f);
    _freq_function->if_nondimensional(true);
    
    
    // tell the physics about boundary conditions
    _fluid_discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _fluid_discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _fluid_discipline->add_side_load(              i, *_far_field);
    
    _small_dist_pressure_function =
    new MAST::SmallDisturbancePressureFunction(*_fluid_sys_init, *_flight_cond);

    
    _k_upper            = infile("k_upper",  0.75);
    _k_lower            = infile("k_lower",  0.05);
    _n_k_divs           = infile("n_k_divs",   10);
    
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE STRUCTURAL DATA
    //////////////////////////////////////////////////////////////////////
    
    x_div_loc.resize     (2);
    x_relative_dx.resize (2);
    x_divs.resize        (1);
    divs.resize          (1);
    
    x_coord_divs.reset   (new MeshInitializer::CoordinateDivisions);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<2; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div+1);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div+1);
    }
    x_divs[0]       = infile( "x_div_nelem", 0, 1);
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(1, x_div_loc, x_relative_dx, x_divs);

    
    // setup length for use in setup of flutter solver
    _length = x_div_loc[1]-x_div_loc[0];
    (*_b_ref) = _length;
    
    // create the mesh
    _structural_mesh       = new libMesh::SerialMesh(__init->comm());
    
    
    MeshInitializer().init(divs, *_structural_mesh, libMesh::EDGE2);
    
    // create the equation system
    _structural_eq_sys    = new  libMesh::EquationSystems(*_structural_mesh);
    
    // create the libmesh system
    _structural_sys       = &(_structural_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _structural_sys->set_eigenproblem_type(libMesh::GHEP);
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys_init  = new MAST::StructuralSystemInitialization(*_structural_sys,
                                                                     _structural_sys->name(),
                                                                     fetype);
    _structural_discipline = new MAST::StructuralDiscipline(*_structural_eq_sys);
    
    
    // create and add the boundary condition and loads
    _dirichlet_left = new MAST::DirichletBoundaryCondition;
    _dirichlet_right= new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 3;  // tx
    _dirichlet_left->init (0, constrained_vars);
    _dirichlet_right->init(1, constrained_vars);
    _structural_discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _structural_discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _structural_discipline->init_system_dirichlet_bc(*_structural_sys);
    
    // initialize the equation system
    _structural_eq_sys->init();
    
    // initialize the motion object
    _panel_shape                = new MAST::PanelShape(n_sin_index,
                                                       t_c,
                                                       _length);
    _steady_motion_function     = new MAST::FunctionSurfaceMotion("motion");
    _unsteady_motion_function   = new MAST::FlexibleSurfaceMotion("small_disturbance_motion",
                                                                  *_structural_sys_init);
    _slip_wall->add(*_steady_motion_function);
    _slip_wall->add(*_unsteady_motion_function);
    
    // initialize the steady motion function object
    _steady_motion_function->init(*_zero_freq_func, *_panel_shape);
    
    
    
    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _structural_sys->set_exchange_A_and_B(true);
    _structural_sys->set_n_requested_eigenvalues(infile("n_modes", 10));
    
    // create the property functions and add them to the
    
    _thy             = new MAST::Parameter("thy",  0.0015);
    _thz             = new MAST::Parameter("thz",    1.00);
    _rho             = new MAST::Parameter("rho",   2.7e3);
    _E               = new MAST::Parameter("E",     72.e9);
    _nu              = new MAST::Parameter("nu",     0.33);
    _zero            = new MAST::Parameter("zero",     0.);
    _mach            = new MAST::Parameter("mach",     3.);
    _rho_air         = new MAST::Parameter("rho" ,   1.05);
    _gamma_air       = new MAST::Parameter("gamma",   1.4);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_thy);
    _params_for_sensitivity.push_back(_thz);
    
    
    
    _thy_f           = new MAST::ConstantFieldFunction("hy",          *_thy);
    _thz_f           = new MAST::ConstantFieldFunction("hz",          *_thz);
    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off",     *_zero);
    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off",     *_zero);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);
    
    // create the material property card
    _m_card          = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_rho_f);
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    
    // create the element property card
    _p_card          = new MAST::Solid1DSectionElementPropertyCard;
    
    // tell the card about the orientation
    libMesh::Point orientation;
    orientation(1) = 1.;
    _p_card->y_vector() = orientation;
    
    // add the section properties to the card
    _p_card->add(*_thy_f);
    _p_card->add(*_thz_f);
    _p_card->add(*_hyoff_f);
    _p_card->add(*_hzoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    
    _p_card->init();
    
    _structural_discipline->set_property_for_subdomain(0, *_p_card);
    
    // pressure boundary condition for the beam
    _pressure    =  new MAST::BoundaryConditionBase(MAST::SMALL_DISTURBANCE_MOTION);
    _pressure->add(*_small_dist_pressure_function);
    _pressure->add(*_unsteady_motion_function);
    _structural_discipline->add_volume_load(0, *_pressure);

    _flutter_solver  = new MAST::UGFlutterSolver;
}






MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis::
~BeamEulerFSIFlutterNonuniformAeroBaseAnalysis() {
    
    delete _fluid_eq_sys;
    delete _fluid_mesh;
    
    delete _fluid_discipline;
    delete _fluid_sys_init;
    
    delete _far_field;
    delete _symm_wall;
    delete _slip_wall;
    
    delete _flight_cond;
    
    delete _omega;
    delete _zero_omega;
    delete _velocity;
    delete _b_ref;
    
    delete _omega_f;
    delete _zero_omega_f;
    delete _velocity_f;
    delete _b_ref_f;
    
    delete _zero_freq_func;
    delete _freq_function;

    delete _small_dist_pressure_function;
    
    delete _panel_shape;
    delete _steady_motion_function;
    delete _unsteady_motion_function;
    delete _pressure;
    
    delete _structural_eq_sys;
    delete _structural_mesh;
    
    
    delete _structural_discipline;
    delete _structural_sys_init;
    
    delete _m_card;
    delete _p_card;
    
    delete _dirichlet_left;
    delete _dirichlet_right;
    
    delete _thy_f;
    delete _thz_f;
    delete _rho_f;
    delete _E_f;
    delete _nu_f;
    delete _hyoff_f;
    delete _hzoff_f;
    delete _mach_f;
    delete _rho_air_f;
    delete _gamma_air_f;
    
    
    delete _thy;
    delete _thz;
    delete _rho;
    delete _E;
    delete _nu;
    delete _zero;
    delete _mach;
    delete _rho_air;
    delete _gamma_air;
    
    
    // delete the basis vectors
    if (_basis.size())
        for (unsigned int i=0; i<_basis.size(); i++)
            if (_basis[i]) delete _basis[i];

    delete _flutter_solver;
    delete _augment_send_list_obj;
}





MAST::Parameter*
MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis::
get_parameter(const std::string &nm) {
    
    MAST::Parameter *rval = nullptr;
    
    // look through the vector of parameters to see if the name is available
    std::vector<MAST::Parameter*>::iterator
    it   =  _params_for_sensitivity.begin(),
    end  =  _params_for_sensitivity.end();
    
    bool
    found = false;
    
    for ( ; it != end; it++) {
        
        if (nm == (*it)->name()) {
            rval    = *it;
            found   = true;
        }
    }
    
    // if the param was not found, then print the message
    if (!found) {
        libMesh::out
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            libMesh::out << "   " << (*it)->name() << std::endl;
        libMesh::out << std::endl;
    }
    
    return rval;
}





Real
MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis::
solve(bool if_write_output,
      const Real tol,
      const unsigned int max_bisection_iters) {
    
    
    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID SOLUTION
    /////////////////////////////////////////////////////////////////
    // the modal and flutter problems are solved on rank 0, while
    // the fluid solution is setup on the global communicator
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(4);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_e();
    
    // create the vector for storing the base solution.
    // we will swap this out with the system solution, initialize and
    // then swap it back.
    libMesh::NumericVector<Real>& base_sol =
    _fluid_sys->add_vector("fluid_base_solution");
    _fluid_sys->solution->swap(base_sol);
    _fluid_sys_init->initialize_solution(s);
    
    /////////////////////////////////////////////////////////////////
    // Fluid steady-state solution
    /////////////////////////////////////////////////////////////////
    
    // create the nonlinear assembly object
    MAST::ConservativeFluidTransientAssembly   steady_fluid_assembly;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver     transient_fluid_solver;
    
    // now solve the system
    steady_fluid_assembly.attach_discipline_and_system(*_fluid_discipline,
                                                       transient_fluid_solver,
                                                       *_fluid_sys_init);
    
    libMesh::NonlinearImplicitSystem&      nonlin_sys   =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_fluid_sys->system());
    
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_fluid_mesh);
    
    // time solver parameters
    unsigned int
    t_step            = 0,
    n_iters_change_dt = 4,
    iter_count_dt     = 0,
    max_time_steps    = 1e4;
    
    Real
    tval           = 0.,
    vel_0          = 0.,
    vel_1          = 1.e+12,
    p              = 0.5,
    factor         = 0.,
    min_factor     = 1.5,
    time_step_size = 1.e-6;
    
    transient_fluid_solver.dt         = time_step_size;
    transient_fluid_solver.beta       = 1.0;
    
    // set the previous state to be same as the current state to account for
    // zero velocity as the initial condition
    transient_fluid_solver.solution(1).zero();
    transient_fluid_solver.solution(1).add(1., transient_fluid_solver.solution());
    transient_fluid_solver.solution(1).close();
    
    
    if (if_write_output)
        libMesh::out << "Writing output to : output.exo" << std::endl;
    
    // loop over time steps
    while ((t_step <= max_time_steps) &&
           (vel_1  >=  1.e-2)) {
        
        // change dt if the iteration count has increased to threshold
        if (iter_count_dt == n_iters_change_dt) {
            
            libMesh::out
            << "Changing dt:  old dt = " << transient_fluid_solver.dt
            << "    new dt = ";
            
            factor                        = std::pow(vel_0/vel_1, p);
            factor                        = std::max(factor, min_factor);
            transient_fluid_solver.dt    *= factor;
            
            libMesh::out << transient_fluid_solver.dt << std::endl;
            
            iter_count_dt = 0;
            vel_0         = vel_1;
        }
        else
            iter_count_dt++;
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << tval
        << " :  xdot-L2 = " << vel_1
        << std::endl;
        
        // write the time-step
        if (if_write_output) {
            
            exodus_writer.write_timestep("output.exo",
                                         *_fluid_eq_sys,
                                         t_step+1,
                                         nonlin_sys.time);
        }
        
        transient_fluid_solver.solve();
        
        transient_fluid_solver.advance_time_step();
        
        // get the velocity L2 norm
        vel_1 = transient_fluid_solver.velocity().l2_norm();
        if (t_step == 0) vel_0 = vel_1;
        
        tval  += transient_fluid_solver.dt;
        t_step++;
    }
    
    steady_fluid_assembly.clear_discipline_and_system();

    
    
    // swap back the calculate solution so that we store it for later use.
    _fluid_sys->solution->swap(base_sol);
    
    // create the nonlinear assembly object
    MAST::FrequencyDomainLinearizedComplexAssembly   complex_fluid_assembly;
    
    // solver for complex solution
    MAST::ComplexSolverBase                          complex_fluid_solver;
    
    // now setup the assembly object
    complex_fluid_assembly.attach_discipline_and_system(*_fluid_discipline,
                                          complex_fluid_solver,
                                          *_fluid_sys_init);
    complex_fluid_assembly.set_base_solution(base_sol);
    complex_fluid_assembly.set_frequency_function(*_freq_function);
    
    
    
    
    ////////////////////////////////////////////////////////////
    // STRUCTURAL MODAL EIGENSOLUTION
    ////////////////////////////////////////////////////////////
    
    // create the nonlinear assembly object
    MAST::StructuralModalEigenproblemAssembly   modal_assembly;
    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
    
    modal_assembly.attach_discipline_and_system(*_structural_discipline,
                                                *_structural_sys_init);
    
    
    MAST::StructuralNearNullVectorSpace nsp;
    _structural_sys->nonlinear_solver->nearnullspace_object = &nsp;
    
    _structural_sys->eigenproblem_solve();
    modal_assembly.clear_discipline_and_system();
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_structural_sys->get_n_converged_eigenvalues(),
                     _structural_sys->get_n_requested_eigenvalues());
    
    if (_basis.size() > 0)
        libmesh_assert(_basis.size() == nconv);
    else {
        _basis.resize(nconv);
        for (unsigned int i=0; i<_basis.size(); i++)
            _basis[i] = nullptr;
    }
    
    libMesh::ExodusII_IO*
    writer = nullptr;
    
    if (if_write_output)
        writer = new libMesh::ExodusII_IO(*_structural_mesh);

    
    for (unsigned int i=0; i<nconv; i++) {
        
        // create a vector to store the basis
        if (_basis[i] == nullptr)
            _basis[i] = _structural_sys->solution->zero_clone().release();
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        _structural_sys->get_eigenpair(i, re, im, *_basis[i]);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        if (if_write_output) {
            
            // We write the file in the ExodusII format.
            // copy the solution for output
            _structural_sys->solution->swap(*_basis[i]);
            writer->write_timestep("modes.exo",
                                   *_structural_eq_sys,
                                   i+1, i);
            _structural_sys->solution->swap(*_basis[i]);
        }
    }

    ///////////////////////////////////////////////////////////////////
    // FLUTTER SOLUTION
    ///////////////////////////////////////////////////////////////////
    // clear flutter solver and set the output file
    _flutter_solver->clear();

    MAST::FSIGeneralizedAeroForceAssembly fsi_assembly;
    fsi_assembly.attach_discipline_and_system(*_structural_discipline,
                                              *_structural_sys_init);
    
    std::ostringstream oss;
    oss << "flutter_output_" << __init->comm().rank() << ".txt";
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(oss.str());
    

    fsi_assembly.init(*_freq_function,
                      &complex_fluid_solver,         // fluid complex solver
                      _small_dist_pressure_function,
                      _unsteady_motion_function);
    _flutter_solver->attach_assembly(fsi_assembly);
    _flutter_solver->initialize(*_omega,
                                *_b_ref,
                                _flight_cond->rho(),
                                _k_lower,         // lower kr
                                _k_upper,         // upper kr
                                _n_k_divs,        // number of divisions
                                _basis);          // basis vectors
    
    
    // find the roots for the specified divisions
    _flutter_solver->scan_for_roots();
    _flutter_solver->print_crossover_points();
    
    // now ask the flutter solver to return the critical flutter root,
    // which is the flutter cross-over point at the lowest velocity
    std::pair<bool, MAST::FlutterRootBase*>
    sol = _flutter_solver->find_critical_root(tol, max_bisection_iters);
    
    
    _flutter_solver->print_sorted_roots();
    fsi_assembly.clear_discipline_and_system();
    _flutter_solver->clear_assembly_object();
    
    // make sure solution was found
    libmesh_assert(sol.first);
    _flutter_root = sol.second;

    
    if (if_write_output) {
        
        
        // now write the flutter mode to an output file.
        // Flutter mode Y = sum_i (X_i * (xi_re + xi_im)_i)
        // using the right eigenvector of the system.
        // where i is the structural mode
        //
        // The time domain simulation assumes the temporal solution to be
        // X(t) = (Y_re + i Y_im) exp(p t)
        //      = (Y_re + i Y_im) exp(p_re t) * (cos(p_im t) + i sin(p_im t))
        //      = exp(p_re t) (Z_re + i Z_im ),
        // where Z_re = Y_re cos(p_im t) - Y_im sin(p_im t), and
        //       Z_im = Y_re sin(p_im t) + Y_im cos(p_im t).
        //
        // We write the simulation of the mode over a period of oscillation
        //
        
        // first write the structural mode
        // first calculate the real and imaginary vectors
        std::auto_ptr<libMesh::NumericVector<Real> >
        re(_structural_sys->solution->zero_clone().release()),
        im(_structural_sys->solution->zero_clone().release());
        
        
        // first the real part
        _structural_sys->solution->zero();
        for (unsigned int i=0; i<_basis.size(); i++) {
            re->add(sol.second->eig_vec_right(i).real(), *_basis[i]);
            im->add(sol.second->eig_vec_right(i).imag(), *_basis[i]);
        }
        re->close();
        im->close();
        
        // now open the output processor for writing
        libMesh::ExodusII_IO flutter_mode_output(*_structural_mesh);
        
        // use N steps in a time-period
        Real
        t_sys = _structural_sys->time,
        pi    = acos(-1.);
        unsigned int
        N_divs = 100;
        
        
        for (unsigned int i=0; i<=N_divs; i++) {
            _structural_sys->time   =  2.*pi*(i*1.)/(N_divs*1.);
            
            _structural_sys->solution->zero();
            _structural_sys->solution->add( cos(_structural_sys->time), *re);
            _structural_sys->solution->add(-sin(_structural_sys->time), *im);
            _structural_sys->solution->close();
            flutter_mode_output.write_timestep("structural_flutter_mode.exo",
                                               *_structural_eq_sys,
                                               i+1,
                                               _structural_sys->time);
            
            // reset the system time
            _structural_sys->time = t_sys;
        }
        
        
        
        // next write the fluid mode
        {
            // get the size of basis from the structural node
            unsigned int
            n_basis = (unsigned int)_basis.size();
            _fluid_sys->comm().broadcast(n_basis);
            
            // first calculate the fluid basis
            std::vector<libMesh::NumericVector<Real>*>
            fluid_basis_re(n_basis),
            fluid_basis_im(n_basis);
            for (unsigned int i=0; i<n_basis; i++) {
                
                // solve the fluid system for the given structural mode and
                // the frequency of the flutter root
                _unsteady_motion_function->init(*_freq_function, *_basis[i]);
                
                complex_fluid_solver.solve_block_matrix();
                
                fluid_basis_re[i] = complex_fluid_solver.real_solution().clone().release();
                fluid_basis_im[i] = complex_fluid_solver.imag_solution().clone().release();
            }
            
            
            // first calculate the real and imaginary vectors
            std::auto_ptr<libMesh::NumericVector<Real> >
            re(_fluid_sys->solution->zero_clone().release()),
            im(_fluid_sys->solution->zero_clone().release());
            
            
            // first the real part
            _fluid_sys->solution->zero();
            for (unsigned int i=0; i<n_basis; i++) {
                
                re->add( sol.second->eig_vec_right(i).real(), *fluid_basis_re[i]);
                re->add(-sol.second->eig_vec_right(i).imag(), *fluid_basis_im[i]);
                
                im->add(sol.second->eig_vec_right(i).imag(), *fluid_basis_re[i]);
                im->add(sol.second->eig_vec_right(i).real(), *fluid_basis_im[i]);
            }
            re->close();
            im->close();
            
            // now open the output processor for writing
            libMesh::ExodusII_IO flutter_mode_output(*_fluid_mesh);
            
            // use N steps in a time-period
            Real
            t_sys = _fluid_sys->time,
            pi    = acos(-1.);
            unsigned int
            N_divs = 100;
            
            
            for (unsigned int i=0; i<=N_divs; i++) {
                _fluid_sys->time   =  2.*pi*(i*1.)/(N_divs*1.);
                
                _fluid_sys->solution->zero();
                _fluid_sys->solution->add( cos(_fluid_sys->time), *re);
                _fluid_sys->solution->add(-sin(_fluid_sys->time), *im);
                _fluid_sys->solution->close();
                flutter_mode_output.write_timestep("fluid_flutter_mode.exo",
                                                   *_fluid_eq_sys,
                                                   i+1,
                                                   _fluid_sys->time);
            }
            
            // reset the system time
            _fluid_sys->time = t_sys;
            
            // now delete the solution
            for (unsigned int i=0; i<n_basis; i++) {
                delete fluid_basis_re[i];
                delete fluid_basis_im[i];
            }
        }
    }
    

    
    
    complex_fluid_assembly.clear_discipline_and_system();
    
    return _flutter_root->V;
}





Real
MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis::
sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(false);
    
    //Make sure that  a solution is available for sensitivity
    libmesh_assert(_flutter_root);

    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID ASSEMBLY/SOLVER OBJECTS
    /////////////////////////////////////////////////////////////////
    libMesh::NumericVector<Real>& base_sol =
    _fluid_sys->get_vector("fluid_base_solution");
    
    // create the nonlinear assembly object
    MAST::FrequencyDomainLinearizedComplexAssembly   assembly;
    
    // Transient solver for time integration
    MAST::ComplexSolverBase                          solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(*_fluid_discipline,
                                          solver,
                                          *_fluid_sys_init);
    assembly.set_base_solution(base_sol);
    assembly.set_frequency_function(*_freq_function);


    // it is assumed that the modal basis is available from the flutter solution
    
    
    ///////////////////////////////////////////////////////////////////
    // FLUTTER SOLUTION SENSITIVITY
    ///////////////////////////////////////////////////////////////////
    // clear flutter solver and set the output file
    
    MAST::FSIGeneralizedAeroForceAssembly fsi_assembly;
    
    fsi_assembly.attach_discipline_and_system(*_structural_discipline,
                                              *_structural_sys_init);
    
    std::ostringstream oss;
    oss << "flutter_output_" << __init->comm().rank() << ".txt";
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(oss.str());
    
    
    fsi_assembly.init(*_freq_function,
                      &solver,                       // fluid complex solver
                      _small_dist_pressure_function,
                      _unsteady_motion_function);
    _flutter_solver->attach_assembly(fsi_assembly);

    // flutter solver will need velocity to be defined as a parameter for
    // sensitivity analysis
    _structural_discipline->add_parameter(p);
    _fluid_discipline->add_parameter(*_omega);
    
    libMesh::ParameterVector params;
    params.resize(1);
    params[0]  =  p.ptr();
    
    // calculate the sensitivity
    _flutter_solver->calculate_sensitivity(*_flutter_root, params, 0);

    
    // clean up before exiting
    fsi_assembly.clear_discipline_and_system();
    _flutter_solver->clear_assembly_object();
    _structural_discipline->remove_parameter(p);
    _fluid_discipline->remove_parameter(*_omega);
    
    
    
    return _flutter_root->V_sens;
}

