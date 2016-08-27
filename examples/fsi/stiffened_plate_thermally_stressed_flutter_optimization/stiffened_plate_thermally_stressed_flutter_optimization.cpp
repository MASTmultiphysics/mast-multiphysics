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


// MAST includes
#include "examples/fsi/stiffened_plate_thermally_stressed_flutter_optimization/stiffened_plate_thermally_stressed_flutter_optimization.h"
#include "examples/fsi/base/gaf_database.h"
#include "driver/driver_base.h"
#include "optimization/optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/stress_output_base.h"
#include "base/nonlinear_system.h"
#include "examples/fluid/meshing/mesh_initializer.h"
#include "examples/fluid/meshing/panel_mesh_3D.h"
#include "examples/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/flight_condition.h"
#include "fluid/small_disturbance_pressure_function.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "solver/complex_solver_base.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/ug_flutter_solver.h"
#include "boundary_condition/flexible_surface_motion.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "aeroelasticity/flutter_root_base.h"
#include "examples/base/augment_ghost_elem_send_list.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_near_null_vector_space.h"


// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/nonlinear_solver.h"


extern
libMesh::LibMeshInit     *__init;
extern
MAST::FunctionEvaluation *__my_func_eval;



void
stiffened_plate_euler_flutter_optim_obj(int*    mode,
                                        int*    n,
                                        double* x,
                                        double* f,
                                        double* g,
                                        int*    nstate) {
    
    
    // make sure that the global variable has been setup
    libmesh_assert(__my_func_eval);
    
    // initialize the local variables
    Real
    obj = 0.;
    
    unsigned int
    n_vars  =  __my_func_eval->n_vars(),
    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
    
    libmesh_assert_equal_to(*n, n_vars);
    
    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (0);
    
    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), false);
    
    // copy the dvars
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];
    
    
    __my_func_eval->evaluate(dvars,
                             obj,
                             true,       // request the derivatives of obj
                             obj_grad,
                             fvals,
                             eval_grads,
                             grads);
    
    
    // now copy them back as necessary
    *f  = obj;
    for (unsigned int i=0; i<n_vars; i++)
        g[i] = obj_grad[i];
}






void
stiffened_plate_euler_flutter_optim_con(int*    mode,
                                        int*    ncnln,
                                        int*    n,
                                        int*    ldJ,
                                        int*    needc,
                                        double* x,
                                        double* c,
                                        double* cJac,
                                        int*    nstate) {
    
    
    // make sure that the global variable has been setup
    libmesh_assert(__my_func_eval);
    
    // initialize the local variables
    Real
    obj = 0.;
    
    unsigned int
    n_vars  =  __my_func_eval->n_vars(),
    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
    
    libmesh_assert_equal_to(    *n, n_vars);
    libmesh_assert_equal_to(*ncnln, n_con);
    
    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (n_vars*n_con, 0.);
    
    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), true);
    
    // copy the dvars
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];
    
    
    __my_func_eval->evaluate(dvars,
                             obj,
                             true,       // request the derivatives of obj
                             obj_grad,
                             fvals,
                             eval_grads,
                             grads);
    
    
    // now copy them back as necessary
    
    // first the constraint functions
    for (unsigned int i=0; i<n_con; i++)
        c[i] = fvals[i];
    
    // next, the constraint gradients
    for (unsigned int i=0; i<n_con*n_vars; i++)
        cJac[i] = grads[i];
    
    
}



MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::
StiffenedPlateThermallyStressedFSIFlutterSizingOptimization(const libMesh::Parallel::Communicator& comm):
MAST::FunctionEvaluation(comm),
_initialized                            (false),
_length                                 (0.),
_width                                  (0.),
_stress_limit                           (0.),
_k_lower                                (0.),
_k_upper                                (0.),
_V0_flutter                             (0.),
_n_eig                                  (0),
_n_elems                                (0),
_n_divs_between_stiff                   (0),
_n_stiff                                (0),
_n_plate_elems                          (0),
_n_elems_per_stiff                      (0),
_n_stations                             (0),
_n_k_divs                               (0),
_n_dv_stations_x                        (0),
_n_load_steps                           (0),
_structural_mesh                        (nullptr),
_fluid_mesh                             (nullptr),
_structural_eq_sys                      (nullptr),
_fluid_eq_sys                           (nullptr),
_structural_sys                         (nullptr),
_fluid_sys                              (nullptr),
_structural_sys_init                    (nullptr),
_structural_discipline                  (nullptr),
_nsp                                    (nullptr),
_fluid_sys_init                         (nullptr),
_fluid_discipline                       (nullptr),
_frequency_domain_fluid_assembly        (nullptr),
_complex_solver                         (nullptr),
_modal_assembly                         (nullptr),
_structural_nonlinear_assembly          (nullptr),
_flight_cond                            (nullptr),
_far_field                              (nullptr),
_symm_wall                              (nullptr),
_slip_wall                              (nullptr),
_pressure                               (nullptr),
_motion_function                        (nullptr),
_small_dist_pressure_function           (nullptr),
_omega                                  (nullptr),
_velocity                               (nullptr),
_b_ref                                  (nullptr),
_omega_f                                (nullptr),
_velocity_f                             (nullptr),
_b_ref_f                                (nullptr),
_freq_function                          (nullptr),
_flutter_solver                         (nullptr),
_gaf_database                           (nullptr),
_E                                      (nullptr),
_nu                                     (nullptr),
_kappa                                  (nullptr),
_alpha                                  (nullptr),
_rho                                    (nullptr),
_temp                                   (nullptr),
_p_cav                                  (nullptr),
_zero                                   (nullptr),
_E_f                                    (nullptr),
_nu_f                                   (nullptr),
_kappa_f                                (nullptr),
_alpha_f                                (nullptr),
_rho_f                                  (nullptr),
_hoff_plate_f                           (nullptr),
_thzoff_stiff_f                         (nullptr),
_temp_f                                 (nullptr),
_ref_temp_f                             (nullptr),
_p_cav_f                                (nullptr),
_weight                                 (nullptr),
_m_card                                 (nullptr),
_p_card_plate                           (nullptr),
_dirichlet_left                         (nullptr),
_dirichlet_right                        (nullptr),
_dirichlet_bottom                       (nullptr),
_dirichlet_top                          (nullptr),
_T_load                                 (nullptr),
_p_load                                 (nullptr),
_augment_send_list_obj                  (nullptr)
{ }


void
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::
init(GetPot &infile,
     libMesh::ElemType etype,
     bool if_nonlin) {
    
    libmesh_assert(!_initialized);
    
    // number of elements
    _n_elems    = infile("n_elems",   20);
    
    // number of stations
    _n_stations = infile("n_stations", 3);
    
    
    // now setup the optimization data
    _n_vars                = _n_stations; // for thickness variable
    _n_eq                  = 0;
    _n_ineq                = 1;           // one flutter velocity constraint
    _max_iters             = 1000;
    
    
    
    // length of domain
    _length                = infile("length", 10.);
    
    
    _V0_flutter            =  infile("V0_flutter",     410.);

    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE FLUID DATA
    //////////////////////////////////////////////////////////////////////
    
    // initialize the libMesh object
    _fluid_mesh              = new libMesh::ParallelMesh(__init->comm());
    _fluid_eq_sys            = new libMesh::EquationSystems(*_fluid_mesh);
    
    
    // add the system to be used for analysis
    _fluid_sys = &(_fluid_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    _fluid_sys->set_init_B_matrix();
    
    const unsigned int
    dim                 = 3,
    nx_divs             = 3,
    ny_divs             = 3,
    nz_divs             = 1,
    panel_bc_id         = 10,
    symmetry_bc_id      = 11;
    
    libMesh::ElemType
    elem_type           =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "HEX8"));
    
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
    y_relative_dx    (ny_divs+1),
    z_div_loc        (nz_divs+1),
    z_relative_dx    (nz_divs+1);
    
    std::vector<unsigned int>
    x_divs           (nx_divs),
    y_divs           (ny_divs),
    z_divs           (nz_divs);
    
    std::auto_ptr<MAST::MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions),
    z_coord_divs    (new MeshInitializer::CoordinateDivisions);
    
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
    
    
    // now read in the values: z-coord
    for (unsigned int i_div=0; i_div<nz_divs+1; i_div++) {
        
        z_div_loc[i_div]     = infile("z_div_loc", 0., i_div);
        z_relative_dx[i_div] = infile( "z_rel_dx", 0., i_div);
        
        if (i_div < nz_divs) //  this is only till ny_divs
            z_divs[i_div]    = infile( "z_div_nelem",  0, i_div);
    }
    
    divs[2] = z_coord_divs.get();
    z_coord_divs->init(nz_divs, z_div_loc, z_relative_dx, z_divs);

    
    
    // initialize the mesh
    MAST::PanelMesh3D().init(0.,               // t/c
                             false,            // if cos bump
                             0,                // n max bumps in x
                             0,                // n max bumps in y
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
    
    _flight_cond->ref_chord       = infile("ref_c",      1.);
    _flight_cond->altitude        = infile( "alt",       0.);
    _flight_cond->mach            = infile("mach",       .5);
    _flight_cond->gas_property.cp = infile(  "cp",    1003.);
    _flight_cond->gas_property.cv = infile(  "cv",     716.);
    _flight_cond->gas_property.T  = infile("temp",     300.);
    _flight_cond->gas_property.rho= infile( "rho_f",   1.05);
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _fluid_discipline->set_flight_condition(*_flight_cond);
    
    // define parameters
    _omega             = new MAST::Parameter("omega",       0.);
    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
    _b_ref             = new MAST::Parameter("b_ref",       0.);
    
    
    // now define the constant field functions based on this
    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
    
    // initialize the frequency function
    _freq_function     = new MAST::FrequencyFunction("freq",
                                                     *_omega_f,
                                                     *_velocity_f,
                                                     *_b_ref_f);
    _freq_function->if_nondimensional(true);
    
    // tell the physics about boundary conditions
    _fluid_discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _fluid_discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=5; i++)
        _fluid_discipline->add_side_load(              i, *_far_field);
    
    _small_dist_pressure_function =
    new MAST::SmallDisturbancePressureFunction(*_fluid_sys_init, *_flight_cond);
    
    
    _k_upper            = infile("k_upper",  0.75);
    _k_lower            = infile("k_lower",  0.05);
    _n_k_divs           = infile("n_k_divs",   10);
    
    
    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID SOLUTION
    /////////////////////////////////////////////////////////////////
    // the modal and flutter problems are solved on rank 0, while
    // the fluid solution is setup on the global communicator
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(5);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_u3();
    s(4) = _flight_cond->rho_e();
    
    // create the vector for storing the base solution.
    // we will swap this out with the system solution, initialize and
    // then swap it back.
    libMesh::NumericVector<Real>& base_sol =
    _fluid_sys->add_vector("fluid_base_solution");
    _fluid_sys->solution->swap(base_sol);
    _fluid_sys_init->initialize_solution(s);
    _fluid_sys->solution->swap(base_sol);
    
    // create the nonlinear assembly object
    _frequency_domain_fluid_assembly = new MAST::FrequencyDomainLinearizedComplexAssembly;
    
    // solver for complex solution
    _complex_solver                  = new MAST::ComplexSolverBase;
    
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE STRUCTURAL DATA
    //////////////////////////////////////////////////////////////////////
    
    _n_divs_x                = infile("n_divs_x",             40);
    _n_divs_between_stiff    = infile("n_divs_between_stiff", 10);
    _n_stiff                 = infile("n_stiffeners",          3);
    
    _n_plate_elems           = _n_divs_x*(_n_stiff+1)*_n_divs_between_stiff;
    _n_elems_per_stiff       = _n_divs_x;
    _n_elems                 = _n_plate_elems + _n_stiff * _n_elems_per_stiff;
    
    
    // number of stations
    _n_dv_stations_x         = infile("n_stations", 4);
    
    // number of load steps
    _n_load_steps            = infile("n_load_steps", 20);
    
    
    // now setup the optimization data
    _n_eig                 =  infile("n_modes", 20);
    _n_vars                = _n_dv_stations_x + 2*_n_dv_stations_x * _n_stiff; // for thickness variable
    _n_eq                  = 0;
    _n_ineq                = _n_eig + 1 + _n_elems; // constraint that each eigenvalue > 0 flutter constraint + one element stress functional per elem
    _max_iters             = 1000;
    
    
    
    // length of domain
    _length        = infile("x_div_loc",   0., 2);
    _length       -= infile("x_div_loc",   0., 1);
    
    _width         = infile("y_div_loc",   0., 2);
    _width        -= infile("y_div_loc",   0., 1);
    
    // limit stress
    _stress_limit  = infile("max_stress", 4.00e8);
    
    // setup length for use in setup of flutter solver
    (*_b_ref) = _length;

    
    // initialize the dv vector data
    const Real
    th_l                   = infile("thickness_lower", 0.0001),
    th_u                   = infile("thickness_upper", 0.01),
    th                     = infile("thickness", 0.0015),
    dx                     = _length/(_n_stations-1);
    
    _dv_init.resize            (_n_vars);
    _dv_scaling.resize         (_n_vars);
    _dv_low.resize             (_n_vars);
    _problem_parameters.resize (_n_vars);
    
    // design variables for the thickness values
    for (unsigned int i=0; i<_n_vars; i++) {
        
        _dv_init[i]    =  infile("dv_init", th/th_u, i);
        _dv_low[i]     = th_l/th_u;
        _dv_scaling[i] =      th_u;
    }
    
    // create the mesh
    _structural_mesh       = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    MAST::StiffenedPanelMesh panel_mesh;
    panel_mesh.init(_n_stiff,
                    _n_divs_x,
                    _n_divs_between_stiff,
                    _length,
                    _width,
                    *_structural_mesh,
                    libMesh::TRI3,
                    true);
    
    
    // create the equation system
    _structural_eq_sys    = new  libMesh::EquationSystems(*_structural_mesh);
    
    // create the libmesh system
    _structural_sys       = &(_structural_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _structural_sys->set_eigenproblem_type(libMesh::GHEP);
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys_init  =
    new MAST::StructuralSystemInitialization(*_structural_sys,
                                             _structural_sys->name(),
                                             fetype);
    _structural_discipline = new MAST::StructuralDiscipline(*_structural_eq_sys);
    
    // create and add the boundary condition and loads
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 5;  // tz
    /*_dirichlet_left->init  (0, constrained_vars);
    _dirichlet_right->init (1, constrained_vars);
    _dirichlet_top->init   (2, constrained_vars);
    _dirichlet_bottom->init(3, constrained_vars);*/
    _dirichlet_left->init  (0, _structural_sys_init->vars());
    _dirichlet_right->init (1, _structural_sys_init->vars());
    _dirichlet_top->init   (2, _structural_sys_init->vars());
    _dirichlet_bottom->init(3, _structural_sys_init->vars());
    _structural_discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _structural_discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _structural_discipline->add_dirichlet_bc(2, *_dirichlet_top);
    _structural_discipline->add_dirichlet_bc(3, *_dirichlet_bottom);
    _structural_discipline->init_system_dirichlet_bc(*_structural_sys);
    
    // initialize the equation system
    _structural_eq_sys->init();
    
    // initialize the motion object
    _motion_function   = new MAST::FlexibleSurfaceMotion("small_disturbance_motion",
                                                         *_structural_sys_init);
    _slip_wall->add(*_motion_function);
    
    
    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _structural_sys->set_exchange_A_and_B(true);
    _structural_sys->set_n_requested_eigenvalues(_n_eig);
    
    
    // create the thickness variables
    _th_station_parameters_plate.resize(_n_vars);
    _th_station_functions_plate.resize(_n_vars);
    
    std::map<Real, MAST::FieldFunction<Real>*>
    th_station_vals,
    thy_station_vals,
    thz_station_vals;
    
    for (unsigned int i=0; i<_n_stations; i++) {
        std::ostringstream oss;
        oss << "h_" << i;
        
        // now we need a parameter that defines the thickness at the
        // specified station and a constant function that defines the
        // field function at that location.
        MAST::Parameter* h               =
        new MAST::Parameter(oss.str(), infile("thickness", 0.0015));
        
        MAST::ConstantFieldFunction* h_f =
        new MAST::ConstantFieldFunction("h", *h);
        
        // add this to the thickness map
        th_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                (i*dx, h_f));
        
        // add the function to the parameter set
        _th_station_parameters_plate[i]          = h;
        _th_station_functions_plate[i]           = h_f;
        
        // tell the assembly system about the sensitvity parameter
        _structural_discipline->add_parameter(*h);
        _problem_parameters[i] = h;
    }
    
    // now create the h function and give it to the property card
    _th_plate_f.reset(new MAST::MultilinearInterpolation("h", th_station_vals));
    
    
    // create the property functions and add them to the
    
    _rho             = new MAST::Parameter("rho",   infile("rho",    2.7e3));
    _E               = new MAST::Parameter("E",     infile(  "E",    72.e9));
    _nu              = new MAST::Parameter("nu",    infile( "nu",     0.33));
    _kappa           = new MAST::Parameter("kappa",    infile("kappa",  5./6.));
    _alpha           = new MAST::Parameter("alpha",infile("alpha",    2.5e-5));
    _zero            = new MAST::Parameter("zero",     0.);
    _temp            = new MAST::Parameter( "temperature",infile("temp", 60.));
    _p_cav           = new MAST::Parameter("p_cav",infile("p_cav",        0.));
    
    _rho_f           = new MAST::ConstantFieldFunction("rho",               *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",                   *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",                 *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",           *_kappa);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    _thzoff_stiff_f  = new MAST::ConstantFieldFunction("hz_off",           *_zero);
    _hoff_plate_f    = new MAST::ConstantFieldFunction("off",              *_zero);
    _temp_f          = new MAST::ConstantFieldFunction("temperature",      *_temp);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature",  *_zero);
    _p_cav_f         = new MAST::ConstantFieldFunction("pressure",        *_p_cav);

    // create the material property card
    _m_card          = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(  *_rho_f);
    _m_card->add(    *_E_f);
    _m_card->add(   *_nu_f);
    _m_card->add(*_kappa_f);
    _m_card->add(*_alpha_f);
    
    // create the element property card
    _p_card_plate      = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card_plate->add(*_th_plate_f);
    _p_card_plate->add(*_hoff_plate_f);
    
    // tell the section property about the material property
    _p_card_plate->set_material(*_m_card);
    //if (if_nonlin) _p_card_plate->set_strain(MAST::VON_KARMAN_STRAIN);

    _structural_discipline->set_property_for_subdomain(0, *_p_card_plate);

    
    // now add the property cards for each stiffener
    // element orientation
    libMesh::Point orientation;
    orientation(2) = 1.;
    
    // property card per stiffener
    _p_card_stiff.resize(_n_stiff);
    
    // thickness per stiffener station
    _thy_station_parameters_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thy_station_functions_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thz_station_parameters_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thz_station_functions_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thy_stiff_f.resize(_n_stiff);
    _thz_stiff_f.resize(_n_stiff);
    _hyoff_stiff_f.resize(_n_stiff);
    
    for (unsigned int i=0; i<_n_stiff; i++) {
        
        // this map is used to store the thickness parameter along length
        thy_station_vals.clear();
        thz_station_vals.clear();
        
        // first define the thickness station parameters and the thickness
        // field function
        for (unsigned int j=0; j<_n_dv_stations_x; j++) {
            std::ostringstream ossy, ossz;
            ossy << "h_y_" << j << "_stiff_" << i;
            ossz << "h_z_" << j << "_stiff_" << i;
            
            // now we need a parameter that defines the thickness at the
            // specified station and a constant function that defines the
            // field function at that location.
            MAST::Parameter
            *h_y  = new MAST::Parameter(ossy.str(), infile("thickness", 0.002)),
            *h_z  = new MAST::Parameter(ossz.str(), infile("thickness", 0.002));
            
            MAST::ConstantFieldFunction
            *h_y_f = new MAST::ConstantFieldFunction(ossy.str(), *h_y),
            *h_z_f = new MAST::ConstantFieldFunction(ossy.str(), *h_z);
            
            // add this to the thickness map
            thy_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                    (j*dx, h_y_f));
            thz_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                    (j*dx, h_z_f));
            
            // add the function to the parameter set
            _thy_station_parameters_stiff[i*_n_dv_stations_x+j]          = h_y;
            _thy_station_functions_stiff [i*_n_dv_stations_x+j]          = h_y_f;
            _thz_station_parameters_stiff[i*_n_dv_stations_x+j]          = h_z;
            _thz_station_functions_stiff [i*_n_dv_stations_x+j]          = h_z_f;
            
            
            // tell the assembly system about the sensitvity parameter
            _structural_discipline->add_parameter(*h_y);
            _structural_discipline->add_parameter(*h_z);
            _problem_parameters[(2*i+1)*_n_dv_stations_x+j] = h_y;
            _problem_parameters[(2*i+2)*_n_dv_stations_x+j] = h_z;
        }
        
        // now create the h_y function and give it to the property card
        _thy_stiff_f[i]   = new MAST::MultilinearInterpolation("hy", thy_station_vals);
        _thz_stiff_f[i]   = new MAST::MultilinearInterpolation("hz", thz_station_vals);
        _hyoff_stiff_f[i] = new MAST::SectionOffset("hy_off",
                                                    *_thy_stiff_f[i],
                                                    -1.);
        
        _p_card_stiff[i]  = new MAST::Solid1DSectionElementPropertyCard;
        
        
        // add the section properties to the card
        _p_card_stiff[i]->add(*_thy_stiff_f[i]);
        _p_card_stiff[i]->add(*_thz_stiff_f[i]);
        _p_card_stiff[i]->add(*_hyoff_stiff_f[i]);
        _p_card_stiff[i]->add(*_thzoff_stiff_f);
        _p_card_stiff[i]->y_vector() = orientation;
        
        // tell the section property about the material property
        _p_card_stiff[i]->set_material(*_m_card);
        //if (if_nonlin) _p_card_stiff[i]->set_strain(MAST::VON_KARMAN_STRAIN);
        
        _p_card_stiff[i]->init();
        
        // the domain ID of the stiffener is 1 plus the stiff number
        _structural_discipline->set_property_for_subdomain(i+1, *_p_card_stiff[i]);
    }
    
    // initialize the null space object and assign it to the structural
    // module
    _nsp = new MAST::StructuralNearNullVectorSpace;
    _structural_sys->nonlinear_solver->nearnullspace_object = _nsp;
    
    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _structural_mesh->local_elements_begin(),
    e_end   = _structural_mesh->local_elements_end();
    
    // points where stress is evaluated
    std::vector<libMesh::Point> pts;
    
    for ( ; e_it != e_end; e_it++) {
        
        pts.clear();
        if ((*e_it)->type() == libMesh::QUAD4 ||
            (*e_it)->type() == libMesh::QUAD8 ||
            (*e_it)->type() == libMesh::QUAD9) {
            
            pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
        }
        else if ((*e_it)->type() == libMesh::TRI3 ||
                 (*e_it)->type() == libMesh::TRI6) {
            
            pts.push_back(libMesh::Point(1./3., 1./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(1./3., 1./3.,-1.)); // lower skin
            pts.push_back(libMesh::Point(2./3., 1./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(2./3., 1./3.,-1.)); // lower skin
            pts.push_back(libMesh::Point(1./3., 2./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(1./3., 2./3.,-1.)); // lower skin
        }
        else if ((*e_it)->type() == libMesh::EDGE2 ||
                 (*e_it)->type() == libMesh::EDGE3) {
            
            pts.push_back(libMesh::Point(-1/sqrt(3), 1., 0.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3),-1., 0.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3), 1., 0.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3),-1., 0.)); // lower skin
        }
        else
            libmesh_assert(false); // should not get here
        
        
        MAST::StressStrainOutputBase * output = new MAST::StressStrainOutputBase;
        
        // tell the object to evaluate the data for this object only
        std::set<const libMesh::Elem*> e_set;
        e_set.insert(*e_it);
        output->set_elements_in_domain(e_set);
        output->set_points_for_evaluation(pts);
        output->set_volume_loads(_structural_discipline->volume_loads());
        _outputs.push_back(output);
        
        _structural_discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }
    
    // make sure that the number of output objects here is the same as the
    // number of local elems in the mesh
    libmesh_assert_equal_to(_outputs.size(), _structural_mesh->n_local_elem());

    
    // initialize the load
    _T_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    _T_load->add(*_temp_f);
    _T_load->add(*_ref_temp_f);
    _structural_discipline->add_volume_load(0, *_T_load);          // for the panel
    for (unsigned int i=0; i<_n_stiff; i++)
        _structural_discipline->add_volume_load(i+1, *_T_load);    // for the stiffeners
    
    // pressure load
    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _p_load->add(*_p_cav_f);
    _structural_discipline->add_volume_load(0, *_p_load);          // for the panel
    
    
    
    // pressure boundary condition for the beam
    _pressure    =  new MAST::BoundaryConditionBase(MAST::SMALL_DISTURBANCE_MOTION);
    _pressure->add(*_small_dist_pressure_function);
    _pressure->add(*_motion_function);
    _structural_discipline->add_volume_load(0, *_pressure);
    
    // modal assembly object for the natural modes eigen-analysis
    _modal_assembly = new MAST::StructuralModalEigenproblemAssembly;
    
    // nonlinear assembly object for the thermoelastic analysis
    _structural_nonlinear_assembly = new MAST::StructuralNonlinearAssembly;
    
    // create the function to calculate weight
    _weight = new MAST::StiffenedPlateWeight(*_structural_discipline);
    
    _flutter_solver  = new MAST::UGFlutterSolver;
    
    ////////////////////////////////////////////////////////////
    // STRUCTURAL MODAL EIGENSOLUTION
    ////////////////////////////////////////////////////////////
    
    libMesh::out << "DVs used for basis generation..." << std::endl;
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "th     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << (*_problem_parameters[i])() << std::endl;

    // create the nonlinear assembly object
    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
    
    _modal_assembly->attach_discipline_and_system(*_structural_discipline,
                                                  *_structural_sys_init);
    
    _structural_sys->eigenproblem_solve();
    _modal_assembly->clear_discipline_and_system();
    
    // now set the nonlinear strain, if that has been requested
    if (if_nonlin) {
     
        _p_card_plate->set_strain(MAST::VON_KARMAN_STRAIN);
        for (unsigned int i=0; i<_n_stiff; i++)
            _p_card_stiff[i]->set_strain(MAST::VON_KARMAN_STRAIN);
    }
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_structural_sys->get_n_converged_eigenvalues(),
                     _structural_sys->get_n_requested_eigenvalues());
    
    if (_basis.size() > 0)
        libmesh_assert(_basis.size() == nconv);
    else {
        _basis.resize(nconv);
        for (unsigned int i=0; i<_basis.size(); i++)
            _basis[i] = NULL;
    }
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // create a vector to store the basis
        if (_basis[i] == NULL)
            _basis[i] = _structural_sys->solution->zero_clone().release();
        
        std::ostringstream file_name;
        
        // get the eigenvalue and eigenvector
        Real
        re = 0.,
        im = 0.;
        _structural_sys->get_eigenpair(i, re, im, *_basis[i]);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
    }


    //////////////////////////////////////////////////////////////////////
    //    CALCULATE AND STORE THE GAFs
    //////////////////////////////////////////////////////////////////////
    // initialize the GAF interpolation assembly object
    _gaf_database = new  MAST::GAFDatabase(_basis.size());
    bool
    calculate_gafs = infile("calculate_gafs", true);
    
    if (calculate_gafs) {
        
        // now setup the assembly object
        _frequency_domain_fluid_assembly->attach_discipline_and_system(*_fluid_discipline,
                                                                       *_complex_solver,
                                                                       *_fluid_sys_init);
        _frequency_domain_fluid_assembly->set_base_solution(base_sol);
        _frequency_domain_fluid_assembly->set_frequency_function(*_freq_function);

        
        _gaf_database->attach_discipline_and_system(*_structural_discipline,
                                                    *_structural_sys_init);
        
        
        _gaf_database->init(*_freq_function,
                            _complex_solver,               // fluid complex solver
                            _small_dist_pressure_function,
                            _motion_function);
        
        _gaf_database->set_evaluate_mode(true);
        
        libMesh::out
        << "Building GAF database..." << std::endl;
        
        // now iterate over the reduced frequencies and calculate the GAF matrices
        for (unsigned int i=0; i<=_n_k_divs; i++) {
            
            Real
            kval = _k_upper + (_k_lower-_k_upper)*(1.*i)/(1.*_n_k_divs);
            
            libMesh::out << " ***********   kr = " << kval
            << "  ***********" << std::endl;
            
            
            // initialize reduced frequency
            (*_omega) = kval;
            
            // first the GAF values, then the sensitivity values
            {
                ComplexMatrixX&
                mat = _gaf_database->add_kr_mat(kval,
                                                ComplexMatrixX::Zero(_basis.size(),
                                                                     _basis.size()),
                                                false);
                
                _gaf_database->assemble_generalized_aerodynamic_force_matrix(_basis, mat);
            }
            
            // now the sensitivity
            {
                ComplexMatrixX&
                mat = _gaf_database->add_kr_mat(kval,
                                                ComplexMatrixX::Zero(_basis.size(),
                                                                     _basis.size()),
                                                true);
                
                _gaf_database->assemble_generalized_aerodynamic_force_matrix(_basis,
                                                                             mat,
                                                                             _omega);
            }
        }
        
        _gaf_database->clear_discipline_and_system();
        _frequency_domain_fluid_assembly->clear_discipline_and_system();
        _gaf_database->set_evaluate_mode(false);
        _gaf_database->write_gaf_file("gaf_database.txt", _basis);
    }
    else
        _gaf_database->read_gaf_file("gaf_database.txt", _basis);
    
    _initialized = true;
}




MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::~StiffenedPlateThermallyStressedFSIFlutterSizingOptimization() {

    if (!_initialized)
        return;

    
    delete _fluid_eq_sys;
    delete _fluid_mesh;
    
    delete _fluid_discipline;
    delete _fluid_sys_init;
    
    delete _frequency_domain_fluid_assembly;
    delete _complex_solver;

    delete _far_field;
    delete _symm_wall;
    delete _slip_wall;
    
    delete _flight_cond;
    
    delete _omega;
    delete _velocity;
    delete _b_ref;
    
    delete _omega_f;
    delete _velocity_f;
    delete _b_ref_f;
    
    delete _freq_function;
    
    delete _small_dist_pressure_function;
    
    delete _motion_function;
    delete _pressure;
    
    delete _T_load;
    delete _p_load;
    
    delete _structural_eq_sys;
    delete _structural_mesh;
    
    delete _modal_assembly;
    delete _structural_nonlinear_assembly;
    
    delete _structural_discipline;
    delete _structural_sys_init;
    
    delete _m_card;
    delete _p_card_plate;
    for (unsigned int i=0; i<_n_stiff; i++) delete _p_card_stiff[i];
   
    delete _dirichlet_left;
    delete _dirichlet_right;
    delete _dirichlet_top;
    delete _dirichlet_bottom;
    
    delete _rho_f;
    delete _E_f;
    delete _p_cav_f;
    delete _nu_f;
    delete _ref_temp_f;
    delete _temp_f;
    delete _kappa_f;
    delete _alpha_f;
    delete _hoff_plate_f;
    
    delete _thzoff_stiff_f;
    for (unsigned int i=0; i<_n_stiff; i++) delete _hyoff_stiff_f[i];
    for (unsigned int i=0; i<_n_stiff; i++) delete _thy_stiff_f[i];
    for (unsigned int i=0; i<_n_stiff; i++) delete _thz_stiff_f[i];
    
    
    delete _rho;
    delete _E;
    delete _p_cav;
    delete _temp;
    delete _nu;
    delete _zero;
    delete _kappa;
    delete _alpha;
    
    
    // delete the basis vectors
    if (_basis.size())
        for (unsigned int i=0; i<_basis.size(); i++)
            if (_basis[i]) delete _basis[i];
    
    // delete the h_y station functions
    {
        std::vector<MAST::ConstantFieldFunction*>::iterator
        it  = _th_station_functions_plate.begin(),
        end = _th_station_functions_plate.end();
        for (; it != end; it++)  delete *it;
    }
    
    {
        std::vector<MAST::ConstantFieldFunction*>::iterator
        it  = _thy_station_functions_stiff.begin(),
        end = _thy_station_functions_stiff.end();
        for (; it != end; it++)  delete *it;
    }
    
    {
        std::vector<MAST::ConstantFieldFunction*>::iterator
        it  = _thz_station_functions_stiff.begin(),
        end = _thz_station_functions_stiff.end();
        for (; it != end; it++)  delete *it;
    }
    
    // delete the h_y station parameters
    {
        std::vector<MAST::Parameter*>::iterator
        it  = _th_station_parameters_plate.begin(),
        end = _th_station_parameters_plate.end();
        for (; it != end; it++)  delete *it;
    }
    
    {
        std::vector<MAST::Parameter*>::iterator
        it  = _thy_station_parameters_stiff.begin(),
        end = _thy_station_parameters_stiff.end();
        for (; it != end; it++)  delete *it;
    }
    
    {
        std::vector<MAST::Parameter*>::iterator
        it  = _thz_station_parameters_stiff.begin(),
        end = _thz_station_parameters_stiff.end();
        for (; it != end; it++)  delete *it;
    }
    
    
    delete _weight;
    
    delete _flutter_solver;
    delete _gaf_database;
    delete _augment_send_list_obj;
}



void
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::
init_dvar(std::vector<Real>& x,
          std::vector<Real>& xmin,
          std::vector<Real>& xmax) {
    
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}


// define the class to provide the interface for nonlinear steady-state
// solution that can be used by the flutter solver to solve for different
// flight velocities
namespace MAST {
    
    class StiffenedPlateThermallyStressedFSISteadySolverInterface:
    public MAST::FlutterSolverBase::SteadySolver {
        
    public:
        StiffenedPlateThermallyStressedFSISteadySolverInterface
        (MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization& obj,
         bool if_output,
         bool if_clear_vector_on_exit,
         unsigned int n_steps):
        MAST::FlutterSolverBase::SteadySolver(),
        _obj(obj),
        _if_write_output(if_output),
        _n_steps(n_steps),
        _if_clear_vector_on_exit(if_clear_vector_on_exit) {
            
            _obj._structural_sys->add_vector("base_solution");
        }
        
        virtual ~StiffenedPlateThermallyStressedFSISteadySolverInterface() {
            
            if (_if_clear_vector_on_exit)
                _obj._structural_sys->remove_vector("base_solution");
        }
        
        
        /*!
         *  solves for the steady state solution, and @returns
         *  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solve() {
            
            libMesh::NumericVector<Real>&
            sol = _obj._structural_sys->get_vector("base_solution");
            *_obj._structural_sys->solution = sol;
            
            // now do the solve
            libmesh_assert(_obj._initialized);
            
            bool if_vk = (_obj._p_card_plate->strain_type() == MAST::VON_KARMAN_STRAIN);
            
            ///////////////////////////////////////////////////////////////
            // first, solve the quasi-steady problem
            ///////////////////////////////////////////////////////////////
            // set the number of load steps
            unsigned int
            n_steps = 1;
            if (if_vk) n_steps = _n_steps;
            
            Real
            T0      = (*_obj._temp)(),
            p0      = (*_obj._p_cav)();
            
            // zero the solution before solving
            _obj.clear_stresss();
            
            _obj._structural_nonlinear_assembly->attach_discipline_and_system
            (*_obj._structural_discipline,
             *_obj._structural_sys_init);
            libMesh::ExodusII_IO writer(_obj._structural_sys->get_mesh());
            
            // now iterate over the load steps
            for (unsigned int i=0; i<n_steps; i++) {
                
                // modify the load
                (*_obj._temp)()      =  T0*(i+1.)/(1.*n_steps);
                (*_obj._p_cav)()     =  p0*(i+1.)/(1.*n_steps);
                
                libMesh::out
                << "Load step: " << i
                << "  : T = " << (*_obj._temp)()
                << "  : p = " << (*_obj._p_cav)()
                << std::endl;
                
                _obj._structural_sys->solve();
                _obj._structural_discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
                writer.write_timestep("output_all_steps.exo",
                                      *_obj._structural_eq_sys,
                                      i+1,
                                      (i+1)/(1.*n_steps));
                
            }
            
            // copy the solution to the base solution vector
            sol = *_obj._structural_sys->solution;
            
            _obj._structural_nonlinear_assembly->clear_discipline_and_system();
            
            return sol;
        }
        
        
        /*!
         *   sets the number of steps to be used for nonlinaer steady analysis.
         *   The default is 25 for noninear and 1 for linear.
         */
        void set_n_load_steps(unsigned int n) {
            _n_steps = n;
        }
        
        
        /*!
         * @returns  a non-const-reference to the solution.
         */
        virtual libMesh::NumericVector<Real>&
        solution() { return _obj._structural_sys->get_vector("base_solution"); }
        
        
        
        /*!
         * @returns  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solution() const { return _obj._structural_sys->get_vector("base_solution"); }
        
        
    protected:
        
        /*!
         *   pointer to the object that hold all the solution data
         */
        MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization& _obj;
        
        /*!
         *   flag to toggle output
         */
        bool _if_write_output;
        
        /*!
         *   number of nonliear load increment steps
         */
        unsigned int _n_steps;
        
        /*!
         *   deletes the solution vector from system when the class is
         *   destructed unless this flag is false.
         */
        bool _if_clear_vector_on_exit;
    };
}



void
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::
evaluate(const std::vector<Real>& dvars,
         Real& obj,
         bool eval_obj_grad,
         std::vector<Real>& obj_grad,
         std::vector<Real>& fvals,
         std::vector<bool>& eval_grads,
         std::vector<Real>& grads) {
    
    libmesh_assert(_initialized);
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // while this function tries to ensure that the gradient and function values
    // returned to the optimizer on all ranks is the same, we may run the
    // possibility of slight differences in the input DVs due to unforeseen
    // reasons. Hence, we will borrow the DV values from rank 0, just to
    // make sure that all processors are evaluating the functions at the
    // same exact values.
    std::vector<Real>
    my_dvars(dvars.begin(), dvars.end());
    __init->comm().broadcast(my_dvars);
    
    
    // set the parameter values equal to the DV value
    // first the plate thickness values
    for (unsigned int i=0; i<_n_vars; i++)
        (*_problem_parameters[i]) = my_dvars[i]*_dv_scaling[i];
    
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    // zero the function evaluations
    std::fill(fvals.begin(), fvals.end(), 0.);
    
    
    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "th     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << (*_problem_parameters[i])() << std::endl;

    // the optimization problem is defined as
    // min weight, subject to constraints on displacement and stresses
    Real
    wt      = 0.,
    pval    = 2.;
    
    bool
    if_write_output = true;

    
    // calculate weight
    (*_weight)(pt, 0., wt);


    ///////////////////////////////////////////////////////////////////
    // THERMOELASTIC SOLUTION
    ///////////////////////////////////////////////////////////////////
    _structural_sys->solution->zero();
    this->clear_stresss();

    MAST::StiffenedPlateThermallyStressedFSISteadySolverInterface
    steady_solve(*this,
                 if_write_output,
                 false,
                 _n_load_steps);
    
    steady_solve.solve();

    
    
    ///////////////////////////////////////////////////////////////////
    // MODAL SOLUTION
    ///////////////////////////////////////////////////////////////////
    // modal analysis is about the thermoelastic base state, and is
    // used to ensure that the natural frequencies are positive about this
    // state
    _modal_assembly->attach_discipline_and_system(*_structural_discipline,
                                                  *_structural_sys_init);
    _modal_assembly->set_base_solution(steady_solve.solution());
    _structural_sys->eigenproblem_solve();
    _modal_assembly->clear_discipline_and_system();
    
    unsigned int
    nconv = std::min(_structural_sys->get_n_converged_eigenvalues(),
                     _structural_sys->get_n_requested_eigenvalues());
    
    
    // vector of eigenvalues
    std::vector<Real> eig_vals(nconv);
    
    bool if_all_eig_positive = true;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    mode_vec(_structural_sys->solution->zero_clone().release());
    
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        _structural_sys->get_eigenpair(i, re, im, *mode_vec);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        eig_vals[i]          = re;
        if_all_eig_positive  = (if_all_eig_positive && (re>0.))?true:false;
        
        if (if_write_output) {
            
            std::ostringstream file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << i
            << ".exo";
            
            libMesh::out
            << "Writing mode " << i << " to : "
            << file_name.str() << std::endl;
            
            // copy the solution for output
            (*_structural_sys->solution) = *mode_vec;
            
            // We write the file in the ExodusII format.
            std::set<std::string> nm;
            nm.insert(_structural_sys->name());
            libMesh::ExodusII_IO(*_structural_mesh).write_equation_systems
            (file_name.str(),
             *_structural_eq_sys,
             &nm);
        }
    }

    

    ///////////////////////////////////////////////////////////////////
    // FLUTTER SOLUTION
    ///////////////////////////////////////////////////////////////////
    
    std::pair<bool, MAST::FlutterRootBase*> sol(false, nullptr);
    
    if (if_all_eig_positive) {
        
        // clear flutter solver and set the output file
        _flutter_solver->clear();
        
        std::ostringstream oss;
        oss << "flutter_output_" << __init->comm().rank() << ".txt";
        _flutter_solver->set_output_file(oss.str());
        
        _gaf_database->attach_discipline_and_system(*_structural_discipline,
                                                    *_structural_sys_init);
        
        libMesh::NumericVector<Real>&
        base_sol = _fluid_sys->get_vector("fluid_base_solution");
        _frequency_domain_fluid_assembly->attach_discipline_and_system(*_fluid_discipline,
                                                                       *_complex_solver,
                                                                       *_fluid_sys_init);
        _frequency_domain_fluid_assembly->set_base_solution(base_sol);
        _frequency_domain_fluid_assembly->set_frequency_function(*_freq_function);
        
        
        _gaf_database->init(*_freq_function,
                            _complex_solver,                       // fluid complex solver
                            _small_dist_pressure_function,
                            _motion_function);
        _flutter_solver->attach_assembly(*_gaf_database);
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
        Real
        tol                 = 1.0e-6;
        unsigned int
        max_bisection_iters = 10;
        
        sol = _flutter_solver->find_critical_root(tol, max_bisection_iters);
        
        
        _flutter_solver->print_sorted_roots();
        _gaf_database->clear_discipline_and_system();
        _flutter_solver->clear_assembly_object();
        _frequency_domain_fluid_assembly->clear_discipline_and_system();
    }
    
    
    _structural_nonlinear_assembly->attach_discipline_and_system(*_structural_discipline,
                                                                 *_structural_sys_init);
    _structural_nonlinear_assembly->calculate_outputs(steady_solve.solution());
    _structural_nonlinear_assembly->clear_discipline_and_system();

    
    
    //////////////////////////////////////////////////////////////////////
    // get the objective and constraints
    //////////////////////////////////////////////////////////////////////
    
    // set the function and objective values
    obj = wt;
    
    // parallel sum of the weight
    this->comm().sum(obj);
    
    
    // set the eigenvalue constraints  -eig <= 0. scale
    // by an arbitrary 1/1.e7 factor
    for (unsigned int i=0; i<nconv; i++)
        fvals[i] = -eig_vals[i]/1.e7;
    
    
    // we need to make sure that the stress functional in a parallel environment
    // correspond to unique spots in the constraint function vector. Below,
    // the output functionals are created for each local element. These will
    // be mapped to unique spots by identifying the number of elements on each
    // subdomain
    std::vector<unsigned int>
    beginning_elem_id(__init->comm().size(), 0);
    for (unsigned int i=1; i<__init->comm().size(); i++)
        beginning_elem_id[i] = _structural_mesh->n_elem_on_proc(i-1);
    
    // now use this info to identify the beginning elem id
    for (unsigned int i=2; i<__init->comm().size(); i++)
        beginning_elem_id[i] += beginning_elem_id[i-1];
    
    const unsigned int
    my_id0 = beginning_elem_id[__init->comm().rank()];
    
    // copy the element von Mises stress values as the functions
    for (unsigned int i=0; i<_outputs.size(); i++)
        fvals[_n_eig+1+i+my_id0] =  -1. +
        _outputs[i]->von_Mises_p_norm_functional_for_all_elems(pval)/_stress_limit;
    
    // now sum the value of the stress constraints, so that all ranks
    // have the same values
    for (unsigned int i=0; i<_n_elems; i++)
        this->comm().sum(fvals[_n_eig+1+i]);
    
    // copy the flutter velocity to the constraint vector
    //     Vf        >= V0
    // or, V0/Vf     <= 1
    // or, V0/Vf - 1 <= 0
    //
    if (!if_all_eig_positive)
        fvals[_n_eig+0]  =  100.;
    else if (sol.second)
        fvals[_n_eig+0]  =  _V0_flutter/sol.second->V - 1.;
    else
        fvals[_n_eig+0]  =  -100.;

    //////////////////////////////////////////////////////////////////
    //   evaluate sensitivity if needed
    //////////////////////////////////////////////////////////////////
    
    // sensitivity of the objective function
    if (eval_obj_grad) {
        
        Real w_sens = 0.;
        
        // set gradient of weight
        for (unsigned int i=0; i<_n_vars; i++) {
            
            _weight->derivative(MAST::PARTIAL_DERIVATIVE,
                                *_problem_parameters[i],
                                pt,
                                0.,
                                w_sens);
            obj_grad[i] = w_sens*_dv_scaling[i];
        }
        
        // parallel sum
        this->comm().sum(obj_grad);
    }

    
    // now check if the sensitivity of constraint function is requested
    bool if_sens = false;
    
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_sens = (if_sens || eval_grads[i]);
    
    if (if_sens) {
        
        // first initialize the gradient vector to zero
        std::fill(grads.begin(), grads.end(), 0.);

        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        
        // we are going to choose to use one parametric sensitivity at a time
        for (unsigned int i=0; i<_n_vars; i++) {
            
            libMesh::ParameterVector params;
            
            params.resize(1);
            params[0]  = _problem_parameters[i]->ptr();
            
            // copy the solution to be used for sensitivity
            // If a flutter solution was found, then this depends on velocity.
            // Otherwise, it is independent of velocity
            *_structural_sys->solution = steady_solve.solution();
            
            // iterate over each dv and calculate the sensitivity
            libMesh::NumericVector<Real>&
            dXdp = _structural_sys->add_sensitivity_solution(0);
            dXdp.zero();
            this->clear_stresss();
            
            // sensitivity analysis
            _structural_nonlinear_assembly->attach_discipline_and_system
            (*_structural_discipline,
             *_structural_sys_init);
            _structural_sys->sensitivity_solve(params);
            
            // evaluate sensitivity of the outputs
            _structural_nonlinear_assembly->calculate_output_sensitivity
            (params,
             true,    // true for total sensitivity
             *(_structural_sys->solution));
            
            _structural_nonlinear_assembly->clear_discipline_and_system();
            
            // copy the sensitivity values in the output. This accounts for the
            // sensitivity of state wrt parameter. However, if a flutter root
            // was found, the state depends on velocity, which depends on the
            // parameter. Hence, the total sensitivity of stress constraint
            // would need to include the latter component, which was added
            // above.
            for (unsigned int j=0; j<_outputs.size(); j++)
                grads[(i*_n_ineq) + (j+_n_eig+1+my_id0)] =
                _dv_scaling[i]/_stress_limit *
                _outputs[j]->von_Mises_p_norm_functional_sensitivity_for_all_elems
                (pval, _problem_parameters[i]);

            
            
            // calculate sensitivity only if the flutter root was found
            // else, set it to zero
            if (sol.second) {
                
                _gaf_database->attach_discipline_and_system(*_structural_discipline,
                                                            *_structural_sys_init);
                
                libMesh::NumericVector<Real>& base_sol =
                _fluid_sys->add_vector("fluid_base_solution");
                _frequency_domain_fluid_assembly->attach_discipline_and_system(*_fluid_discipline,
                                                                               *_complex_solver,
                                                                               *_fluid_sys_init);
                _frequency_domain_fluid_assembly->set_base_solution(base_sol);
                _frequency_domain_fluid_assembly->set_frequency_function(*_freq_function);
                
                _gaf_database->init(*_freq_function,
                                    _complex_solver,                       // fluid complex solver
                                    _small_dist_pressure_function,
                                    _motion_function);
                _flutter_solver->attach_assembly(*_gaf_database);
                _flutter_solver->initialize(*_omega,
                                            *_b_ref,
                                            _flight_cond->rho(),
                                            _k_lower,         // lower kr
                                            _k_upper,         // upper kr
                                            _n_k_divs,        // number of divisions
                                            _basis);          // basis vectors

                _flutter_solver->calculate_sensitivity(*sol.second,
                                                       params,
                                                       0,
                                                       &dXdp);

                _gaf_database->clear_discipline_and_system();
                _flutter_solver->clear_assembly_object();
                _frequency_domain_fluid_assembly->clear_discipline_and_system();

                // copy the sensitivity values in the output
                grads[(i*_n_ineq) + (_n_eig+0)] = -_dv_scaling[i] *
                _V0_flutter / pow(sol.second->V,2) * sol.second->V_sens;
            }
            else
                grads[(i*_n_ineq) + (_n_eig+0)] = 0.;

            // now, sum the sensitivity of the stress function gradients
            // so that all processors have the same values
            for (unsigned int j=0; j<_n_elems; j++)
                __init->comm().sum(grads[(i*_n_ineq) + (j+_n_eig+1)]);
            
            
            // calculate the sensitivity of the eigenvalues
            std::vector<Real> eig_sens(nconv);
            _modal_assembly->set_base_solution(steady_solve.solution());
            _modal_assembly->set_base_solution(dXdp, true);
            _modal_assembly->attach_discipline_and_system(*_structural_discipline,
                                                          *_structural_sys_init);
            // this should not be necessary, but currently the eigenproblem sensitivity
            // depends on availability of matrices before sensitivity
            _structural_sys->assemble_eigensystem();
            _structural_sys->eigenproblem_sensitivity_solve(params, eig_sens);
            _modal_assembly->clear_discipline_and_system();
            
            for (unsigned int j=0; j<nconv; j++)
                grads[(i*_n_ineq) + j] = -_dv_scaling[i]*eig_sens[j]/1.e7;
        }
        
        // tell all ranks about the gradients
        this->comm().broadcast(grads);
    }
}



void
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::clear_stresss() {
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}





void
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::
output(unsigned int iter,
       const std::vector<Real>& x,
       Real obj,
       const std::vector<Real>& fval,
       bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
}


MAST::FunctionEvaluation::funobj
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::get_objective_evaluation_function() {
    
    return stiffened_plate_euler_flutter_optim_obj;
}



MAST::FunctionEvaluation::funcon
MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization::get_constraint_evaluation_function() {
    
    return stiffened_plate_euler_flutter_optim_con;
}

