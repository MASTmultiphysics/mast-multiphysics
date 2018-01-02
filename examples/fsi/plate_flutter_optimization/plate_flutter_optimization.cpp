/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include <map>

// MAST includes
#include "examples/fsi/plate_flutter_optimization/plate_flutter_optimization.h"
#include "examples/fsi/base/gaf_database.h"
#include "optimization/optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "base/eigenproblem_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/stress_output_base.h"
#include "base/nonlinear_system.h"
#include "examples/fluid/meshing/mesh_initializer.h"
#include "examples/fluid/meshing/panel_mesh_3D.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/flight_condition.h"
#include "fluid/pressure_function.h"
#include "fluid/frequency_domain_pressure_function.h"
#include "base/complex_assembly_base.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "solver/complex_solver_base.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/ug_flutter_solver.h"
#include "base/complex_mesh_field_function.h"
#include "elasticity/complex_normal_rotation_mesh_function.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "aeroelasticity/flutter_root_base.h"
#include "examples/base/augment_ghost_elem_send_list.h"
#include "solver/slepc_eigen_solver.h"


// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"


extern
libMesh::LibMeshInit     *__init;
extern
MAST::FunctionEvaluation *__my_func_eval;



void
plate_euler_flutter_optim_obj(int*    mode,
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
plate_euler_flutter_optim_con(int*    mode,
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




MAST::PlateFSIFlutterSizingOptimization::
PlateFSIFlutterSizingOptimization(const libMesh::Parallel::Communicator& comm):
MAST::FunctionEvaluation(comm),
_initialized                            (false),
_length                                 (0.),
_width                                  (0.),
_k_lower                                (0.),
_k_upper                                (0.),
_V0_flutter                             (0.),
_n_elems                                (0),
_n_stations                             (0),
_n_k_divs                               (0.),
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
_frequency_domain_fluid_assembly        (nullptr),
_frequency_domain_elem_ops              (nullptr),
_complex_solver                         (nullptr),
_modal_assembly                         (nullptr),
_modal_elem_ops                         (nullptr),
_flight_cond                            (nullptr),
_far_field                              (nullptr),
_symm_wall                              (nullptr),
_slip_wall                              (nullptr),
_pressure                               (nullptr),
_displ                                  (nullptr),
_normal_rot                             (nullptr),
_pressure_function                      (nullptr),
_freq_domain_pressure_function          (nullptr),
_omega                                  (nullptr),
_velocity                               (nullptr),
_b_ref                                  (nullptr),
_omega_f                                (nullptr),
_velocity_f                             (nullptr),
_b_ref_f                                (nullptr),
_flutter_solver                         (nullptr),
_gaf_database                           (nullptr),
_E                                      (nullptr),
_nu                                     (nullptr),
_kappa                                  (nullptr),
_rho                                    (nullptr),
_press                                  (nullptr),
_zero                                   (nullptr),
_E_f                                    (nullptr),
_nu_f                                   (nullptr),
_kappa_f                                (nullptr),
_rho_f                                  (nullptr),
_press_f                                (nullptr),
_weight                                 (nullptr),
_m_card                                 (nullptr),
_p_card                                 (nullptr),
_dirichlet_left                         (nullptr),
_dirichlet_right                        (nullptr),
_dirichlet_bottom                       (nullptr),
_dirichlet_top                          (nullptr),
_augment_send_list_obj                  (nullptr)
 { }


void
MAST::PlateFSIFlutterSizingOptimization::init(GetPot &infile,
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
    
    std::unique_ptr<MAST::MeshInitializer::CoordinateDivisions>
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
    
    _pressure_function =
    new MAST::PressureFunction(*_fluid_sys_init, *_flight_cond);
    _freq_domain_pressure_function =
    new MAST::FrequencyDomainPressureFunction(*_fluid_sys_init, *_flight_cond);
    
    _pressure_function->set_calculate_cp(true);
    _freq_domain_pressure_function->set_calculate_cp(true);

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
    _frequency_domain_fluid_assembly = new MAST::ComplexAssemblyBase;
    _frequency_domain_elem_ops       =
    new MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations;
    
    // solver for complex solution
    _complex_solver                  = new MAST::ComplexSolverBase;
    
    // now setup the assembly object
    _frequency_domain_fluid_assembly->attach_discipline_and_system(*_frequency_domain_elem_ops,
                                                                   *_fluid_discipline,
                                                                   *_complex_solver,
                                                                   *_fluid_sys_init);
    _frequency_domain_fluid_assembly->set_base_solution(base_sol);
    _frequency_domain_elem_ops->set_frequency_function(*_freq_function);

    
    
    //////////////////////////////////////////////////////////////////////
    //    SETUP THE STRUCTURAL DATA
    //////////////////////////////////////////////////////////////////////
    
    x_div_loc.resize     (2);
    y_div_loc.resize     (2);
    x_relative_dx.resize (2);
    y_relative_dx.resize (2);
    x_divs.resize        (1);
    y_divs.resize        (1);
    divs.resize          (2);
    
    x_coord_divs.reset   (new MeshInitializer::CoordinateDivisions);
    y_coord_divs.reset   (new MeshInitializer::CoordinateDivisions);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<2; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div+1);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div+1);
    }
    x_divs[0]       = infile( "x_div_nelem", 0, 1);
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(1, x_div_loc, x_relative_dx, x_divs);
    
    // now read in the values: y-coord
    for (unsigned int i_div=0; i_div<2; i_div++) {
        
        y_div_loc[i_div]        = infile("y_div_loc",   0., i_div+1);
        y_relative_dx[i_div]    = infile( "y_rel_dx",   0., i_div+1);
    }
    y_divs[0]       = infile( "y_div_nelem", 0, 1);
    
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(1, y_div_loc, y_relative_dx, y_divs);

    
    // setup length for use in setup of flutter solver
    _length = x_div_loc[1]-x_div_loc[0];
    (*_b_ref) = _length;
    
    // initialize the dv vector data
    const Real
    th_l                   = infile("thickness_lower", 0.0001),
    th_u                   = infile("thickness_upper", 0.01),
    th                     = infile("thickness", 0.0015),
    dx                     = _length/(_n_stations-1);
    
    _dv_init.resize    (_n_vars);
    _dv_scaling.resize (_n_vars);
    _dv_low.resize     (_n_vars);
    
    // design variables for the thickness values
    for (unsigned int i=0; i<_n_vars; i++) {
        
        _dv_init[i]    =  infile("dv_init", th/th_u, i);
        _dv_low[i]     = th_l/th_u;
        _dv_scaling[i] =      th_u;
    }
    
    // create the mesh
    _structural_mesh       = new libMesh::SerialMesh(__init->comm());
    
    MeshInitializer().init(divs, *_structural_mesh, libMesh::QUAD4);
    
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
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 5;  // tz
    _dirichlet_left->init  (0, constrained_vars);
    _dirichlet_right->init (1, constrained_vars);
    _dirichlet_top->init   (2, constrained_vars);
    _dirichlet_bottom->init(3, constrained_vars);
    _structural_discipline->add_dirichlet_bc(0, *_dirichlet_left);
    _structural_discipline->add_dirichlet_bc(1, *_dirichlet_right);
    _structural_discipline->add_dirichlet_bc(2, *_dirichlet_top);
    _structural_discipline->add_dirichlet_bc(3, *_dirichlet_bottom);
    _structural_discipline->init_system_dirichlet_bc(*_structural_sys);
    
    // initialize the equation system
    _structural_eq_sys->init();
    
    // initialize the motion object
    _displ        = new MAST::ComplexMeshFieldFunction(*_structural_sys_init,
                                                       "frequency_domain_displacement");
    _normal_rot   = new MAST::ComplexNormalRotationMeshFunction("frequency_domain_normal_rotation",
                                                                *_displ);
    _slip_wall->add(*_displ);
    _slip_wall->add(*_normal_rot);
    
    
    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _structural_sys->set_exchange_A_and_B(true);
    _structural_sys->set_n_requested_eigenvalues(infile("n_modes", 16));
    
    
    // create the thickness variables
    _th_station_parameters.resize(_n_vars);
    _th_station_functions.resize(_n_vars);
    
    std::map<Real, MAST::FieldFunction<Real>*> th_station_vals;
    
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
        _th_station_parameters[i]          = h;
        _th_station_functions[i]           = h_f;
        
        // tell the assembly system about the sensitvity parameter
        _structural_discipline->add_parameter(*h);
    }
    
    // now create the h function and give it to the property card
    _th_f.reset(new MAST::MultilinearInterpolation("h", th_station_vals));
    
    
    // create the property functions and add them to the
    
    _rho             = new MAST::Parameter("rho",   infile("rho",    2.7e3));
    _E               = new MAST::Parameter("E",     infile(  "E",    72.e9));
    _nu              = new MAST::Parameter("nu",    infile( "nu",     0.33));
    _kappa           = new MAST::Parameter("kappa",    infile("kappa",  5./6.));
    _zero            = new MAST::Parameter("zero",     0.);
    
    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",     *_kappa);
    _hoff_f          = new MAST::ConstantFieldFunction("off",        *_zero);
    
    // create the material property card
    _m_card          = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(*_rho_f);
    _m_card->add(*_E_f);
    _m_card->add(*_nu_f);
    _m_card->add(*_kappa_f);
    
    // create the element property card
    _p_card          = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    
    _structural_discipline->set_property_for_subdomain(0, *_p_card);
    
    
    
    // pressure boundary condition for the beam
    _pressure    =  new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _pressure->add(*_pressure_function);
    _pressure->add(*_freq_domain_pressure_function);
    _structural_discipline->add_volume_load(0, *_pressure);
    
    // modal assembly object for the natural modes eigen-analysis
    _modal_assembly = new MAST::EigenproblemAssembly;
    _modal_elem_ops = new MAST::StructuralModalEigenproblemAssemblyElemOperations;
    
    
    // create the function to calculate weight
    _weight = new MAST::PlateWeight(*_structural_discipline);
    
    _flutter_solver  = new MAST::UGFlutterSolver;
    
    ////////////////////////////////////////////////////////////
    // STRUCTURAL MODAL EIGENSOLUTION
    ////////////////////////////////////////////////////////////
    
    // create the nonlinear assembly object
    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
    
    _modal_assembly->attach_discipline_and_system(*_modal_elem_ops,
                                                  *_structural_discipline,
                                                  *_structural_sys_init);
    
    
    _structural_sys->eigenproblem_solve();
    _modal_assembly->clear_discipline_and_system();
    
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
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // create a vector to store the basis
        if (_basis[i] == nullptr)
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

    _gaf_database->attach_discipline_and_system(*_gaf_database,
                                                *_structural_discipline,
                                                *_structural_sys_init);
    
    
    _gaf_database->init(_freq_function,
                        _complex_solver,               // fluid complex solver
                        _pressure_function,
                        _freq_domain_pressure_function,
                        _displ);
    
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
    
    _initialized = true;
}




MAST::PlateFSIFlutterSizingOptimization::~PlateFSIFlutterSizingOptimization() {

    if (!_initialized)
        return;

    
    delete _fluid_eq_sys;
    delete _fluid_mesh;
    
    delete _fluid_discipline;
    delete _fluid_sys_init;
    
    delete _frequency_domain_fluid_assembly;
    delete _frequency_domain_elem_ops;
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
    
    delete _pressure_function;
    delete _freq_domain_pressure_function;
    
    delete _displ;
    delete _normal_rot;
    delete _pressure;
    
    delete _structural_eq_sys;
    delete _structural_mesh;
    
    delete _modal_assembly;
    delete _modal_elem_ops;
    
    delete _structural_discipline;
    delete _structural_sys_init;
    
    delete _m_card;
    delete _p_card;
    
    delete _dirichlet_left;
    delete _dirichlet_right;
    delete _dirichlet_top;
    delete _dirichlet_bottom;
    
    delete _rho_f;
    delete _E_f;
    delete _nu_f;
    delete _kappa_f;
    delete _hoff_f;
    
    
    delete _rho;
    delete _E;
    delete _nu;
    delete _zero;
    delete _kappa;
    
    
    // delete the basis vectors
    if (_basis.size())
        for (unsigned int i=0; i<_basis.size(); i++)
            if (_basis[i]) delete _basis[i];
    
    // delete the h_y station functions
    {
        std::vector<MAST::ConstantFieldFunction*>::iterator
        it  = _th_station_functions.begin(),
        end = _th_station_functions.end();
        for (; it != end; it++)  delete *it;
    }
    
    
    // delete the h_y station parameters
    {
        std::vector<MAST::Parameter*>::iterator
        it  = _th_station_parameters.begin(),
        end = _th_station_parameters.end();
        for (; it != end; it++)  delete *it;
    }
    
    delete _weight;
    
    delete _flutter_solver;
    delete _gaf_database;
    delete _augment_send_list_obj;
}



void
MAST::PlateFSIFlutterSizingOptimization::init_dvar(std::vector<Real>& x,
                                                  std::vector<Real>& xmin,
                                                  std::vector<Real>& xmax) {
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}



void
MAST::PlateFSIFlutterSizingOptimization::evaluate(const std::vector<Real>& dvars,
                                                 Real& obj,
                                                 bool eval_obj_grad,
                                                 std::vector<Real>& obj_grad,
                                                 std::vector<Real>& fvals,
                                                 std::vector<bool>& eval_grads,
                                                 std::vector<Real>& grads) {
    
    libmesh_assert(_initialized);
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;
    
    // the optimization problem is defined as
    // min weight, subject to constraints on displacement and stresses
    Real
    wt      = 0.,
    pval    = 2.;
    
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        (*_th_station_parameters[i]) = dvars[i]*_dv_scaling[i];
    
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "th     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << (*_th_station_parameters[i])() << std::endl;
    
    // calculate weight
    (*_weight)(pt, 0., wt);



    ///////////////////////////////////////////////////////////////////
    // FLUTTER SOLUTION
    ///////////////////////////////////////////////////////////////////
    // clear flutter solver and set the output file
    _flutter_solver->clear();
    
    std::ostringstream oss;
    oss << "flutter_output_" << __init->comm().rank() << ".txt";
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(oss.str());
    
    _gaf_database->attach_discipline_and_system(*_gaf_database,
                                                *_structural_discipline,
                                                *_structural_sys_init);
    
    libMesh::NumericVector<Real>&
    base_sol = _fluid_sys->get_vector("fluid_base_solution");
    _frequency_domain_fluid_assembly->attach_discipline_and_system(*_frequency_domain_elem_ops,
                                                                   *_fluid_discipline,
                                                                   *_complex_solver,
                                                                   *_fluid_sys_init);
    _frequency_domain_fluid_assembly->set_base_solution(base_sol);
    _frequency_domain_elem_ops->set_frequency_function(*_freq_function);
    _pressure_function->init(base_sol);
    

    _gaf_database->init(_freq_function,
                        _complex_solver,                       // fluid complex solver
                        _pressure_function,
                        _freq_domain_pressure_function,
                        _displ);
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
    
    std::pair<bool, MAST::FlutterRootBase*>
    sol = _flutter_solver->find_critical_root(tol, max_bisection_iters);
    
    
    _flutter_solver->print_sorted_roots();
    _gaf_database->clear_discipline_and_system();
    _flutter_solver->clear_assembly_object();
    _frequency_domain_fluid_assembly->clear_discipline_and_system();

    // copy the flutter velocity to the contraint vector
    //     Vf        >= V0
    // or, V0/Vf     <= 1
    // or, V0/Vf - 1 <= 0
    //
    if (sol.second)
        fvals[0]  =  _V0_flutter/sol.second->V - 1.;
    else
        fvals[0]  =  -100.;

    // tell all ranks about the constraint function values
    this->comm().broadcast(fvals);

    
    //////////////////////////////////////////////////////////////////////
    // get the objective and constraints
    //////////////////////////////////////////////////////////////////////
    
    // set the function and objective values
    obj = wt;

    // tell all ranks about the obj function
    this->comm().sum(obj);
    
    //////////////////////////////////////////////////////////////////
    //   evaluate sensitivity if needed
    //////////////////////////////////////////////////////////////////
    
    // sensitivity of the objective function
    if (eval_obj_grad) {
        
        Real w_sens = 0.;
        std::fill(obj_grad.begin(), obj_grad.end(), 0.);
        
        // set gradient of weight
        for (unsigned int i=0; i<_n_vars; i++) {
            _weight->derivative(*_th_station_parameters[i],
                                pt,
                                0.,
                                w_sens);
            obj_grad[i] = w_sens*_dv_scaling[i];
        }

        // tell all processors about the sens values
        this->comm().sum(obj_grad);
    }


    
    // now check if the sensitivity of constraint function is requested
    bool if_sens = false;
    
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_sens = (if_sens || eval_grads[i]);
    
    if (if_sens) {
        
        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        
        // we are going to choose to use one parametric sensitivity at a time
        for (unsigned int i=0; i<_n_vars; i++) {
            
            libMesh::ParameterVector params;
            
            params.resize(1);
            params[0]  = _th_station_parameters[i]->ptr();
            
            // calculate sensitivity only if the flutter root was found
            // else, set it to zero
            if (sol.second) {
                
                _gaf_database->attach_discipline_and_system(*_gaf_database,
                                                            *_structural_discipline,
                                                            *_structural_sys_init);
                
                _frequency_domain_fluid_assembly->attach_discipline_and_system
                (*_frequency_domain_elem_ops,
                 *_fluid_discipline,
                 *_complex_solver,
                 *_fluid_sys_init);
                _frequency_domain_fluid_assembly->set_base_solution(base_sol);
                _frequency_domain_elem_ops->set_frequency_function(*_freq_function);
                
                _gaf_database->init(_freq_function,
                                    _complex_solver,                       // fluid complex solver
                                    _pressure_function,
                                    _freq_domain_pressure_function,
                                    _displ);
                _flutter_solver->attach_assembly(*_gaf_database);
                _flutter_solver->initialize(*_omega,
                                            *_b_ref,
                                            _flight_cond->rho(),
                                            _k_lower,         // lower kr
                                            _k_upper,         // upper kr
                                            _n_k_divs,        // number of divisions
                                            _basis);          // basis vectors

                _flutter_solver->calculate_sensitivity(*sol.second, params, 0);

                _gaf_database->clear_discipline_and_system();
                _flutter_solver->clear_assembly_object();
                _frequency_domain_fluid_assembly->clear_discipline_and_system();

                // copy the sensitivity values in the output
                grads[i] = -_dv_scaling[i] *
                _V0_flutter / pow(sol.second->V,2) * sol.second->V_sens;
            }
            else
                grads[i] = 0.;
        }
        
        // tell all ranks about the gradients
        this->comm().broadcast(grads);
    }
}








void
MAST::PlateFSIFlutterSizingOptimization::output(unsigned int iter,
                                               const std::vector<Real>& x,
                                               Real obj,
                                               const std::vector<Real>& fval,
                                               bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
}


MAST::FunctionEvaluation::funobj
MAST::PlateFSIFlutterSizingOptimization::get_objective_evaluation_function() {
    
    return plate_euler_flutter_optim_obj;
}



MAST::FunctionEvaluation::funcon
MAST::PlateFSIFlutterSizingOptimization::get_constraint_evaluation_function() {
    
    return plate_euler_flutter_optim_con;
}

