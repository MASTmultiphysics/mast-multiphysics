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
#include <ostream>


// MAST includes
#include "examples/fsi/plate_flutter_solution/plate_euler_fsi_flutter_solution.h"
#include "examples/fluid/meshing/panel_mesh_3D.h"
#include "base/nonlinear_system.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "fluid/pressure_function.h"
#include "fluid/frequency_domain_pressure_function.h"
#include "base/complex_assembly_base.h"
#include "solver/complex_solver_base.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "base/complex_mesh_field_function.h"
#include "elasticity/complex_normal_rotation_mesh_function.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/ug_flutter_root.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_element_base.h"
#include "base/eigenproblem_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "aeroelasticity/ug_flutter_solver.h"
#include "examples/base/augment_ghost_elem_send_list.h"
#include "solver/slepc_eigen_solver.h"
#include "examples/base/plot_results.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_solver.h"


extern libMesh::LibMeshInit* __init;



MAST::PlateEulerFSIFlutterAnalysis::PlateEulerFSIFlutterAnalysis():
_structural_mesh                       (nullptr),
_fluid_mesh                            (nullptr),
_structural_eq_sys                     (nullptr),
_fluid_eq_sys                          (nullptr),
_structural_sys                        (nullptr),
_fluid_sys                             (nullptr),
_structural_sys_init                   (nullptr),
_structural_discipline                 (nullptr),
_fluid_sys_init                        (nullptr),
_fluid_discipline                      (nullptr),
_flight_cond                           (nullptr),
_far_field                             (nullptr),
_symm_wall                             (nullptr),
_slip_wall                             (nullptr),
_pressure                              (nullptr),
_displ                                 (nullptr),
_normal_rot                            (nullptr),
_pressure_function                     (nullptr),
_freq_domain_pressure_function         (nullptr),
_omega                                 (nullptr),
_velocity                              (nullptr),
_b_ref                                 (nullptr),
_omega_f                               (nullptr),
_velocity_f                            (nullptr),
_b_ref_f                               (nullptr),
_freq_function                         (nullptr),
_length                                (0.),
_k_lower                               (0.),
_k_upper                               (0.),
_n_k_divs                              (0),
_th                                    (nullptr),
_rho                                   (nullptr),
_E                                     (nullptr),
_nu                                    (nullptr),
_kappa                                 (nullptr),
_zero                                  (nullptr),
_mach                                  (nullptr),
_rho_air                               (nullptr),
_gamma_air                             (nullptr),
_th_f                                  (nullptr),
_rho_f                                 (nullptr),
_E_f                                   (nullptr),
_nu_f                                  (nullptr),
_kappa_f                               (nullptr),
_hoff_f                                (nullptr),
_mach_f                                (nullptr),
_rho_air_f                             (nullptr),
_gamma_air_f                           (nullptr),
_flutter_solver                        (nullptr),
_flutter_root                          (nullptr),
_m_card                                (nullptr),
_p_card                                (nullptr),
_dirichlet_left                        (nullptr),
_dirichlet_right                       (nullptr),
_dirichlet_bottom                      (nullptr),
_dirichlet_top                         (nullptr),
_augment_send_list_obj                 (nullptr) {

    
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
    
    std::unique_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MeshInitializer::CoordinateDivisions),
    z_coord_divs    (new MeshInitializer::CoordinateDivisions);
    
    std::vector<MeshInitializer::CoordinateDivisions*>
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
    
    _augment_send_list_obj  = new MAST::AugmentGhostElementSendListObj(*_fluid_sys);
    _fluid_sys->get_dof_map().attach_extra_send_list_object(*_augment_send_list_obj);
    
    // initialize the equation system for analysis
    _fluid_eq_sys->init();
    
    // print the information
    _fluid_eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD);
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
    _structural_discipline = new MAST::PhysicsDisciplineBase(*_structural_eq_sys);
    
    
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
    
    // create the property functions and add them to the
    
    _th              = new MAST::Parameter("th", 0.0015);
    _rho             = new MAST::Parameter("rho", 2.7e3);
    _E               = new MAST::Parameter("E",   72.e9);
    _nu              = new MAST::Parameter("nu",   0.33);
    _kappa           = new MAST::Parameter("kappa",  5./6.);
    _zero            = new MAST::Parameter("zero",   0.);
    _mach            = new MAST::Parameter("mach",   3.);
    _rho_air         = new MAST::Parameter("rho" , 1.05);
    _gamma_air       = new MAST::Parameter("gamma", 1.4);
    
    
    
    // prepare the vector of parameters with respect to which the sensitivity
    // needs to be benchmarked
    _params_for_sensitivity.push_back(_E);
    _params_for_sensitivity.push_back(_nu);
    _params_for_sensitivity.push_back(_th);
    
    
    
    _th_f            = new MAST::ConstantFieldFunction("h",            *_th);
    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",     *_kappa);
    _hoff_f          = new MAST::ConstantFieldFunction("off",        *_zero);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);
    
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
    
    _flutter_solver  = new MAST::UGFlutterSolver;
    std::ostringstream oss;
    oss << "flutter_output_" << __init->comm().rank() << ".txt";
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(oss.str());
}






MAST::PlateEulerFSIFlutterAnalysis::~PlateEulerFSIFlutterAnalysis() {
    
    delete _fluid_eq_sys;
    delete _fluid_mesh;
    
    delete _fluid_discipline;
    delete _fluid_sys_init;
    
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
    
    delete _structural_discipline;
    delete _structural_sys_init;
    
    delete _m_card;
    delete _p_card;
    
    delete _dirichlet_left;
    delete _dirichlet_right;
    delete _dirichlet_top;
    delete _dirichlet_bottom;
    
    delete _th_f;
    delete _rho_f;
    delete _E_f;
    delete _nu_f;
    delete _kappa_f;
    delete _hoff_f;
    delete _mach_f;
    delete _rho_air_f;
    delete _gamma_air_f;
    
    
    delete _th;
    delete _rho;
    delete _E;
    delete _nu;
    delete _kappa;
    delete _zero;
    delete _mach;
    delete _rho_air;
    delete _gamma_air;
    
    
    // delete the basis vectors
    if (_basis.size())
        for (unsigned int i=0; i<_basis.size(); i++)
            delete _basis[i];
    
    delete _flutter_solver;
    delete _augment_send_list_obj;
}



Real
MAST::PlateEulerFSIFlutterAnalysis::solve(bool if_write_output,
                                          const Real tol,
                                          const unsigned int max_bisection_iters) {
    
    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID SOLUTION
    /////////////////////////////////////////////////////////////////
    
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
    MAST::ComplexAssemblyBase                                     assembly;
    MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations  elem_ops;
    
    // Transient solver for time integration
    MAST::ComplexSolverBase                          solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(elem_ops,
                                          *_fluid_discipline,
                                          solver,
                                          *_fluid_sys_init);
    assembly.set_base_solution(base_sol);
    elem_ops.set_frequency_function(*_freq_function);
    _pressure_function->init(base_sol);


    
    
    ////////////////////////////////////////////////////////////
    // STRUCTURAL MODAL EIGENSOLUTION
    ////////////////////////////////////////////////////////////
    
    // create the nonlinear assembly object
    MAST::EigenproblemAssembly   modal_assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations modal_elem_ops;
    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
    
    modal_assembly.attach_discipline_and_system(modal_elem_ops,
                                                *_structural_discipline,
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
    MAST::FSIGeneralizedAeroForceAssembly fsi_assembly;
    _flutter_solver->clear();
    
    fsi_assembly.attach_discipline_and_system(fsi_assembly,
                                              *_structural_discipline,
                                              *_structural_sys_init);
    std::ostringstream oss;
    oss << "flutter_output_" << __init->comm().rank() << ".txt";
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(oss.str());
    
    fsi_assembly.init(&solver,                       // fluid complex solver
                      _pressure_function,
                      _freq_domain_pressure_function,
                      _displ);
    _flutter_solver->attach_assembly(fsi_assembly);
    _flutter_solver->initialize(*_omega,
                                *_b_ref,
                                _flight_cond->rho(),
                                _k_lower,            // lower kr
                                _k_upper,            // upper kr
                                _n_k_divs,           // number of divisions
                                _basis);             // basis vectors
    
    
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

    
    if (sol.first && if_write_output) {
        
        MAST::plot_structural_flutter_solution("structural_flutter_mode.exo",
                                               *_structural_sys,
                                               sol.second->eig_vec_right,
                                               _basis);
        MAST::plot_fluid_flutter_solution("fluid_flutter_mode.exo",
                                          *_structural_sys,
                                          *_fluid_sys,
                                          *_displ,
                                          solver,
                                          sol.second->eig_vec_right,
                                          _basis);
    }
    

    
    
    assembly.clear_discipline_and_system();
    
    
    return _flutter_root->V;
}





Real
MAST::PlateEulerFSIFlutterAnalysis::sensitivity_solve(MAST::Parameter& p) {
    
    /*
     // create the nonlinear assembly object
     MAST::StructuralNonlinearAssembly   assembly;
     
     assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
     
     MAST::NonlinearSystem&  nonlin_sys   =
     assembly.system();
     
     libMesh::ParameterVector params;
     params.resize(1);
     params[0]  =  p.ptr();
     
     // zero the solution before solving
     nonlin_sys.add_sensitivity_solution(0).zero();
     this->clear_stresss();
     
     nonlin_sys.sensitivity_solve(params);
     
     // evaluate sensitivity of the outputs
     assembly.calculate_output_sensitivity(params,
     true,    // true for total sensitivity
     *(_sys->solution));
     
     
     assembly.clear_discipline_and_system();
     
     // write the solution for visualization
     if (if_write_output) {
     
     std::ostringstream oss1, oss2;
     oss1 << "output_" << p.name() << ".exo";
     oss2 << "output_" << p.name() << ".exo";
     
     libMesh::out
     << "Writing sensitivity output to : " << oss1.str()
     << "  and stress/strain sensitivity to : " << oss2.str()
     << std::endl;
     
     
     _sys->solution->swap(_sys->get_sensitivity_solution(0));
     
     // write the solution for visualization
     _discipline->update_stress_strain_data( &p);
     libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
     *_eq_sys);
     
     
     _sys->solution->swap(_sys->get_sensitivity_solution(0));
     }
     */
    return 0.;
}

