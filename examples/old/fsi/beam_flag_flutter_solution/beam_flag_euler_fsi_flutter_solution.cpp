///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2019  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//// C++ includes
//#include <iostream>
//#include <ostream>
//
//
//// MAST includes
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_euler_fsi_flutter_solution.h"
//#include "examples/fluid/meshing/flag_mesh_2D.h"
//#include "examples/fsi/beam_flutter_solution/constrain_beam_dofs.h"
//#include "base/nonlinear_system.h"
//#include "fluid/conservative_fluid_system_initialization.h"
//#include "fluid/conservative_fluid_discipline.h"
//#include "base/complex_assembly_base.h"
//#include "fluid/frequency_domain_linearized_complex_assembly.h"
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_pressure_function.h"
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_frequency_domain_pressure_function.h"
//#include "solver/complex_solver_base.h"
//#include "fluid/flight_condition.h"
//#include "base/parameter.h"
//#include "base/constant_field_function.h"
//#include "base/boundary_condition_base.h"
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_frequency_domain_displacement.h"
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_frequency_domain_normal_rotation.h"
//#include "aeroelasticity/frequency_function.h"
//#include "aeroelasticity/ug_flutter_root.h"
//#include "elasticity/structural_system_initialization.h"
//#include "elasticity/structural_element_base.h"
//#include "base/eigenproblem_assembly.h"
//#include "elasticity/structural_modal_eigenproblem_assembly.h"
//#include "elasticity/fsi_generalized_aero_force_assembly.h"
//#include "elasticity/structural_near_null_vector_space.h"
//#include "property_cards/solid_1d_section_element_property_card.h"
//#include "property_cards/isotropic_material_property_card.h"
//#include "boundary_condition/dirichlet_boundary_condition.h"
//#include "aeroelasticity/ug_flutter_solver.h"
//#include "examples/base/augment_ghost_elem_send_list.h"
//#include "solver/slepc_eigen_solver.h"
//#include "examples/base/plot_results.h"
//
//
//// libMesh includes
//#include "libmesh/mesh_generation.h"
//#include "libmesh/exodusII_io.h"
//#include "libmesh/numeric_vector.h"
//#include "libmesh/getpot.h"
//#include "libmesh/string_to_enum.h"
//#include "libmesh/nonlinear_solver.h"
//#include "libmesh/dof_map.h"
//
//
//
//
//MAST::BeamFlagEulerFSIFlutterAnalysis::
//BeamFlagEulerFSIFlutterAnalysis(const libMesh::Parallel::Communicator& comm_in,
//                                unsigned int order_increment,
//                                unsigned int n_refine):
//libMesh::ParallelObject                (comm_in),
//_structural_mesh                       (nullptr),
//_fluid_mesh                            (nullptr),
//_structural_eq_sys                     (nullptr),
//_fluid_eq_sys                          (nullptr),
//_structural_sys                        (nullptr),
//_fluid_sys                             (nullptr),
//_structural_sys_init                   (nullptr),
//_structural_discipline                 (nullptr),
//_fluid_sys_init                        (nullptr),
//_fluid_discipline                      (nullptr),
//_flight_cond                           (nullptr),
//_far_field                             (nullptr),
//_slip_wall                             (nullptr),
//_pressure                              (nullptr),
//_displ                                 (nullptr),
//_normal_rot                            (nullptr),
//_pressure_function                     (nullptr),
//_freq_domain_pressure_function         (nullptr),
//_omega                                 (nullptr),
//_velocity                              (nullptr),
//_b_ref                                 (nullptr),
//_omega_f                               (nullptr),
//_velocity_f                            (nullptr),
//_b_ref_f                               (nullptr),
//_freq_function                         (nullptr),
//_length                                (0.),
//_k_lower                               (0.),
//_k_upper                               (0.),
//_n_k_divs                              (0),
//_thy                                   (nullptr),
//_thz                                   (nullptr),
//_rho                                   (nullptr),
//_E                                     (nullptr),
//_nu                                    (nullptr),
//_kappa                                 (nullptr),
//_zero                                  (nullptr),
//_thy_f                                 (nullptr),
//_thz_f                                 (nullptr),
//_rho_f                                 (nullptr),
//_E_f                                   (nullptr),
//_nu_f                                  (nullptr),
//_kappa_f                               (nullptr),
//_hyoff_f                               (nullptr),
//_hzoff_f                               (nullptr),
//_flutter_solver                        (nullptr),
//_flutter_root                          (nullptr),
//_m_card                                (nullptr),
//_p_card                                (nullptr),
//_dirichlet_left                        (nullptr),
//_augment_send_list_obj                 (nullptr),
//_constraint_beam_dofs                  (nullptr) {
//    
//    //////////////////////////////////////////////////////////////////////
//    //    SETUP THE FLUID DATA
//    //////////////////////////////////////////////////////////////////////
//    
//    // initialize the libMesh object
//    _fluid_mesh              = new libMesh::ParallelMesh(this->comm());
//    _fluid_eq_sys            = new libMesh::EquationSystems(*_fluid_mesh);
//    
//    // add the system to be used for analysis
//    _fluid_sys = &(_fluid_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
//    _fluid_sys->set_init_B_matrix();
//    
//    
//    // initialize the flow conditions
//    GetPot infile("input.in");
//    
//    
//    const unsigned int
//    dim                 = 2,
//    nx_divs             = infile("nx_divs", 0),
//    ny_divs             = infile("ny_divs", 0),
//    panel_bc_id         = 10,
//    n_flags             = infile("n_flags", 0);
//    
//    const bool
//    bottom_slip_wall = infile("if_bottom_slip_wall", false),
//    top_slip_wall    = infile(   "if_top_slip_wall", false);
//    
//    
//    Real
//    tol                 = 1.e-6,
//    thickness           = 0.;
//    
//    std::string
//    etype;
//    
//    libMesh::ElemType
//    elem_type           =
//    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "QUAD4"));
//    
//    libMesh::FEFamily
//    fe_family           =
//    libMesh::Utility::string_to_enum<libMesh::FEFamily>(infile("fe_family", "SZABAB"));
//    
//    libMesh::Order
//    fe_order            =
//    libMesh::Utility::string_to_enum<libMesh::Order>(infile("fe_order", "FIRST"));
//    
//    // change these types depending on the order
//    if (order_increment > 0) {
//        
//        elem_type = libMesh::QUAD9;
//        fe_family = libMesh::SZABAB;
//        fe_order  = static_cast<libMesh::Order>(fe_order + order_increment);
//    }
//    
//    
//    std::vector<Real>
//    x_div_loc        (nx_divs+1),
//    x_relative_dx    (nx_divs+1),
//    y_div_loc        (ny_divs+1),
//    y_relative_dx    (ny_divs+1);
//    
//    std::vector<unsigned int>
//    x_divs           (nx_divs),
//    y_divs           (ny_divs);
//    
//    std::unique_ptr<MAST::MeshInitializer::CoordinateDivisions>
//    x_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions),
//    y_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions);
//    
//    std::vector<MAST::MeshInitializer::CoordinateDivisions*>
//    divs(dim);
//    
//    
//    // now read in the values: x-coord
//    for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
//        
//        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
//        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
//        
//        if (i_div < nx_divs) //  this is only till nx_divs
//            x_divs[i_div]       = pow(2, n_refine)*infile( "x_div_nelem", 0, i_div);
//    }
//    
//    divs[0] = x_coord_divs.get();
//    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
//    
//    
//    // now read in the values: y-coord
//    for (unsigned int i_div=0; i_div<ny_divs+1; i_div++) {
//        
//        y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
//        y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
//        
//        if (i_div < ny_divs) //  this is only till ny_divs
//            y_divs[i_div]    = pow(2, n_refine)*infile( "y_div_nelem",  0, i_div);
//    }
//    
//    divs[1] = y_coord_divs.get();
//    y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
//    
//    
//    
//    
//    // initialize the mesh
//    MAST::FlagMesh2D
//    flag_mesh  = MAST::FlagMesh2D(n_flags,
//                                  panel_bc_id,
//                                  divs);
//    flag_mesh.init_fluid_mesh(*_fluid_mesh, elem_type);
//    
//    // get the thickness and make sure that all flags are the same thickness
//    thickness           = flag_mesh.thickness(0);
//    
//    for (unsigned int i=0; i<n_flags; i++) {
//     
//        libmesh_assert_less_equal(std::abs(thickness-flag_mesh.thickness(i)),
//                                  tol);
//    }
//    
//    _fluid_discipline   = new MAST::ConservativeFluidDiscipline(*_fluid_eq_sys);
//    _fluid_sys_init     = new MAST::ConservativeFluidSystemInitialization(*_fluid_sys,
//                                                                          _fluid_sys->name(),
//                                                                          libMesh::FEType(fe_order, fe_family),
//                                                                          dim);
//    _augment_send_list_obj = new MAST::AugmentGhostElementSendListObj(*_fluid_sys);
//    
//    _fluid_sys->get_dof_map().attach_extra_send_list_object(*_augment_send_list_obj);
//    
//    // initialize the equation system for analysis
//    _fluid_eq_sys->init();
//    
//    // print the information
//    _fluid_mesh->print_info();
//    _fluid_eq_sys->print_info();
//    
//    // create the oundary conditions for slip-wall and far-field
//    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD);
//    _slip_wall     = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
//    _symm_wall     = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
//    
//    _flight_cond    =  new MAST::FlightCondition;
//    for (unsigned int i=0; i<3; i++) {
//        
//        _flight_cond->body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
//        _flight_cond->body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
//        _flight_cond->body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
//        _flight_cond->body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
//        _flight_cond->body_angular_rates(i) = infile("body_angular_rates", 0., i);
//    }
//    
//    _flight_cond->ref_chord       = infile("ref_c",    1.);
//    _flight_cond->mach            = infile("mach",     .5);
//    _flight_cond->gas_property.cp = infile(  "cp",  1003.);
//    _flight_cond->gas_property.cv = infile(  "cv",   716.);
//    _flight_cond->gas_property.T  = infile("temp",   300.);
//    _flight_cond->gas_property.rho= infile( "rho",   1.05);
//    
//    _flight_cond->init();
//    
//    // tell the discipline about the fluid values
//    _fluid_discipline->set_flight_condition(*_flight_cond);
//    
//    // define parameters
//    _omega             = new MAST::Parameter("omega",       0.);
//    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
//    _b_ref             = new MAST::Parameter("b_ref",       1.);
//    
//    
//    // now define the constant field functions based on this
//    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
//    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
//    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
//    
//    // initialize the frequency function
//    _freq_function     = new MAST::FrequencyFunction("freq",
//                                                     *_omega_f,
//                                                     *_velocity_f,
//                                                     *_b_ref_f);
//    _freq_function->if_nondimensional(true);
//    
//    // tell the physics about boundary conditions
//    _fluid_discipline->add_side_load(    panel_bc_id, *_slip_wall);
//    // left and right boundaries are far-field
//    _fluid_discipline->add_side_load(              1, *_far_field); // right
//    _fluid_discipline->add_side_load(              3, *_far_field); // left
//    
//    // upper and lower boundary types are identified from the input file
//    if (bottom_slip_wall)
//        _fluid_discipline->add_side_load(              0, *_symm_wall); // bottom
//    else
//        _fluid_discipline->add_side_load(              0, *_far_field); // bottom
//    
//    if (top_slip_wall)
//        _fluid_discipline->add_side_load(              2, *_symm_wall); // top
//    else
//        _fluid_discipline->add_side_load(              2, *_far_field); // top
//    
//    _pressure_function =
//    new MAST::BeamFlagPressureFunction(*_fluid_sys_init,
//                                       *_flight_cond,
//                                       thickness);
//    _freq_domain_pressure_function =
//    new MAST::BeamFlagFrequencyDomainPressureFunction(*_fluid_sys_init,
//                                                      *_flight_cond,
//                                                      thickness);
//    
//    _pressure_function->set_calculate_cp(true);
//    _freq_domain_pressure_function->set_calculate_cp(true);
//    
//    _k_upper            = infile("k_upper",  0.75);
//    _k_lower            = infile("k_lower",  0.05);
//    _n_k_divs           = infile("n_k_divs",   10);
//    
//    
//    //////////////////////////////////////////////////////////////////////
//    //    SETUP THE STRUCTURAL DATA
//    //////////////////////////////////////////////////////////////////////
//    
//    // setup length for use in setup of flutter solver
//    _length = x_div_loc[2]-x_div_loc[1];
//    (*_b_ref) = _length;
//    
//    // create the mesh
//    _structural_mesh       = new libMesh::SerialMesh(this->comm());
//    
//    
//    // change these types depending on the order
//    if (order_increment > 0)
//        elem_type = libMesh::EDGE3;
//    else
//        elem_type = libMesh::EDGE2;
//    
//    flag_mesh.init_structural_mesh(*_structural_mesh, elem_type);
//    
//    // create the equation system
//    _structural_eq_sys    = new  libMesh::EquationSystems(*_structural_mesh);
//    
//    // create the libmesh system
//    _structural_sys       = &(_structural_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
//    _structural_sys->set_eigenproblem_type(libMesh::GHEP);
//    
//    
//    
//    // initialize the system to the right set of variables
//    _structural_sys_init  = new MAST::StructuralSystemInitialization(*_structural_sys,
//                                                                     _structural_sys->name(),
//                                                                     libMesh::FEType(fe_order, fe_family));
//    _structural_discipline = new MAST::PhysicsDisciplineBase(*_structural_eq_sys);
//    
//    
//    // create and add the boundary condition and loads
//    _dirichlet_left = new MAST::DirichletBoundaryCondition;
//    _dirichlet_left->init (0, _structural_sys_init->vars());
//    _structural_discipline->add_dirichlet_bc(0, *_dirichlet_left);
//    _structural_discipline->init_system_dirichlet_bc(*_structural_sys);
//    _constraint_beam_dofs = new MAST::ConstrainBeamDofs(*_structural_sys);
//    _structural_sys->attach_constraint_object(*_constraint_beam_dofs);
//
//    // initialize the equation system
//    _structural_eq_sys->init();
//    _structural_mesh->print_info();
//    _structural_eq_sys->print_info();
//    
//    // initialize the motion object
//    _displ        =
//    new MAST::BeamFlagFrequencyDomainDisplacement(*_structural_sys_init,
//                                                  "frequency_domain_displacement",
//                                                  flag_mesh.midplane_coordinates());
//    _normal_rot   =
//    new MAST::BeamFlagFrequencyDomainNormalRotation("frequency_domain_normal_rotation",
//                                                    *_displ,
//                                                    flag_mesh.midplane_coordinates());
//    _slip_wall->add(*_displ);
//    _slip_wall->add(*_normal_rot);
//    
//    
//    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
//    _structural_sys->set_exchange_A_and_B(true);
//    _structural_sys->set_n_requested_eigenvalues(infile("n_modes", 20));
//    
//    // create the property functions and add them to the
//    
//    _thy             = new MAST::Parameter("thy",              thickness);
//    _thz             = new MAST::Parameter("thz",                   1.00);
//    _rho             = new MAST::Parameter("rho", infile("rho_s", 2.7e3));
//    _E               = new MAST::Parameter("E",   infile("E_s",   72.e9));
//    _nu              = new MAST::Parameter("nu",  infile("nu_s",   0.33));
//    _kappa           = new MAST::Parameter("kappa",   5./6.);
//    _zero            = new MAST::Parameter("zero",       0.);
//    
//    
//    
//    // prepare the vector of parameters with respect to which the sensitivity
//    // needs to be benchmarked
//    _params_for_sensitivity.push_back(_E);
//    _params_for_sensitivity.push_back(_nu);
//    _params_for_sensitivity.push_back(_thy);
//    _params_for_sensitivity.push_back(_thz);
//    
//    
//    
//    _thy_f           = new MAST::ConstantFieldFunction("hy",          *_thy);
//    _thz_f           = new MAST::ConstantFieldFunction("hz",          *_thz);
//    _rho_f           = new MAST::ConstantFieldFunction("rho",         *_rho);
//    _E_f             = new MAST::ConstantFieldFunction("E",             *_E);
//    _nu_f            = new MAST::ConstantFieldFunction("nu",           *_nu);
//    _kappa_f         = new MAST::ConstantFieldFunction("kappa",     *_kappa);
//    _hyoff_f         = new MAST::ConstantFieldFunction("hy_off",     *_zero);
//    _hzoff_f         = new MAST::ConstantFieldFunction("hz_off",     *_zero);
//    
//    // create the material property card
//    _m_card          = new MAST::IsotropicMaterialPropertyCard;
//    
//    // add the material properties to the card
//    _m_card->add(*_rho_f);
//    _m_card->add(*_E_f);
//    _m_card->add(*_nu_f);
//    _m_card->add(*_kappa_f);
//    
//    // create the element property card
//    _p_card          = new MAST::Solid1DSectionElementPropertyCard;
//    //_p_card->set_bending_model(MAST::TIMOSHENKO);
//    
//    // tell the card about the orientation
//    libMesh::Point orientation;
//    orientation(1) = 1.;
//    _p_card->y_vector() = orientation;
//    
//    // add the section properties to the card
//    _p_card->add(*_thy_f);
//    _p_card->add(*_thz_f);
//    _p_card->add(*_hyoff_f);
//    _p_card->add(*_hzoff_f);
//    
//    // tell the section property about the material property
//    _p_card->set_material(*_m_card);
//    
//    _p_card->init();
//    
//    _structural_discipline->set_property_for_subdomain(0, *_p_card);
//    
//    // pressure boundary condition for the beam
//    _pressure    =  new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
//    _pressure->add(*_pressure_function);
//    _pressure->add(*_freq_domain_pressure_function);
//    _structural_discipline->add_volume_load(0, *_pressure);
//    
//    _flutter_solver  = new MAST::UGFlutterSolver;
//}
//
//
//
//
//
//
//MAST::BeamFlagEulerFSIFlutterAnalysis::~BeamFlagEulerFSIFlutterAnalysis() {
//    
//    delete _fluid_eq_sys;
//    delete _fluid_mesh;
//    
//    delete _fluid_discipline;
//    delete _fluid_sys_init;
//    
//    delete _far_field;
//    delete _slip_wall;
//    delete _symm_wall;
//    
//    delete _flight_cond;
//    
//    delete _omega;
//    delete _velocity;
//    delete _b_ref;
//    
//    delete _omega_f;
//    delete _velocity_f;
//    delete _b_ref_f;
//    
//    delete _freq_function;
//    
//    delete _pressure_function;
//    delete _freq_domain_pressure_function;
//    
//    
//    delete _displ;
//    delete _normal_rot;
//    delete _pressure;
//    
//    delete _structural_eq_sys;
//    delete _structural_mesh;
//    delete _constraint_beam_dofs;
//
//    
//    delete _structural_discipline;
//    delete _structural_sys_init;
//    
//    delete _m_card;
//    delete _p_card;
//    
//    delete _dirichlet_left;
//    
//    delete _thy_f;
//    delete _thz_f;
//    delete _rho_f;
//    delete _E_f;
//    delete _nu_f;
//    delete _kappa_f;
//    delete _hyoff_f;
//    delete _hzoff_f;
//    
//    
//    delete _thy;
//    delete _thz;
//    delete _rho;
//    delete _E;
//    delete _nu;
//    delete _kappa;
//    delete _zero;
//    
//    
//    // delete the basis vectors
//    if (_basis.size())
//        for (unsigned int i=0; i<_basis.size(); i++)
//            if (_basis[i]) delete _basis[i];
//    
//    delete _flutter_solver;
//    
//    delete _augment_send_list_obj;
//}
//
//
//
//
//Real
//MAST::BeamFlagEulerFSIFlutterAnalysis::solve(bool if_write_output,
//                                             const Real tol,
//                                             const unsigned int max_bisection_iters) {
//    
//    
//    /////////////////////////////////////////////////////////////////
//    //  INITIALIZE FLUID SOLUTION
//    /////////////////////////////////////////////////////////////////
//    // the modal and flutter problems are solved on rank 0, while
//    // the fluid solution is setup on the global communicator
//    
//    // initialize the solution
//    RealVectorX s = RealVectorX::Zero(4);
//    s(0) = _flight_cond->rho();
//    s(1) = _flight_cond->rho_u1();
//    s(2) = _flight_cond->rho_u2();
//    s(3) = _flight_cond->rho_e();
//    
//    // create the vector for storing the base solution.
//    // we will swap this out with the system solution, initialize and
//    // then swap it back.
//    libMesh::NumericVector<Real>& base_sol =
//    _fluid_sys->add_vector("fluid_base_solution");
//    _fluid_sys->solution->swap(base_sol);
//    _fluid_sys_init->initialize_solution(s);
//    _fluid_sys->solution->swap(base_sol);
//    
//    // create the nonlinear assembly object
//    MAST::ComplexAssemblyBase                                      assembly;
//    MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations   elem_ops;
//    
//    // solver for complex solution
//    MAST::ComplexSolverBase                          solver;
//    
//    // now setup the assembly object
//    assembly.set_discipline_and_system(elem_ops,
//                                          *_fluid_discipline,
//                                          solver,
//                                          *_fluid_sys_init);
//    assembly.set_base_solution(base_sol);
//    elem_ops.set_frequency_function(*_freq_function);
//    _pressure_function->init(base_sol);
//    
//    
//    
//    
//    ////////////////////////////////////////////////////////////
//    // STRUCTURAL MODAL EIGENSOLUTION
//    ////////////////////////////////////////////////////////////
//    
//    // create the nonlinear assembly object
//    MAST::EigenproblemAssembly                                modal_assembly;
//    MAST::StructuralModalEigenproblemAssemblyElemOperations   modal_elem_ops;
//    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
//    
//    modal_assembly.set_discipline_and_system(modal_elem_ops,
//                                                *_structural_discipline,
//                                                *_structural_sys_init);
//    
//    
//    MAST::StructuralNearNullVectorSpace nsp;
//    _structural_sys->nonlinear_solver->nearnullspace_object = &nsp;
//    
//    _structural_sys->eigenproblem_solve();
//    modal_assembly.clear_discipline_and_system();
//    
//    // Get the number of converged eigen pairs.
//    unsigned int
//    nconv = std::min(_structural_sys->get_n_converged_eigenvalues(),
//                     _structural_sys->get_n_requested_eigenvalues());
//    
//    if (_basis.size() > 0)
//        libmesh_assert(_basis.size() == nconv);
//    else {
//        _basis.resize(nconv);
//        for (unsigned int i=0; i<_basis.size(); i++)
//            _basis[i] = nullptr;
//    }
//    
//    libMesh::ExodusII_IO*
//    writer = nullptr;
//    
//    if (if_write_output)
//        writer = new libMesh::ExodusII_IO(*_structural_mesh);
//    
//    for (unsigned int i=0; i<nconv; i++) {
//        
//        // create a vector to store the basis
//        if (_basis[i] == nullptr)
//            _basis[i] = _structural_sys->solution->zero_clone().release();
//        
//        // now write the eigenvalue
//        Real
//        re = 0.,
//        im = 0.;
//        _structural_sys->get_eigenpair(i, re, im, *_basis[i]);
//        
//        libMesh::out
//        << std::setw(35) << std::fixed << std::setprecision(15)
//        << re << std::endl;
//        
//        if (if_write_output) {
//            
//            // We write the file in the ExodusII format.
//            // copy the solution for output
//            _structural_sys->solution->swap(*_basis[i]);
//            writer->write_timestep("modes.exo",
//                                   *_structural_eq_sys,
//                                   i+1, i);
//            _structural_sys->solution->swap(*_basis[i]);
//        }
//    }
//    
//    ///////////////////////////////////////////////////////////////////
//    // FLUTTER SOLUTION
//    ///////////////////////////////////////////////////////////////////
//    // clear flutter solver and set the output file
//    _flutter_solver->clear();
//    
//    MAST::FSIGeneralizedAeroForceAssembly fsi_assembly;
//    fsi_assembly.set_discipline_and_system(fsi_assembly,
//                                              *_structural_discipline,
//                                              *_structural_sys_init);
//    
//    std::ostringstream oss;
//    oss << "flutter_output_" << this->comm().rank() << ".txt";
//    if (this->comm().rank() == 0)
//        _flutter_solver->set_output_file(oss.str());
//    
//    
//    fsi_assembly.init(&solver,                       // fluid complex solver
//                      _pressure_function,
//                      _freq_domain_pressure_function,
//                      _displ);
//    _flutter_solver->attach_assembly(fsi_assembly);
//    _flutter_solver->initialize(*_omega,
//                                *_b_ref,
//                                _flight_cond->rho(),
//                                _k_lower,         // lower kr
//                                _k_upper,         // upper kr
//                                _n_k_divs,        // number of divisions
//                                _basis);          // basis vectors
//    
//    
//    // find the roots for the specified divisions
//    _flutter_solver->scan_for_roots();
//    _flutter_solver->print_crossover_points();
//    
//    // now ask the flutter solver to return the critical flutter root,
//    // which is the flutter cross-over point at the lowest velocity
//    std::pair<bool, MAST::FlutterRootBase*>
//    sol = _flutter_solver->find_critical_root(tol, max_bisection_iters);
//    
//    
//    _flutter_solver->print_sorted_roots();
//    fsi_assembly.clear_discipline_and_system();
//    _flutter_solver->clear_assembly_object();
//    
//    // make sure solution was found
//    libmesh_assert(sol.first);
//    _flutter_root = sol.second;
//    
//    
//    if (sol.first && if_write_output) {
//        
//        MAST::plot_structural_flutter_solution("structural_flutter_mode.exo",
//                                               *_structural_sys,
//                                               sol.second->eig_vec_right,
//                                               _basis);
//        MAST::plot_fluid_flutter_solution("fluid_flutter_mode.exo",
//                                          *_structural_sys,
//                                          *_fluid_sys,
//                                          *_displ,
//                                          solver,
//                                          sol.second->eig_vec_right,
//                                          _basis);
//    }
//    
//    
//    
//    
//    assembly.clear_discipline_and_system();
//    
//    if (sol.first)
//        return _flutter_root->V;
//    else
//        return 0.;
//}
//
//
//
//
//
//Real
//MAST::BeamFlagEulerFSIFlutterAnalysis::sensitivity_solve(MAST::Parameter& p) {
//    
//    //Make sure that  a solution is available for sensitivity
//    libmesh_assert(_flutter_root);
//    
//    /////////////////////////////////////////////////////////////////
//    //  INITIALIZE FLUID ASSEMBLY/SOLVER OBJECTS
//    /////////////////////////////////////////////////////////////////
//    libMesh::NumericVector<Real>& base_sol =
//    _fluid_sys->get_vector("fluid_base_solution");
//    
//    // create the nonlinear assembly object
//    MAST::ComplexAssemblyBase                                     assembly;
//    MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations  elem_ops;
//    
//    // Transient solver for time integration
//    MAST::ComplexSolverBase                          solver;
//    
//    // now solve the system
//    assembly.set_discipline_and_system(elem_ops,
//                                          *_fluid_discipline,
//                                          solver,
//                                          *_fluid_sys_init);
//    assembly.set_base_solution(base_sol);
//    elem_ops.set_frequency_function(*_freq_function);
//    
//    
//    // it is assumed that the modal basis is available from the flutter solution
//    
//    
//    ///////////////////////////////////////////////////////////////////
//    // FLUTTER SOLUTION SENSITIVITY
//    ///////////////////////////////////////////////////////////////////
//    // clear flutter solver and set the output file
//    
//    MAST::FSIGeneralizedAeroForceAssembly fsi_assembly;
//    fsi_assembly.set_discipline_and_system(fsi_assembly,
//                                              *_structural_discipline,
//                                              *_structural_sys_init);
//    
//    std::ostringstream oss;
//    oss << "flutter_output_" << this->comm().rank() << ".txt";
//    if (this->comm().rank() == 0)
//        _flutter_solver->set_output_file(oss.str());
//    
//    
//    fsi_assembly.init(&solver,                       // fluid complex solver
//                      _pressure_function,
//                      _freq_domain_pressure_function,
//                      _displ);
//    _flutter_solver->attach_assembly(fsi_assembly);
//    
//    // flutter solver will need velocity to be defined as a parameter for
//    // sensitivity analysis
//    // calculate the sensitivity
//    _flutter_solver->calculate_sensitivity(*_flutter_root, p);
//    
//    
//    // clean up before exiting
//    fsi_assembly.clear_discipline_and_system();
//    _flutter_solver->clear_assembly_object();
//    
//    
//    
//    return _flutter_root->V_sens;
//}
//
