/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include <vector>

// MAST includes
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "examples/fluid/meshing/mesh_initializer.h"
#include "examples/base/augment_ghost_elem_send_list.h"
#include "examples/base/plot_results.h"
#include "examples/base/input_wrapper.h"
#include "constrain_beam_dofs.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/complex_mesh_field_function.h"
#include "base/complex_assembly_base.h"
#include "base/eigenproblem_assembly.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/flight_condition.h"
#include "fluid/pressure_function.h"
#include "fluid/frequency_domain_pressure_function.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/ug_flutter_solver.h"
#include "aeroelasticity/ug_flutter_root.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/complex_normal_rotation_mesh_function.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/fluid_structure_assembly_elem_operations.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/isotropic_material_property_card.h"
#include "solver/slepc_eigen_solver.h"
#include "solver/complex_solver_base.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_elem_type.h"    // ElemType
#include "libmesh/fe_type.h"           // FEFamily, Order
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"           
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nonlinear_solver.h"

int main(int argc, char* argv[]) {

    // BEGIN_TRANSLATE Panel flutter using V-g solver and Euler aerodynamics
    //
    //   \tableofcontents
    //
    // This example solves for the flutter speed of a semi-infinite (2D)
    // panel with
    //
    // Initialize libMesh library.
    libMesh::LibMeshInit init(argc, argv);

    // initialize the wrapper to read input parameters. This will check
    // if the executable parameters included a parameter of type
    // `input=${filename}`. If included, then the input parameters will be read
    // from this filename. Otherwise, the parameters will be read from the
    // executable arguments. The wrapper uses default values for parameters
    // if none are provided.
    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    const unsigned int
    dim                 = 2,
    nx_divs             = 3,
    ny_divs             = 1,
    panel_bc_id         = 4,     // this is the id of the boundary set for panel
    symmetry_bc_id      = 5,     // this is the id of the boundary set for remaining region on the bottom
    n_modes             = input("n_modes", "number of structural modes to use for creation of reduced order flutter eigenproblem", 8),
    n_divs_ff_to_panel  = input("n_divs_farfield_to_panel", "number of element divisions from far-field to panel", 30),
    n_divs_panel        = input("n_divs_panel", "number of element divisions on panel", 10),
    n_k_divs            = input("n_k_divs", "number of divisions between upper and lower reduced frequencies for search of flutter root", 10),
    max_bisection_iters = input("flutter_max_bisection_it", "maximum number of bisection iterations to search flutter root after initial sweep of reduced frequencies", 10);
    
    const Real
    k_upper             = input("k_upper", "upper value of reduced frequency for flutter solver", 0.75),
    k_lower             = input("k_lower", "lower value of reduced frequency for flutter solver", 0.00),
    tol                 = input(  "g_tol", "tolerace for convergence of damping of flutter root", 1e-4),
    length              = input("panel_l",                                     "length of panel",  0.3),
    ff_to_panel_l       = input("farfield_to_l_ratio", "Ratio of distance of farfield boundary to panel length",  5.0),
    ff_to_panel_e_size  = input("farfield_to_panel_elem_size_ratio", "Ratio of element size at far-field to element size at panel",  20.0);
    
    
    //////////////////////////////////////////////////////////////////////
    //  \section fluid_init Initialize Fluid Solver
    //////////////////////////////////////////////////////////////////////
    std::string s;
    s                   = input("fe_order", "order of finite element shape basis functions",     "first");
    libMesh::Order
    fe_order            = libMesh::Utility::string_to_enum<libMesh::Order>(s);
    s                   = input("fluid_elem_type",  "type of geometric element in the fluid mesh",     "quad4");
    libMesh::ElemType
    elem_type           = libMesh::Utility::string_to_enum<libMesh::ElemType>(s);
    s                   = input("fe_family",      "family of finite element shape functions", "lagrange");
    libMesh::FEFamily
    fe_family           = libMesh::Utility::string_to_enum<libMesh::FEFamily>(s);
    
    // setting up fluid mesh
    std::vector<Real>
    x_div_loc           = {-length*ff_to_panel_l, 0., length, length*(1.+ff_to_panel_l)},
    x_relative_dx       = {ff_to_panel_e_size, 1., 1., ff_to_panel_e_size},
    y_div_loc           = {0., length*ff_to_panel_l},
    y_relative_dx       = {1., ff_to_panel_e_size};
    
    std::vector<unsigned int>
    x_divs              = {n_divs_ff_to_panel, n_divs_panel, n_divs_ff_to_panel},
    y_divs              = {n_divs_ff_to_panel};
    
    MAST::MeshInitializer::CoordinateDivisions
    x_coord_divs,
    y_coord_divs;
    
    x_coord_divs.init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    y_coord_divs.init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    std::vector<MAST::MeshInitializer::CoordinateDivisions*>
    divs = {&x_coord_divs, &y_coord_divs};
    
    // initialize the fluid mesh
    libMesh::ParallelMesh
    fluid_mesh(init.comm());
    
    // initialize the mesh
    MAST::PanelMesh2D().init(0.,               // t/c
                             false,            // if cos bump
                             0,                // n max bumps
                             divs,
                             fluid_mesh,
                             elem_type);
    
    // initialize equation system
    libMesh::EquationSystems
    fluid_eq_sys(fluid_mesh);
    
    // add the system to be used for fluid analysis
    MAST::NonlinearSystem
    &fluid_sys = fluid_eq_sys.add_system<MAST::NonlinearSystem>("fluid");
    
    // fluid properties, boundary conditions
    MAST::ConservativeFluidDiscipline
    fluid_discipline(fluid_eq_sys);
    
    // variables
    MAST::ConservativeFluidSystemInitialization
    fluid_sys_init(fluid_sys,
                   fluid_sys.name(),
                   libMesh::FEType(fe_order, fe_family), dim);
    
    // keep fluid boundary on all partitions
    MAST::AugmentGhostElementSendListObj augment_send_list_obj(fluid_sys);
    fluid_sys.get_dof_map().attach_extra_send_list_object(augment_send_list_obj);
    
    fluid_eq_sys.init();
    
    // print the information
    fluid_mesh.print_info();
    fluid_eq_sys.print_info();
    
    // create the boundary conditions for slip-wall and far-field
    MAST::BoundaryConditionBase
    far_field(MAST::FAR_FIELD),
    symm_wall(MAST::SYMMETRY_WALL),
    slip_wall(MAST::SLIP_WALL);
    
    // tell the physics about boundary conditions
    fluid_discipline.add_side_load(   panel_bc_id, slip_wall);
    fluid_discipline.add_side_load(symmetry_bc_id, symm_wall);
    for (unsigned int i=1; i<=3; i++)
        fluid_discipline.add_side_load(i, far_field);
    
    // set fluid properties
    MAST::FlightCondition flight_cond;
    flight_cond.flow_unit_vector << 1, 0, 0;
    flight_cond.ref_chord        = length;
    flight_cond.mach             = input("fluid_mach", "fluid Mach number",                           0.5);
    flight_cond.gas_property.cp  = input("fluid_cp",   "fluid specific heat at constant pressure",  1003.);
    flight_cond.gas_property.cv  = input("fluid_cv",   "fluid specific heat at constant volume",     716.);
    flight_cond.gas_property.T   = input("fluid_T",    "fluid absolute temperature",                 300.);
    flight_cond.gas_property.rho = input("fluid_rho",  "fluid density",                              1.35);
    flight_cond.init();
    
    // tell the discipline about the fluid values
    fluid_discipline.set_flight_condition(flight_cond);
    
    // define parameters
    MAST::Parameter
    omega    (   "omega", 0.),
    velocity ("velocity", flight_cond.velocity_magnitude()),
    b_ref    (   "b_ref", length);
    
    
    // now define the constant field functions based on this
    MAST::ConstantFieldFunction
    omega_f    (   "omega", omega),
    velocity_f ("velocity", velocity),
    b_ref_f    (   "b_ref", b_ref);
    
    // TODO: replace with interpolation
    MAST::PressureFunction
    pressure_function(fluid_sys_init, flight_cond);
    pressure_function.set_calculate_cp(true);
    
    MAST::FrequencyDomainPressureFunction
    freq_domain_pressure_function(fluid_sys_init, flight_cond);
    freq_domain_pressure_function.set_calculate_cp(true);
    
    //////////////////////////////////////////////////////////////////////
    //  \section structural_init Initialize Structural Solver
    //////////////////////////////////////////////////////////////////////
    
    x_div_loc = {0.0, length};
    x_relative_dx = {1., 1.};
    x_divs = {n_divs_panel};
    
    MAST::MeshInitializer::CoordinateDivisions
    x_coord_divs_struct;
    x_coord_divs_struct.init(1, x_div_loc, x_relative_dx, x_divs);
    
    divs = {&x_coord_divs_struct};
    
    // create the mesh
    libMesh::SerialMesh structural_mesh(init.comm());
    
    if (elem_type == libMesh::QUAD4)
        MAST::MeshInitializer().init(divs, structural_mesh, libMesh::EDGE2);
    else if (elem_type == libMesh::QUAD8 || elem_type == libMesh::QUAD9)
        MAST::MeshInitializer().init(divs, structural_mesh, libMesh::EDGE3);
    
    // create the equation system
    libMesh::EquationSystems structural_eq_sys(structural_mesh);
    
    // create the libmesh system
    MAST::NonlinearSystem
    &structural_sys = structural_eq_sys.add_system<MAST::NonlinearSystem>("structural");
    structural_sys.set_eigenproblem_type(libMesh::GHEP);
    
    // initialize the system to the right set of variables
    MAST::PhysicsDisciplineBase
    structural_discipline(structural_eq_sys);
    
    MAST::StructuralSystemInitialization
    structural_sys_init(structural_sys,
                        structural_sys.name(),
                        libMesh::FEType(fe_order, fe_family));
    
    // create and add the boundary condition and loads
    MAST::DirichletBoundaryCondition
    dirichlet_left,
    dirichlet_right;
    std::vector<unsigned int> constrained_vars = {0, 1, 2, 3}; // u, v, w, tx
    dirichlet_left.init (0, constrained_vars);
    dirichlet_right.init(1, constrained_vars);
    
    structural_discipline.add_dirichlet_bc(0, dirichlet_left);
    structural_discipline.add_dirichlet_bc(1, dirichlet_right);
    structural_discipline.init_system_dirichlet_bc(structural_sys);
    
    MAST::ConstrainBeamDofs constraint_beam_dofs(structural_sys);
    structural_sys.attach_constraint_object(constraint_beam_dofs);
    
    // initialize the equation system
    structural_eq_sys.init();
    
    // print information
    structural_mesh.print_info();
    structural_eq_sys.print_info();
    
    // TODO: replace this functionality in conservative_fluid_element
    MAST::ComplexMeshFieldFunction
    displ(structural_sys_init, "frequency_domain_displacement");
    
    MAST::ComplexNormalRotationMeshFunction
    normal_rot("frequency_domain_normal_rotation", displ);
    
    slip_wall.add(displ);
    slip_wall.add(normal_rot);
    
    // set up structural eigenvalue problem
    structural_sys.eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    structural_sys.set_exchange_A_and_B(true);
    structural_sys.set_n_requested_eigenvalues(n_modes);
    
    // create the property functions and add them to the material property card
    MAST::Parameter
    thy        ("thy", input("str_thickness",  "thickness of beam",            0.0015)),
    thz        ("thz", 1.00),
    rho        ("rho", input("str_rho",  "structural material density",         2.7e3)),
    E          ("E",   input("str_E",    "structural material Young's modulus", 72.e9)),
    nu         ("nu",  input("str_nu",   "structural material Poisson's ratio",  0.33)),
    kappa_yy   ("kappa_yy", 5./6.),
    kappa_zz   ("kappa_zz", 5./6.),
    zero       ("zero", 0.);
    
    MAST::ConstantFieldFunction
    thy_f      ("hy", thy),
    thz_f      ("hz", thz),
    rho_f      ("rho", rho),
    E_f        ("E", E),
    nu_f       ("nu", nu),
    kappa_yy_f ("Kappayy", kappa_yy),
    kappa_zz_f ("Kappazz", kappa_zz),
    hyoff_f    ("hy_off", zero),
    hzoff_f    ("hz_off", zero);
    
    MAST::IsotropicMaterialPropertyCard
    m_card;
    m_card.add(rho_f);
    m_card.add(E_f);
    m_card.add(nu_f);
    
    MAST::Solid1DSectionElementPropertyCard
    p_card;
    p_card.set_bending_model(MAST::TIMOSHENKO);
    
    // tell the card about the orientation
    p_card.y_vector()    = RealVectorX::Zero(3);
    p_card.y_vector()(1) = 1.;
    
    // add the section properties to the card
    p_card.add(thy_f);
    p_card.add(thz_f);
    p_card.add(hyoff_f);
    p_card.add(hzoff_f);
    p_card.add(kappa_yy_f);
    p_card.add(kappa_zz_f);
    
    // tell the section property about the material property
    p_card.set_material(m_card);
    
    p_card.init();
    
    structural_discipline.set_property_for_subdomain(0, p_card);
    
    // pressure boundary condition for the beam
    MAST::BoundaryConditionBase pressure(MAST::SURFACE_PRESSURE);
    pressure.add(pressure_function);
    pressure.add(freq_domain_pressure_function);
    structural_discipline.add_volume_load(0, pressure);
    
    //////////////////////////////////////////////////////////////////////
    //  \section flutter_init Initialize Flutter Solver
    //////////////////////////////////////////////////////////////////////
    
    // initialize the frequency function
    MAST::FrequencyFunction
    freq_function("freq", omega_f, velocity_f, b_ref_f);
    freq_function.if_nondimensional(true);
    
    /////////////////////////////////////////////////////////////////
    //  INITIALIZE FLUID SOLUTION
    /////////////////////////////////////////////////////////////////
    // the modal and flutter problems are solved on rank 0, while
    // the fluid solution is setup on the global communicator
    
    // initialize the solution
    RealVectorX fluid_ff_vars(4);
    fluid_ff_vars << flight_cond.rho(),
    flight_cond.rho_u1(),
    flight_cond.rho_u2(),
    flight_cond.rho_e();
    
    // create the vector for storing the base solution.
    // we will swap this out with the system solution, initialize and
    // then swap it back.
    libMesh::NumericVector<Real>& base_sol =
    fluid_sys.add_vector("fluid_base_solution");
    fluid_sys.solution->swap(base_sol);
    fluid_sys_init.initialize_solution(fluid_ff_vars);
    fluid_sys.solution->swap(base_sol);
    
    // create the nonlinear assembly object
    MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations fluid_elem_ops;
    MAST::ComplexAssemblyBase                                    complex_assembly;
    
    // now set up the assembly objects
    fluid_elem_ops.set_discipline_and_system(fluid_discipline, fluid_sys_init);
    complex_assembly.set_discipline_and_system(fluid_discipline, fluid_sys_init);
    complex_assembly.set_base_solution(base_sol);
    fluid_elem_ops.set_frequency_function(freq_function);
    
    pressure_function.init(base_sol);
    
    ////////////////////////////////////////////////////////////
    //  \section solution Solution
    //  \subsection modal_sol Structural Modal Solution
    ////////////////////////////////////////////////////////////
    
    MAST::EigenproblemAssembly modal_assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations modal_elem_ops;
    
    modal_assembly.set_discipline_and_system(structural_discipline, structural_sys_init);
    modal_elem_ops.set_discipline_and_system(structural_discipline, structural_sys_init);
    
    structural_sys.initialize_condensed_dofs(structural_discipline);
    
    MAST::StructuralNearNullVectorSpace nsp;
    structural_sys.nonlinear_solver->nearnullspace_object = &nsp;
    
    structural_sys.eigenproblem_solve(modal_elem_ops, modal_assembly);
    modal_assembly.clear_discipline_and_system();
    modal_elem_ops.clear_discipline_and_system();
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(structural_sys.get_n_converged_eigenvalues(),
                     structural_sys.get_n_requested_eigenvalues());
    
    std::vector<libMesh::NumericVector<Real>*> basis(nconv, nullptr);
    
    libMesh::ExodusII_IO writer(structural_mesh);
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // create a vector to store the basis
        basis[i] = structural_sys.solution->zero_clone().release();
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        structural_sys.get_eigenpair(i, re, im, *basis[i]);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        // We write the file in the ExodusII format.
        // copy the solution for output
        structural_sys.solution->swap(*basis[i]);
        writer.write_timestep("modes.exo",
                              structural_eq_sys,
                              i+1, i);
        structural_sys.solution->swap(*basis[i]);
    }
    
    ///////////////////////////////////////////////////////////////////
    //  \subsection flutter_sol Flutter Solution
    ///////////////////////////////////////////////////////////////////
    
    MAST::UGFlutterSolver flutter_solver;
    
    // solver for complex solution
    MAST::ComplexSolverBase solver;
    
    MAST::FSIGeneralizedAeroForceAssembly      fsi_assembly;
    MAST::FluidStructureAssemblyElemOperations fsi_elem_ops;
    fsi_assembly.set_discipline_and_system(structural_discipline, structural_sys_init);
    fsi_elem_ops.set_discipline_and_system(structural_discipline, structural_sys_init);

    fsi_assembly.init(fsi_elem_ops,
                      solver,                       // fluid complex solver
                      complex_assembly,
                      fluid_elem_ops,
                      pressure_function,
                      freq_domain_pressure_function,
                      displ);
    
    flutter_solver.attach_assembly(fsi_assembly);
    flutter_solver.initialize(omega,
                              b_ref,
                              flight_cond.rho(),
                              k_lower,            // lower kr
                              k_upper,            // upper kr
                              n_k_divs,           // number of divisions
                              basis);             // basis vectors
    
    std::ostringstream oss;
    oss << "flutter_output_" << init.comm().rank() << ".txt";
    if (init.comm().rank() == 0)
        flutter_solver.set_output_file(oss.str());
    
    // find the roots for the specified divisions
    flutter_solver.scan_for_roots();
    flutter_solver.print_crossover_points();
    
    // now ask the flutter solver to return the critical flutter root,
    // which is the flutter cross-over point at the lowest velocity
    std::pair<bool, MAST::FlutterRootBase*>
    sol = flutter_solver.find_critical_root(tol, max_bisection_iters);
    
    flutter_solver.print_sorted_roots();
    
    // make sure solution was found
    libmesh_assert(sol.first);
    
    //  \subsection plot_flutter Plot Flutter Solution
    MAST::plot_structural_flutter_solution("structural_flutter_mode.exo",
                                           structural_sys,
                                           sol.second->eig_vec_right,
                                           basis);
    
    MAST::plot_fluid_flutter_solution("fluid_flutter_mode.exo",
                                      structural_sys,
                                      fluid_sys,
                                      displ,
                                      solver,
                                      sol.second->eig_vec_right,
                                      basis);
    
    // cleanup the data structures
    fsi_assembly.clear_discipline_and_system();
    fsi_elem_ops.clear_discipline_and_system();
    flutter_solver.clear_assembly_object();
    complex_assembly.clear_discipline_and_system();
    
    // delete the basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        if (basis[i]) delete basis[i];
    
    // END_TRANSLATE
    return 0;
}
