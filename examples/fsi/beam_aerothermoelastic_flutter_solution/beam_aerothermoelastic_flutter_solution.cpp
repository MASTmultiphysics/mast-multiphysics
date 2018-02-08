///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2018  Manav Bhatia
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
//#include "examples/fsi/beam_aerothermoelastic_flutter_solution/beam_aerothermoelastic_flutter_solution.h"
//#include "examples/base/plot_results.h"
//#include "base/nonlinear_system.h"
//#include "base/complex_assembly_base.h"
//#include "base/parameter.h"
//#include "base/constant_field_function.h"
//#include "base/boundary_condition_base.h"
//#include "base/complex_mesh_field_function.h"
//#include "aeroelasticity/ug_flutter_solver.h"
//#include "aeroelasticity/frequency_function.h"
//#include "aeroelasticity/flutter_root_base.h"
//#include "base/eigenproblem_assembly.h"
//#include "elasticity/structural_modal_eigenproblem_assembly.h"
//#include "elasticity/structural_near_null_vector_space.h"
//#include "elasticity/structural_system_initialization.h"
//#include "elasticity/fsi_generalized_aero_force_assembly.h"
//#include "elasticity/complex_normal_rotation_mesh_function.h"
//#include "solver/complex_solver_base.h"
//#include "fluid/frequency_domain_linearized_complex_assembly.h"
//#include "fluid/frequency_domain_pressure_function.h"
//#include "fluid/flight_condition.h"
//#include "fluid/conservative_fluid_discipline.h"
//#include "fluid/conservative_fluid_system_initialization.h"
//#include "solver/slepc_eigen_solver.h"
//
//
//// libMesh includes
//#include "libmesh/numeric_vector.h"
//#include "libmesh/nonlinear_solver.h"
//#include "libmesh/getpot.h"
//#include "libmesh/exodusII_io.h"
//
//
//
//
//
//MAST::BeamAerothermoelasticFlutterSolution::
//BeamAerothermoelasticFlutterSolution(const libMesh::Parallel::Communicator& comm_in):
//MAST::BeamEulerFSIAnalysis                    (comm_in),
//_k_lower                                      (0.),
//_k_upper                                      (0.),
//_n_k_divs                                     (0),
//_displ_perturb                                (nullptr),
//_normal_rot_perturb                           (nullptr),
//_freq_domain_pressure_function                (nullptr),
//_omega                                        (nullptr),
//_b_ref                                        (nullptr),
//_zero_omega                                   (nullptr),
//_omega_f                                      (nullptr),
//_b_ref_f                                      (nullptr),
//_zero_omega_f                                 (nullptr),
//_zero_freq_func                               (nullptr),
//_freq_function                                (nullptr),
//_flutter_solver                               (nullptr),
//_flutter_root                                 (nullptr){
// 
//
//    GetPot infile("input.in");
//    
//    
//    _k_upper            = infile("k_upper",  0.75);
//    _k_lower            = infile("k_lower",  0.05);
//    _n_k_divs           = infile("n_k_divs",   10);
//
//    
//    _omega             = new MAST::Parameter("omega",       0.);
//    _zero_omega        = new MAST::Parameter("omega",       0.);
//    _b_ref             = new MAST::Parameter("b_ref",       1.);
//    
//    
//    // now define the constant field functions based on this
//    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
//    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
//    _zero_omega_f      = new MAST::ConstantFieldFunction("omega",  *_zero_omega);
//    
//    // initialize the frequency function
//    _zero_freq_func    = new MAST::FrequencyFunction("freq",
//                                                     *_zero_omega_f,
//                                                     *_velocity_f,
//                                                     *_b_ref_f);
//    
//    _freq_function     = new MAST::FrequencyFunction("freq",
//                                                     *_omega_f,
//                                                     *_velocity_f,
//                                                     *_b_ref_f);
//    _freq_function->if_nondimensional(true);
//
//    _freq_domain_pressure_function =
//    new MAST::FrequencyDomainPressureFunction(*_fluid_sys_init, *_flight_cond);
//    
//    _displ_perturb      = new MAST::ComplexMeshFieldFunction(*_structural_sys_init,
//                                                             "frequency_domain_displacement");
//    _normal_rot_perturb = new MAST::ComplexNormalRotationMeshFunction("frequency_domain_normal_rotation",
//                                                                      *_displ_perturb);
//    _slip_wall->add(*_displ_perturb);
//    _slip_wall->add(*_normal_rot_perturb);
//    
//    
//    // define parameters
//    _omega             = new MAST::Parameter("omega",       0.);
//    _zero_omega        = new MAST::Parameter("omega",       0.);
//    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
//    _b_ref             = new MAST::Parameter("b_ref",       1.);
//    
//    
//    // now define the constant field functions based on this
//    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
//    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
//    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
//    _zero_omega_f      = new MAST::ConstantFieldFunction("omega",  *_zero_omega);
//
//    
//    _structural_sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
//    _structural_sys->set_exchange_A_and_B(true);
//    _structural_sys->set_n_requested_eigenvalues(infile("n_modes", 10));
//
//    
//    // initialize the frequency function
//    _zero_freq_func    = new MAST::FrequencyFunction("freq",
//                                                     *_zero_omega_f,
//                                                     *_velocity_f,
//                                                     *_b_ref_f);
//    
//    _freq_function     = new MAST::FrequencyFunction("freq",
//                                                     *_omega_f,
//                                                     *_velocity_f,
//                                                     *_b_ref_f);
//    _freq_function->if_nondimensional(true);
//    
//
//    _flutter_solver  = new MAST::UGFlutterSolver;
//
//}
//
//
//
//
//MAST::BeamAerothermoelasticFlutterSolution::~BeamAerothermoelasticFlutterSolution() {
// 
//    delete _displ_perturb;
//    delete _normal_rot_perturb;
//    delete _freq_domain_pressure_function;
//    delete _omega;
//    delete _b_ref;
//    delete _zero_omega;
//    delete _omega_f;
//    delete _b_ref_f;
//    delete _zero_omega_f;
//    delete _zero_freq_func;
//    delete _freq_function;
//    delete _flutter_solver;
//    
//    // delete the basis vectors
//    if (_basis.size())
//        for (unsigned int i=0; i<_basis.size(); i++)
//            if (_basis[i]) delete _basis[i];
//}
//
//
//
//
//
//
//Real
//MAST::BeamAerothermoelasticFlutterSolution::solve(bool if_write_output,
//                                                  const Real tol,
//                                                  const unsigned int max_bisection_iters) {
//    
//    // first solve for the FSI steady-state solution
//    MAST::BeamEulerFSIAnalysis::solve(if_write_output,
//                                      tol,
//                                      max_bisection_iters);
//    
//    // swap back the calculate solution so that we store it for later use.
//    libMesh::NumericVector<Real>
//    &fluid_base_sol      = _fluid_sys->add_vector("fluid_base_solution"),
//    &structural_base_sol = _structural_sys->add_vector("structural_base_solution");
//    _fluid_sys->solution->close();
//    _structural_sys->solution->close();
//    fluid_base_sol       = *_fluid_sys->solution;
//    structural_base_sol  = *_structural_sys->solution;
//    
//    
//    // create the nonlinear assembly object
//    MAST::ComplexAssemblyBase                                      complex_fluid_assembly;
//    MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations   complex_elem_ops;
//
//    // solver for complex solution
//    MAST::ComplexSolverBase                          complex_fluid_solver;
//    
//    // now setup the assembly object
//    complex_fluid_assembly.set_discipline_and_system(complex_elem_ops,
//                                                        *_fluid_discipline,
//                                                        complex_fluid_solver,
//                                                        *_fluid_sys_init);
//    complex_fluid_assembly.set_base_solution(fluid_base_sol);
//    complex_elem_ops.set_frequency_function(*_freq_function);
//
//    
//    ////////////////////////////////////////////////////////////
//    // STRUCTURAL MODAL EIGENSOLUTION
//    ////////////////////////////////////////////////////////////
//    
//    // create the nonlinear assembly object
//    MAST::EigenproblemAssembly   modal_assembly;
//    MAST::StructuralModalEigenproblemAssemblyElemOperations   modal_elem_ops;
//    _structural_sys->initialize_condensed_dofs(*_structural_discipline);
//    
//    modal_assembly.set_discipline_and_system(modal_elem_ops,
//                                                *_structural_discipline,
//                                                *_structural_sys_init);
//    modal_assembly.set_base_solution(structural_base_sol);
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
//    
//    // now, solve a linearized stability problem about this steady-state
//    // solution
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
//    fsi_assembly.init(&complex_fluid_solver,         // fluid complex solver
//                      _pressure_function,
//                      _freq_domain_pressure_function,
//                      _displ_perturb);
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
//                                          *_displ_perturb,
//                                          complex_fluid_solver,
//                                          sol.second->eig_vec_right,
//                                          _basis);
//    }
//    
//    complex_fluid_assembly.clear_discipline_and_system();
//
//    return _flutter_root->V;
//}
//
//
//
