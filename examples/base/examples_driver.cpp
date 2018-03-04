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

// MAST includes
#include "examples/base/input_wrapper.h"
#include "examples/structural/bar_extension/bar_extension.h"
//#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
//#include "examples/structural/beam_buckling_prestress/beam_column_buckling.h"
#include "examples/structural/beam_bending/beam_bending.h"
#include "examples/structural/beam_oscillating_load/beam_oscillating_load.h"
//#include "examples/structural/beam_bending_with_offset/beam_bending_with_offset.h"
//#include "examples/structural/beam_bending_thermal_stress_with_offset/beam_bending_thermal_stress.h"
//#include "examples/structural/beam_optimization/beam_optimization.h"
//#include "examples/structural/beam_optimization_single_stress_functional/beam_optimization.h"
//#include "examples/structural/beam_optimization_section_offset/beam_optimization_section_offset.h"
//#include "examples/structural/beam_optimization_thermal_stress/beam_optimization_thermal_stress.h"
//#include "examples/structural/beam_piston_theory_flutter/beam_piston_theory_flutter.h"
//#include "examples/structural/beam_piston_theory_time_accurate/beam_piston_theory_time_accurate.h"
//#include "examples/structural/membrane_extension_uniaxial_stress/membrane_extension_uniaxial.h"
//#include "examples/structural/membrane_extension_biaxial_stress/membrane_extension_biaxial.h"
#include "examples/structural/plate_bending/plate_bending.h"
//#include "examples/structural/plate_bending_level_set/plate_bending_level_set.h"
//#include "examples/structural/stiffened_plate_bending_thermal_stress/stiffened_plate_bending_thermal_stress.h"
//#include "examples/structural/plate_oscillating_load/plate_oscillating_load.h"
//#include "examples/structural/plate_bending_section_offset/plate_bending_section_offset.h"
//#include "examples/structural/plate_bending_thermal_stress/plate_bending_thermal_stress.h"
//#include "examples/structural/plate_modal_analysis/plate_modal_analysis.h"
//#include "examples/structural/plate_thermally_stressed_modal_analysis/plate_thermally_stressed_modal_analysis.h"
//#include "examples/structural/nastran_model_analysis/nastran_model_analysis.h"
//#include "examples/structural/plate_buckling_prestress/plate_buckling_prestress.h"
//#include "examples/structural/plate_piston_theory_flutter/plate_piston_theory_flutter.h"
//#include "examples/structural/plate_thermally_stressed_piston_theory_flutter/plate_thermally_stressed_piston_theory_flutter.h"
//#include "examples/structural/plate_optimization/plate_optimization.h"
//#include "examples/structural/plate_optimization_single_stress_functional/plate_optimization_single_functional.h"
//#include "examples/structural/plate_optimization_section_offset/plate_section_offset_optimization.h"
//#include "examples/structural/plate_optimization_thermal_stress/plate_thermal_stress_optimization.h"
//#include "examples/structural/stiffened_plate_optimization_thermal_stress/stiffened_plate_thermal_stress_optimization.h"
//#include "examples/structural/stiffened_plate_optimization_piston_theory_flutter/stiffened_plate_piston_theory_flutter_optimization.h"
//#include "examples/structural/stiffened_plate_optimization_thermally_stressed_piston_theory_flutter/stiffened_plate_thermally_stressed_piston_theory_flutter_optimization.h"
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"
//#include "optimization/npsol_optimization_interface.h"
#include "optimization/dot_optimization_interface.h"
#include "optimization/gcmma_optimization_interface.h"
//#include "examples/fluid/panel_inviscid_analysis_2D/panel_inviscid_analysis_2d.h"
//#include "examples/fluid/ramp_laminar_analysis_2D/ramp_viscous_analysis_2d.h"
//#include "examples/fluid/panel_inviscid_analysis_3D_half_domain/panel_inviscid_analysis_3D_half_domain.h"
//#include "examples/fluid/panel_small_disturbance_frequency_domain_analysis_2D/panel_small_disturbance_frequency_domain_analysis_2d.h"
//#include "examples/fluid/panel_small_disturbance_frequency_domain_3D/panel_small_disturbance_frequency_domain_inviscid_analysis_3D.h"
//#include "examples/fluid/panel_small_disturbance_frequency_domain_3D_half_domain/panel_small_disturbance_frequency_domain_inviscid_analysis_3D_half_domain.h"
//#include "examples/fsi/beam_flutter_solution/beam_euler_fsi_flutter_solution.h"
//#include "examples/fsi/beam_flag_flutter_solution/beam_flag_euler_fsi_flutter_solution.h"
//#include "examples/fsi/beam_fsi_flutter_high_order_convergence/beam_fsi_high_order_convergence.h"
//#include "examples/fsi/beam_flutter_nonuniform_aero_base_solution/beam_euler_fsi_flutter_nonuniform_aero_base_solution.h"
//#include "examples/fsi/beam_flutter_optimization/beam_flutter_optimization.h"
//#include "examples/fsi/plate_flutter_solution/plate_euler_fsi_flutter_solution.h"
//#include "examples/fsi/plate_flag_flutter_solution/plate_flag_euler_fsi_flutter_solution.h"
//#include "examples/fsi/plate_flutter_optimization/plate_flutter_optimization.h"
//#include "examples/fsi/stiffened_plate_thermally_stressed_flutter_optimization/stiffened_plate_thermally_stressed_flutter_optimization.h"
//#include "examples/fsi/plate_flutter_solution_half_domain/plate_euler_fsi_half_domain_flutter_solution.h"
//#include "examples/fsi/beam_fsi_solution/beam_euler_fsi_solution.h"
////#include "examples/fsi/beam_aerothermoelastic_flutter_solution/beam_aerothermoelastic_flutter_solution.h"
#include "examples/thermal/base/thermal_example_1d.h"
#include "examples/thermal/base/thermal_example_2d.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"




//template <typename ValType>
//void flutter_analysis(const std::string& case_name,
//                      libMesh::ElemType etype,
//                      const bool nonlinear,
//                      bool with_sens,
//                      std::string& par_name)  {
//
//    ValType run_case;
//    run_case.init(etype, nonlinear);
//
//    libMesh::out << "Running case: " << case_name << std::endl;
//    run_case.solve(true);
//    if (with_sens) {
//        MAST::Parameter* p = run_case.get_parameter(par_name);
//        if (p) {
//
//            std::vector<Real> eig;
//            libMesh::out
//            << "Running sensitivity for case: " << case_name
//            << "  wrt  " << par_name << std::endl;
//            run_case.sensitivity_solve(*p);
//        }
//    }
//
//}



//template <typename ValType>
//void fluid_analysis(const std::string& case_name)  {
//
//    ValType run_case;
//
//    libMesh::out << "Running case: " << case_name << std::endl;
//    run_case.solve(true);
//    /*if (with_sens) {
//     MAST::Parameter* p = run_case.get_parameter(par_name);
//     if (p) {
//
//     libMesh::out
//     << "Running sensitivity for case: " << case_name
//     << "  wrt  " << par_name << std::endl;
//     run_case.sensitivity_solve(*p, true);
//     }
//     }*/
//}


int main(int argc, char* const argv[]) {
    
    libMesh::LibMeshInit init(argc, argv);
    
    // use to get arguments from the command line
    GetPot command_line(argc, argv);
    
    // look for the name that the user has requested to run
    std::string
    input_name   = command_line(   "input",    ""),
    prefix       = command_line(  "prefix",    ""),
    case_name    = command_line("run_case",    ""),
    par_name     = command_line(   "param",    "");
    bool
    print_param  = command_line("print_params",        false),
    with_sens    = command_line("with_sensitivity",    false),
    if_nonlin    = command_line("nonlinear",           false),
    verify_grads = command_line("verify_grads",        false);
    
    // initialize the input file for this case
    std::unique_ptr<MAST::Examples::GetPotWrapper> input;
    if (input_name == "")
        input.reset(new MAST::Examples::GetPotWrapper(argc, argv));
    else
        input.reset(new MAST::Examples::GetPotWrapper(input_name));
    input->set_print(print_param);
    
    if (case_name == "bar_extension") {
        
        MAST::Examples::BarExtension example(init.comm());
        example.init(*input, prefix);
        example.static_solve();
        example.static_sensitivity_solve(example.get_parameter("thy"));
        example.static_adjoint_sensitivity_solve(example.get_parameter("thy"));
    }
//    else if (case_name == "beam_modal_analysis")
//        eigenvalue_analysis<MAST::BeamModalAnalysis>(case_name,
//                                                     libMesh::EDGE2,
//                                                     false,
//                                                     with_sens,
//                                                     par_name);
//    else if (case_name == "beam_thermally_stressed_modal_analysis")
//        eigenvalue_analysis<MAST::BeamThermallyStressedModalAnalysis>(case_name,
//                                                                      libMesh::EDGE2,
//                                                                      if_nonlin,
//                                                                      with_sens,
//                                                                      par_name);
//    else if (case_name == "beam_prestress_buckling_analysis")
//        eigenvalue_analysis<MAST::BeamColumnBucklingAnalysis>(case_name,
//                                                              libMesh::EDGE2,
//                                                              false,
//                                                              with_sens,
//                                                              par_name);
    else if (case_name == "beam_bending") {
        
        MAST::Examples::BeamBending example(init.comm());
        example.init(*input, prefix);
        example.static_solve();
        example.static_sensitivity_solve(example.get_parameter("thy"));
        example.static_adjoint_sensitivity_solve(example.get_parameter("thy"));
    }
    else if (case_name == "beam_oscillating_load") {
        
        MAST::Examples::BeamOscillatingLoad example(init.comm());
        example.init(*input, prefix);
        example.transient_solve();
        example.transient_sensitivity_solve(example.get_parameter("thy"));
    }
//    else if (case_name == "beam_bending_with_offset")
//        analysis<MAST::BeamBendingWithOffset>(case_name,
//                                              libMesh::EDGE2,
//                                              if_nonlin,
//                                              with_sens,
//                                              par_name);
//    else if (case_name == "beam_bending_thermal_stress")
//        analysis<MAST::BeamBendingThermalStress>(case_name,
//                                                 libMesh::EDGE2,
//                                                 if_nonlin,
//                                                 with_sens,
//                                                 par_name);
//    else if (case_name == "beam_bending_optimization")
//        optimization<MAST::BeamBendingSizingOptimization>(case_name,
//                                                          libMesh::EDGE2,
//                                                          if_nonlin,
//                                                          verify_grads);
//    else if (case_name == "beam_bending_single_functional_optimization")
//        optimization<MAST::BeamBendingSingleFunctionalSizingOptimization>
//        (case_name,
//         libMesh::EDGE2,
//         if_nonlin,
//         verify_grads);
//    else if (case_name == "beam_bending_section_offset_optimization")
//        optimization<MAST::BeamBendingSectionOffsetSizingOptimization>
//        (case_name,
//         libMesh::EDGE2,
//         if_nonlin,
//         verify_grads);
//    else if (case_name == "beam_bending_thermal_stress_optimization")
//        optimization<MAST::BeamBendingThermalStressSizingOptimization>
//        (case_name,
//         libMesh::EDGE2,
//         if_nonlin,
//         verify_grads);
//    else if (case_name == "beam_piston_theory_flutter_analysis")
//        flutter_analysis<MAST::BeamPistonTheoryFlutterAnalysis>(case_name,
//                                                                libMesh::EDGE2,
//                                                                if_nonlin,
//                                                                with_sens,
//                                                                par_name);
//    else if (case_name == "beam_piston_theory_time_accurate_analysis")
//        flutter_analysis<MAST::BeamPistonTheoryTimeAccurateAnalysis>(case_name,
//                                                                     libMesh::EDGE2,
//                                                                     if_nonlin,
//                                                                     with_sens,
//                                                                     par_name);
//    else if (case_name == "membrane_extension_uniaxial")
//        analysis<MAST::MembraneExtensionUniaxial>(case_name,
//                                                  libMesh::EDGE2,
//                                                  if_nonlin,
//                                                  with_sens,
//                                                  par_name);
//    else if (case_name == "membrane_extension_biaxial")
//        analysis<MAST::MembraneExtensionBiaxial>(case_name,
//                                                 libMesh::EDGE2,
//                                                 if_nonlin,
//                                                 with_sens,
//                                                 par_name);
//    else if (case_name == "plate_modal_analysis")
//        eigenvalue_analysis<MAST::PlateModalAnalysis>(case_name,
//                                                      libMesh::QUAD4,
//                                                      if_nonlin,
//                                                      with_sens,
//                                                      par_name);
//    else if (case_name == "plate_thermally_stressed_modal_analysis")
//        eigenvalue_analysis<MAST::PlateThermallyStressedModalAnalysis>(case_name,
//                                                                       libMesh::QUAD4,
//                                                                       if_nonlin,
//                                                                       with_sens,
//                                                                       par_name);
//    else if (case_name == "plate_prestress_buckling_analysis")
//        eigenvalue_analysis<MAST::PlateBucklingPrestress>(case_name,
//                                                          libMesh::QUAD4,
//                                                          true,
//                                                          with_sens,
//                                                          par_name);
    else if (case_name == "plate_bending") {
        
        MAST::Examples::PlateBending example(init.comm());
        example.init(*input, prefix);
        example.static_solve();
        example.static_sensitivity_solve(example.get_parameter("th"));
        example.static_adjoint_sensitivity_solve(example.get_parameter("th"));
    }
//    else if (case_name == "plate_bending_level_set")
//        analysis<MAST::PlateBendingLevelSet>(case_name,
//                                             libMesh::QUAD4,
//                                             if_nonlin,
//                                             with_sens,
//                                             par_name);
//    else if (case_name == "plate_oscillating_load")
//        analysis<MAST::PlateOscillatingLoad>(case_name,
//                                             libMesh::QUAD4,
//                                             if_nonlin,
//                                             with_sens,
//                                             par_name);
//    else if (case_name == "plate_bending_section_offset")
//        analysis<MAST::PlateBendingWithOffset>(case_name,
//                                               libMesh::QUAD4,
//                                               if_nonlin,
//                                               with_sens,
//                                               par_name);
//    else if (case_name == "plate_bending_thermal_stress")
//        analysis<MAST::PlateBendingThermalStress>(case_name,
//                                                  libMesh::QUAD4,
//                                                  if_nonlin,
//                                                  with_sens,
//                                                  par_name);
//    else if (case_name == "stiffened_plate_bending_thermal_stress")
//        analysis<MAST::StiffenedPlateBendingThermalStress>(case_name,
//                                                           libMesh::QUAD4,
//                                                           if_nonlin,
//                                                           with_sens,
//                                                           par_name);
//    else if (case_name == "plate_piston_theory_flutter_analysis")
//        flutter_analysis<MAST::PlatePistonTheoryFlutterAnalysis>(case_name,
//                                                                 libMesh::QUAD4,
//                                                                 false,
//                                                                 with_sens,
//                                                                 par_name);
//    else if (case_name == "plate_thermally_stressed_piston_theory_flutter_analysis")
//        flutter_analysis<MAST::PlateThermallyStressedPistonTheoryFlutterAnalysis>
//        (case_name,
//         libMesh::QUAD4,
//         if_nonlin,
//         with_sens,
//         par_name);
//    else if (case_name == "nastran_model_analysis") {
//#if MAST_ENABLE_CYTHON == 1
//        analysis<MAST::NastranModelAnalysis>(case_name,
//                                             libMesh::INVALID_ELEM,
//                                             if_nonlin,
//                                             with_sens,
//                                             par_name);
//#else
//    libMesh::out
//    << "MAST not configured with Cython. Working with NASTRAN files requires Pynastran and Cython."
//    << std::endl;
//    libmesh_error();
//#endif
//}
//    else if (case_name == "plate_bending_sizing_optimization")
//        optimization<MAST::PlateBendingSizingOptimization>(case_name,
//                                                           libMesh::QUAD4,
//                                                           verify_grads,
//                                                           if_nonlin);
//    else if (case_name == "plate_bending_single_functional_sizing_optimization")
//        optimization<MAST::PlateBendingSingleStressFunctionalSizingOptimization>
//        (case_name,
//         libMesh::QUAD4,
//         verify_grads,
//         if_nonlin);
//    else if (case_name == "plate_bending_section_offset_optimization")
//        optimization<MAST::PlateBendingSectionOffsetSizingOptimization>(case_name,
//                                                                        libMesh::QUAD4,
//                                                                        verify_grads,
//                                                                        if_nonlin);
//    else if (case_name == "plate_bending_thermal_stress_optimization")
//        optimization<MAST::PlateBendingThermalStressSizingOptimization>(case_name,
//                                                                        libMesh::QUAD4,
//                                                                        verify_grads,
//                                                                        if_nonlin);
//    else if (case_name == "stiffened_plate_bending_thermal_stress_optimization")
//        optimization<MAST::StiffenedPlateBendingThermalStressSizingOptimization>
//        (case_name,
//         libMesh::QUAD4,
//         verify_grads,
//         if_nonlin);
//    else if (case_name == "stiffened_plate_piston_theory_optimization")
//        optimization<MAST::StiffenedPlatePistonTheorySizingOptimization>(case_name,
//                                                                         libMesh::QUAD4,
//                                                                         verify_grads,
//                                                                         if_nonlin);
//    else if (case_name == "stiffened_plate_thermally_stressed_piston_theory_optimization")
//        optimization<MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization>
//        (case_name,
//         libMesh::QUAD4,
//         verify_grads,
//         if_nonlin);
    else if (case_name == "topology_optimization_2D") {

        MAST::GCMMAOptimizationInterface optimizer;
        MAST::Examples::TopologyOptimizationLevelSet2D example(init.comm());
        example.init(*input, prefix);
        optimizer.attach_function_evaluation_object(example);
        //std::vector<Real> dvals(example.n_vars()), dummy(example.n_vars());
        //example.init_dvar(dvals, dummy, dummy);
        //example.verify_gradients(dvals);
        //return 0;
        optimizer.optimize();
    }
//    else if (case_name == "panel_inviscid_analysis_2d")
//        fluid_analysis<MAST::PanelInviscidAnalysis2D>(case_name);
//    else if (case_name == "ramp_laminar_analysis_2d")
//        fluid_analysis<MAST::RampLaminarAnalysis2D>(case_name);
//    else if (case_name == "panel_inviscid_analysis_3d_half_domain")
//        fluid_analysis<MAST::PanelInviscidAnalysis3DHalfDomain>(case_name);
//    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_2d")
//        fluid_analysis<MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis>(case_name);
//    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_3d")
//        fluid_analysis<MAST::PanelSmallDisturbanceFrequencyDomainInviscidAnalysis3D>(case_name);
//    else if (case_name == "panel_inviscid_small_disturbance_frequency_domain_analysis_3d_half_domain")
//        fluid_analysis<MAST::PanelSmallDisturbanceFrequencyDomainInviscidAnalysis3DHalfDomain>(case_name);
//    else if (case_name == "beam_fsi_analysis")
//        fluid_analysis<MAST::BeamEulerFSIAnalysis>(case_name);
//    else if (case_name == "beam_fsi_flutter_analysis")
//        fluid_analysis<MAST::BeamEulerFSIFlutterAnalysis>(case_name);
//    else if (case_name == "beam_flag_fsi_flutter_analysis")
//        fluid_analysis<MAST::BeamFlagEulerFSIFlutterAnalysis>(case_name);
//    else if (case_name == "beam_fsi_ho_convergence")
//        fluid_analysis<MAST::BeamFSIFlutterHighOrderConvergence>(case_name);
//    else if (case_name == "beam_fsi_nonuniform_aero_base_flutter_analysis")
//        fluid_analysis<MAST::BeamEulerFSIFlutterNonuniformAeroBaseAnalysis>(case_name);
////    else if (case_name == "beam_aerothermoelastic_flutter_analysis")
////        fluid_analysis<MAST::BeamAerothermoelasticFlutterSolution>(case_name);
//    else if (case_name == "beam_fsi_flutter_optimization")
//        optimization<MAST::BeamFSIFlutterSizingOptimization>(case_name,
//                                                             libMesh::EDGE2,
//                                                             verify_grads,
//                                                             if_nonlin);
//    else if (case_name == "plate_fsi_flutter_analysis")
//        fluid_analysis<MAST::PlateEulerFSIFlutterAnalysis>(case_name);
//    else if (case_name == "plate_flag_fsi_flutter_analysis")
//        fluid_analysis<MAST::PlateFlagEulerFSIFlutterAnalysis>(case_name);
//    else if (case_name == "plate_fsi_half_domain_flutter_analysis")
//        fluid_analysis<MAST::PlateEulerFSIHalfDomainFlutterAnalysis>(case_name);
//    else if (case_name == "plate_fsi_flutter_optimization")
//        optimization<MAST::PlateFSIFlutterSizingOptimization>(case_name,
//                                                              libMesh::QUAD4,
//                                                              verify_grads,
//                                                              if_nonlin);
//    else if (case_name == "stiffened_plate_fsi_thermally_stressed_flutter_optimization")
//        optimization<MAST::StiffenedPlateThermallyStressedFSIFlutterSizingOptimization>
//        (case_name,
//         libMesh::QUAD4,
//         verify_grads,
//         if_nonlin);
    else if (case_name == "steady_state_conduction_1d") {
        
        MAST::Examples::ThermalExample1D example(init.comm());
        example.init(*input, prefix);
        example.steady_solve();
        example.steady_sensitivity_solve(example.get_parameter("thy"));
    }
    else if (case_name == "transient_conduction_1d") {

        MAST::Examples::ThermalExample1D example(init.comm());
        example.init(*input, prefix);
        example.transient_solve();
        example.transient_sensitivity_solve(example.get_parameter("thy"));
    }
    else if (case_name == "steady_state_conduction_2d") {
        
        MAST::Examples::ThermalExample2D example(init.comm());
        example.init(*input, prefix);
        example.steady_solve();
        example.steady_sensitivity_solve(example.get_parameter("th"));
    }
    else if (case_name == "transient_conduction_2d") {
        
        MAST::Examples::ThermalExample2D example(init.comm());
        example.init(*input, prefix);
        example.transient_solve();
        example.transient_sensitivity_solve(example.get_parameter("th"));
    }
    else {
        libMesh::out
        << "Please run the driver with the name of example specified as: \n"
        << "   run_case=<name> (required, otherwise this message is printed)\n"
        << "   input=<file_name> (optional input file name, otherwise default values of parameters are used) \n"
        << "   prefix=<name> (optional string that will prepend each parameter name) \n"
        << "   print_params=<true/false> (optional flag, if true, the list of problem parameters and their default values will be printed) \n"
        << "Possible values for run_case are:\n\n\n"
        << "**********************************\n"
        << "*********   STRUCTURAL   *********\n"
        << "**********************************\n"
        << "  bar_extension \n"
        << "  beam_bending \n"
        << "  beam_oscillating_load \n"
        << "  beam_modal_analysis\n"
//        << "  beam_thermally_stressed_modal_analysis\n"
//        << "  beam_prestress_buckling_analysis\n"
//        << "  beam_bending_with_offset \n"
//        << "  beam_bending_thermal_stress \n"
//        << "  beam_bending_optimization \n"
//        << "  beam_bending_single_functional_optimization \n"
//        << "  beam_bending_section_offset_optimization \n"
//        << "  beam_bending_thermal_stress_optimization \n"
//        << "  beam_piston_theory_flutter_analysis\n"
//        << "  beam_piston_theory_time_accurate_analysis\n"
//        << "  membrane_extension_uniaxial \n"
//        << "  membrane_extension_biaxial \n"
        << "  plate_bending \n"
//        << "  plate_bending_level_set \n"
//        << "  plate_oscillating_load \n"
//        << "  plate_bending_section_offset \n"
//        << "  plate_bending_thermal_stress \n"
//        << "  stiffened_plate_bending_thermal_stress \n"
//        << "  plate_bending_sizing_optimization \n"
//        << "  plate_bending_single_functional_sizing_optimization \n"
//        << "  plate_bending_section_offset_optimization \n"
//        << "  plate_bending_thermal_stress_optimization \n"
//        << "  plate_modal_analysis\n"
//        << "  plate_thermally_stressed_modal_analysis\n"
//        << "  plate_piston_theory_flutter_analysis\n"
//        << "  plate_thermally_stressed_piston_theory_flutter_analysis\n"
//        << "  plate_prestress_buckling_analysis\n"
//        << "  nastran_model_analysis\n"
//        << "  stiffened_plate_bending_thermal_stress_optimization \n"
//        << "  stiffened_plate_piston_theory_optimization \n"
//        << "  stiffened_plate_thermally_stressed_piston_theory_optimization \n"
        << "  topology_optimization_2D \n"
//        << "*  The default for with_sensitivity is: false.\n"
//        << "*  param is used to specify the parameter name for which sensitivity is desired.\n"
//        << "*  nonlinear is used to turn on/off nonlinear stiffening in the problem.\n"
//        << "*  verify_grads=true will verify the gradients of the optimization problem before calling the optimizer.\n"
        << "\n\n\n"
//        << "**********************************\n"
//        << "***********   FLUID   ************\n"
//        << "**********************************\n"
//        << "  panel_inviscid_analysis_2d \n"
//        << "  ramp_laminar_analysis_2d \n"
//        << "  panel_inviscid_analysis_3d_half_domain \n"
//        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_2d\n"
//        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_3d\n"
//        << "  panel_inviscid_small_disturbance_frequency_domain_analysis_3d_half_domain\n"
//        << "\n\n\n"
        << "**********************************\n"
        << "***********   CONDUCTION   *******\n"
        << "**********************************\n"
        << "  steady_state_conduction_1d \n"
        << "  transient_conduction_1d \n"
        << "  steady_state_conduction_2d \n"
        << "  transient_conduction_2d \n"
        << "\n\n\n"
//        << "**********************************\n"
//        << "***********   FSI     ************\n"
//        << "**********************************\n"
//        << "  beam_fsi_analysis \n"
//        << "  beam_fsi_flutter_analysis \n"
//        << "  beam_flag_fsi_flutter_analysis \n"
//        << "  beam_fsi_ho_convergence \n"
//        << "  beam_fsi_nonuniform_aero_base_flutter_analysis\n"
//        << "  beam_aerothermoelastic_flutter_analysis \n"
//        << "  beam_fsi_flutter_optimization \n"
//        << "  plate_fsi_flutter_analysis \n"
//        << "  plate_flag_fsi_flutter_analysis \n"
//        << "  plate_fsi_half_domain_flutter_analysis \n"
//        << "  plate_fsi_flutter_optimization \n"
//        << "  stiffened_plate_fsi_thermally_stressed_flutter_optimization \n"
//        << "\n\n\n"
        << std::endl;
    }
    
    return 0;
}
