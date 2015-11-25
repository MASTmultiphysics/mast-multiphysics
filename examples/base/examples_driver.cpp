/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "examples/structural/bar_extension/bar_extension.h"
#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
#include "examples/structural/beam_bending/beam_bending.h"
#include "examples/structural/beam_bending_with_offset/beam_bending_with_offset.h"
#include "examples/structural/beam_bending_thermal_stress_with_offset/beam_bending_thermal_stress.h"
#include "examples/structural/beam_optimization/beam_optimization.h"
#include "examples/structural/beam_optimization_single_stress_functional/beam_optimization.h"
#include "examples/structural/beam_optimization_section_offset/beam_optimization_section_offset.h"
#include "examples/structural/beam_optimization_thermal_stress/beam_optimization_thermal_stress.h"
#include "examples/structural/membrane_extension_uniaxial_stress/membrane_extension_uniaxial.h"
#include "examples/structural/membrane_extension_biaxial_stress/membrane_extension_biaxial.h"
#include "examples/structural/plate_bending/plate_bending.h"
#include "examples/structural/plate_bending_section_offset/plate_bending_section_offset.h"
#include "examples/structural/plate_bending_thermal_stress/plate_bending_thermal_stress.h"
#include "examples/structural/plate_optimization/plate_optimization.h"
#include "examples/structural/plate_optimization_single_stress_functional/plate_optimization_single_functional.h"
#include "examples/structural/plate_optimization_section_offset/plate_section_offset_optimization.h"
#include "examples/structural/plate_optimization_thermal_stress/plate_thermal_stress_optimization.h"
#include "optimization/npsol_optimization_interface.h"
#include "optimization/dot_optimization_interface.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"


libMesh::LibMeshInit     *_init          = NULL;
MAST::FunctionEvaluation *__my_func_eval = NULL;



template <typename ValType>
void analysis(const std::string& case_name,
              bool with_sens,
              std::string& par_name)  {
    
    ValType run_case;
    
    std::cout << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            std::cout
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, true);
        }
    }

}




template <typename ValType>
void plate_analysis(const std::string& case_name,
                    bool with_sens,
                    std::string& par_name)  {
    
    ValType run_case;
    run_case.init(libMesh::QUAD4, true);
    
    std::cout << "Running case: " << case_name << std::endl;
    run_case.solve(true);
    if (with_sens) {
        MAST::Parameter* p = run_case.get_parameter(par_name);
        if (p) {
            
            std::cout
            << "Running sensitivity for case: " << case_name
            << "  wrt  " << par_name << std::endl;
            run_case.sensitivity_solve(*p, true);
        }
    }
    
}



template <typename ValType>
void optimization(const std::string& case_name)  {

    
    std::cout
    << case_name << std::endl
    << "  input.in should be provided in the working directory with"
    << " desired parameter values."
    << "  In absence of a parameter value, its default value will be used."
    << std::endl
    << "  Output per iteration is written to optimization_output.txt."
    << std::endl;
    
    GetPot infile("input.in");
    std::ofstream output;
    output.open("optimization_output.txt", std::ofstream::out);
    
    MAST::GCMMAOptimizationInterface optimizer;
    
    // create and attach sizing optimization object
    ValType func_eval(infile, output);
    __my_func_eval = &func_eval;
    
    std::vector<Real> dvals(func_eval.n_vars());
    std::fill(dvals.begin(), dvals.end(), 0.05);
    std::cout << "******* Begin: Verifying gradients ***********" << std::endl;
    //func_eval.verify_gradients(dvals);
    std::cout << "******* End: Verifying gradients ***********" << std::endl;
    
    // attach and optimize
    optimizer.attach_function_evaluation_object(func_eval);
    optimizer.optimize();
    
    output.close();
}




int main(int argc, char* const argv[]) {

    libMesh::LibMeshInit init(argc, argv);
    _init  = &init;
    
    // use to get arguments from the command line
    GetPot command_line(argc, argv);
    
    // look for the name that the user has requested to run
    std::string
    case_name = command_line("run_case", ""),
    par_name  = command_line(   "param", "");
    bool
    with_sens = command_line("with_sensitivity", false);
    

    
    if (case_name == "bar_extension")
        analysis<MAST::BarExtension>(case_name, with_sens, par_name);
    else if (case_name == "beam_modal_analysis") {
        libmesh_error_msg("this case has a memory issue that we need to debug");
        //analysis<MAST::BeamModalAnalysis>(case_name, with_sens, par_name);
    }
    else if (case_name == "beam_bending")
        analysis<MAST::BeamBending>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_with_offset")
        analysis<MAST::BeamBendingWithOffset>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_thermal_stress")
        analysis<MAST::BeamBendingThermalStress>(case_name, with_sens, par_name);
    else if (case_name == "beam_bending_optimization")
        optimization<MAST::BeamBendingSizingOptimization>(case_name);
    else if (case_name == "beam_bending_single_functional_optimization")
        optimization<MAST::BeamBendingSingleFunctionalSizingOptimization>(case_name);
    else if (case_name == "beam_bending_section_offset_optimization")
        optimization<MAST::BeamBendingSectionOffsetSizingOptimization>(case_name);
    else if (case_name == "beam_bending_thermal_stress_optimization")
        optimization<MAST::BeamBendingThermalStressSizingOptimization>(case_name);
    else if (case_name == "membrane_extension_uniaxial")
        analysis<MAST::MembraneExtensionUniaxial>(case_name, with_sens, par_name);
    else if (case_name == "membrane_extension_biaxial")
        analysis<MAST::MembraneExtensionBiaxial>(case_name, with_sens, par_name);
    else if (case_name == "plate_bending")
        plate_analysis<MAST::PlateBending>(case_name, with_sens, par_name);
    else if (case_name == "plate_bending_section_offset")
        plate_analysis<MAST::PlateBendingWithOffset>(case_name, with_sens, par_name);
    else if (case_name == "plate_bending_thermal_stress")
        plate_analysis<MAST::PlateBendingThermalStress>(case_name, with_sens, par_name);
    else if (case_name == "plate_bending_sizing_optimization")
        optimization<MAST::PlateBendingSizingOptimization>(case_name);
    else if (case_name == "plate_bending_single_functional_sizing_optimization")
        optimization<MAST::PlateBendingSingleStressFunctionalSizingOptimization>(case_name);
    else if (case_name == "plate_bending_section_offset_optimization")
        optimization<MAST::PlateBendingSectionOffsetSizingOptimization>(case_name);
    else if (case_name == "plate_bending_thermal_stress_optimization")
        optimization<MAST::PlateBendingThermalStressSizingOptimization>(case_name);
    else {
        std::cout
        << "Please run the driver with the name of example specified as: \n"
        << "   run_case=<name>"
        << "   with_sensitivity=<true/false>"
        << "   param=<name>\n\n"
        << "Possible values are:\n"
        << "  bar_extension \n"
        << "  beam_bending \n"
        << "  beam_bending_with_offset \n"
        << "  beam_bending_thermal_stress \n"
        << "  beam_bending_optimization \n"
        << "  beam_bending_single_functional_optimization \n"
        << "  beam_bending_section_offset_optimization \n"
        << "  beam_bending_thermal_stress_optimization \n"
        << "  membrane_extension_uniaxial \n"
        << "  membrane_extension_biaxial \n"
        << "  plate_bending \n"
        << "  plate_bending_section_offset \n"
        << "  plate_bending_thermal_stress \n"
        << "  plate_bending_optimization \n"
        << "  plate_bending_single_functional_sizing_optimization \n"
        << "  plate_bending_section_offset_optimization \n"
        << "  plate_bending_thermal_stress_optimization \n"
        << "*  The default for --with_sensitivity is: false.\n"
        << "*  param is used to specify the parameter name for which sensitivity is desired."
        << std::endl;
    }
    
    return 0;
}
