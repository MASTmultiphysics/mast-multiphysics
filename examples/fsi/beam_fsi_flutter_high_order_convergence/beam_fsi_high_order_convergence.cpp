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
//#include <fstream>
//#include <iomanip>
//#include <ctime>
//
//
//// MAST includes
//#include "examples/fsi/beam_fsi_flutter_high_order_convergence/beam_fsi_high_order_convergence.h"
//#include "examples/fsi/beam_flutter_solution/beam_euler_fsi_flutter_solution.h"
//#include "base/nonlinear_system.h"
//#include "aeroelasticity/flutter_root_base.h"
//
//
//
//
//
//MAST::BeamFSIFlutterHighOrderConvergence::
//BeamFSIFlutterHighOrderConvergence(const libMesh::Parallel::Communicator& comm_in):
//libMesh::ParallelObject (comm_in) {
//    
//}
//
//
//
//
//
//
//MAST::BeamFSIFlutterHighOrderConvergence::~BeamFSIFlutterHighOrderConvergence() {
//    
//}
//
//
//
//
//Real
//MAST::BeamFSIFlutterHighOrderConvergence::solve(bool if_write_output,
//                                                const Real tol,
//                                                const unsigned int max_bisection_iters) {
//
//    const unsigned int
//    n_cpu = this->comm().size(),
//    rank  = this->comm().rank();
//    
//    std::ofstream output;
//    
//    if (rank == 0) {
//        output.open("convergence_output.txt", std::ofstream::out);
//        
//        output
//        << "N Procs : " << n_cpu << std::endl << std::endl
//        << std::setw(5) <<  "O"
//        << std::setw(5) <<  "Mesh"
//        << std::setw(15) << "N_Dofs"
//        << std::setw(35) << "g"
//        << std::setw(35) << "omega(rad/s)"
//        << std::setw(35) << "Velocity"
//        << std::setw(35) << "Wall_Time(sec)"
//        << std::setw(35) << "CPU_Time(sec)" << std::endl;
//    }
//    
//    for (unsigned int o=0; o<1; o++) { // order loop
//    
//        unsigned int low = 0;
//        if (o == 0) low = 2;
//        for (unsigned int r=4; r<5; r++) {  // h-refinement loop
//
//            libMesh::out
//            << "**************************************************************************" << std::endl
//            << "             Analysis for p = " << o+1 << "  mesh = " << r+1 << std::endl
//            << "**************************************************************************" << std::endl;
//            
//
//            std::clock_t start;
//            Real duration;
//            
//            start = std::clock();
//
//            // create the analysis object with the specified level of interpolation
//            MAST::BeamEulerFSIFlutterAnalysis solver(this->comm(), o, r);
//            unsigned int n_fluid_dofs = solver._fluid_sys->n_dofs();
//            
//            solver.solve(false , 1.e-8, 100);
//
//            duration = ( std::clock() - start ) / (Real) CLOCKS_PER_SEC;
//
//            if (rank == 0 && solver._flutter_root)
//                output
//                << std::setw(5) << (o+1)
//                << std::setw(5) << (r)
//                << std::setw(15) << n_fluid_dofs
//                << std::setw(35) << std::setprecision(15) << solver._flutter_root->g
//                << std::setw(35) << std::setprecision(15) << solver._flutter_root->omega
//                << std::setw(35) << std::setprecision(15) << solver._flutter_root->V
//                << std::setw(35) << std::setprecision(15) << duration
//                << std::setw(35) << std::setprecision(15) << duration*n_cpu << std::endl;
//            
//        }
//        
//        if (rank == 0)
//            output << std::endl << std::endl;
//
//    }
//
//    return 0.;
//}
//
//
//
//
//
