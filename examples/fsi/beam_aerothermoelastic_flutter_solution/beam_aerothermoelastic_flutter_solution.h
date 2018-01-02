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

#ifndef __mast_beam_aerothermoelastic_euler_flutter_solution_h__
#define __mast_beam_aerothermoelastic_euler_flutter_solution_h__


// C++ includes
#include <memory>
#include <vector>

// MAST includes
#include "examples/fsi/beam_fsi_solution/beam_euler_fsi_solution.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    // Forward declerations
    class UGFlutterSolver;
    class FlutterRootBase;
    class ComplexMeshFieldFunction;
    class ComplexNormalRotationMeshFunction;

    
    struct BeamAerothermoelasticFlutterSolution:
    public MAST::BeamEulerFSIAnalysis {
    
        
        BeamAerothermoelasticFlutterSolution();
        
        
        ~BeamAerothermoelasticFlutterSolution();
        
        
        
        /*!
         *  solves the system and returns the flutter velocity
         */
        Real solve(bool if_write_output = false,
                   const Real tol = 1.e-1,
                   const unsigned int max_bisection_iters = 20);
        
        
        /*!
         *   parameters for flutter analysis
         */
        Real
        _k_lower,
        _k_upper;
        
        
        unsigned int
        _n_k_divs;

        
        MAST::ComplexMeshFieldFunction              *_displ_perturb;
        MAST::ComplexNormalRotationMeshFunction     *_normal_rot_perturb;

        MAST::FrequencyDomainPressureFunction *_freq_domain_pressure_function;

        
        // parameters used in the system
        MAST::Parameter
        *_omega,
        *_velocity,
        *_b_ref,
        *_zero_omega;
        
        
        MAST::ConstantFieldFunction
        *_omega_f,
        *_velocity_f,
        *_b_ref_f,
        *_zero_omega_f;

        
        /*!
         *   frequency object
         */
        MAST::FrequencyFunction
        *_zero_freq_func,
        *_freq_function;

        /*!
         *   flutter solver to find the bifurcation point
         */
        MAST::UGFlutterSolver*                   _flutter_solver;
        
        /*!
         *   flutter root from the analysis
         */
        MAST::FlutterRootBase*                   _flutter_root;
        
        // vector of basis vectors from modal analysis
        std::vector<libMesh::NumericVector<Real>*> _basis;
    };
}



#endif //  __mast_beam_aerothermoelastic_euler_flutter_solution_h__

