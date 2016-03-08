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

#ifndef __mast__fsi_generalized_aerodynamic_force_matrix_driver_h__
#define __mast__fsi_generalized_aerodynamic_force_matrix_driver_h__


// MAST includes
#include "elasticity/structural_fluid_interaction_assembly.h"



namespace MAST {
    
    // forward declerations
    class ComplexSolverBase;
    class SmallDisturbancePressureFunction;
    class FlexibleSurfaceMotion;
    class StructuralFluidInteractionAssembly;
    class FrequencyFunction;
    
    
    class FSIGeneralizedAeroForceAssembly:
    public MAST::StructuralFluidInteractionAssembly {
    
    public:
        
        /*!
         *   default constructor
         */
        FSIGeneralizedAeroForceAssembly();
        
        
        /*!
         *   destructor
         */
        ~FSIGeneralizedAeroForceAssembly();
        
        
        /*!
         *    initializes for the given fluid and structural components
         */
        void init(MAST::FrequencyFunction& freq,
                  MAST::ComplexSolverBase& complex_solver,
                  MAST::SmallDisturbancePressureFunction& pressure_func,
                  MAST::FlexibleSurfaceMotion& motion_func);
        
        
        /*!
         *   calculates the reduced order matrix given the basis provided in
         *   \par basis. \par X is the steady state solution about which
         *   the quantity is calculated.
         */
        virtual void
        assemble_generalized_aerodynamic_force_matrix
        (const libMesh::NumericVector<Real>& X,
         std::vector<libMesh::NumericVector<Real>*>& basis,
         ComplexMatrixX& mat);
        
        
    protected:
        
        
        /*!
         *   frequency function
         */
        MAST::FrequencyFunction                      *_freq;

        
        /*!
         *   complex solver
         */
        MAST::ComplexSolverBase                      *_fluid_complex_solver;
        

        
        /*!
         *  small disturbance pressure function boundary condition for
         *  structures
         */
        MAST::SmallDisturbancePressureFunction       *_pressure_function;
        
        
        /*!
         *   flexible surface motion for fluid and structure
         */
        MAST::FlexibleSurfaceMotion                  *_motion;
    };
}


#endif // __mast__fsi_generalized_aerodynamic_force_matrix_driver_h__

