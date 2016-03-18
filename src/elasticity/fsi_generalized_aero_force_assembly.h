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
         *    initializes for the given fluid and structural components. The
         *    structural and fluid communicator objects should provide valid 
         *    MPI communicators on ranks that store the respective 
         *    disciplinary data structures. \par complex_solver and
         *    \par pressure_function should be non-null pointers only on nodes
         *    with a valid fluid communicator. \par motion_func should be
         *   a non-null pointer only when structural_comm is a valid 
         *   communicator.
         */
        void init(MAST::FrequencyFunction&                freq,
                  libMesh::Parallel::Communicator&        structural_comm,
                  libMesh::Parallel::Communicator&        fluid_comm,
                  MAST::ComplexSolverBase*                complex_solver,
                  MAST::SmallDisturbancePressureFunction* pressure_func,
                  MAST::FlexibleSurfaceMotion*            motion_func);
        
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
        /*!
         *   calculates the reduced order matrix given the basis provided in
         *   \par basis. \par X is the steady state solution about which
         *   the quantity is calculated.
         */
        virtual void
        assemble_generalized_aerodynamic_force_matrix
        (std::vector<libMesh::NumericVector<Real>*>& basis,
         ComplexMatrixX& mat);
        

        virtual void
        assemble_reduced_order_quantity
        (std::vector<libMesh::NumericVector<Real>*>& basis,
         std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map);
        
        
    protected:
        
        
        /*!
         *   frequency function
         */
        MAST::FrequencyFunction                      *_freq;

        
        /*!
         *   communicator for structural model
         */
        libMesh::Parallel::Communicator             *_structural_comm;

        
        /*!
         *   communicator for fluid model
         */
        libMesh::Parallel::Communicator             *_fluid_comm;

        
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

