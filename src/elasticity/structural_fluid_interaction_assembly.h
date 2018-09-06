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

#ifndef __mast__structural_fluid_interaction_assembly_h__
#define __mast__structural_fluid_interaction_assembly_h__


// MAST includes
#include "base/assembly_base.h"

namespace MAST {
    
    // Forward declerations
    class RealOutputFunction;
    class FunctionBase;

    enum StructuralQuantityType {
        MASS,                    // mass matrix
        DAMPING,                 // velocity proportional term
        STIFFNESS,               // tangent stiffess matrix
        FORCE,                   // force vector
        INVALID_QTY
    };
    
    
    class StructuralFluidInteractionAssembly:
    public MAST::AssemblyBase {
        
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralFluidInteractionAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralFluidInteractionAssembly();
        
        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
        /*!
         *   if the eigenproblem is defined about a non-zero base solution,
         *   then this method provides the object with the base solution.
         *   The flag \par if_sens tells the method if \par sol
         *   is the sensitivity of the base solution for the current parameter
         *   being solved for
         */
        void set_base_solution(const libMesh::NumericVector<Real>& sol,
                               bool if_sens = false);

        
        /*!
         *   Clears the pointer to base solution.
         *   The flag \par if_sens tells the method to clear the pointer to the
         *   sensitivity vector instead.
         */
        void clear_base_solution(bool if_sens = false);

        /*!
         *   @returns true if a nonzero base solution is used to linearize the
         *   Eigen problem, false otherwise
         */
        bool if_linearized_about_nonzero_solution() const;
        
        
        /*!
         *   @returns a const reference to the base solution (or
         *   its sensitivity when \par if_sens is true) about which
         *   the Eigen problem was linearized.
         */
        const libMesh::NumericVector<Real>&
        base_sol(bool if_sens = false) const;
        
        
        /*!
         *   calculates the reduced order matrix given the basis provided in
         *   \par basis. \par X is the steady state solution about which
         *   the quantity is calculated.
         */
        virtual void
        assemble_reduced_order_quantity
        (std::vector<libMesh::NumericVector<Real>*>& basis,
         std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map);

        
        
        /*!
         *   calculates the sensitivity of reduced order matrix given the basis 
         *   provided in \par basis. \par X is the steady state solution about which
         *   the quantity is calculated, and \par dX_dp is the sensitivity of 
         *   \par X wrt the parameter identified as \par parameters[i]
         */
        virtual void
        assemble_reduced_order_quantity_sensitivity
        (const MAST::FunctionBase& f,
         std::vector<libMesh::NumericVector<Real>*>& basis,
         std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map);

        

    protected:
        
        
        /*!
         *   base solution about which this eigenproblem is defined. This
         *   vector stores the localized values necessary to perform element
         *   calculations.
         */
        const libMesh::NumericVector<Real> * _base_sol;
        
        /*!
         *   sensitivity of base solution may be needed for sensitivity
         *   analysis. This vector stores the localized values necessary to
         *   perform element calculations.
         */
        const libMesh::NumericVector<Real> * _base_sol_sensitivity;
    };
}



#endif // __mast__structural_fluid_interaction_assembly_h__
