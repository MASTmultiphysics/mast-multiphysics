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
#include "base/assembly_elem_operation.h"

namespace MAST {
    
    // Forward declerations
    class RealOutputFunction;
    class FunctionBase;

    enum StructuralQuantityType {
        MASS,                    // mass matrix
        DAMPING,                 // velocity proportional term
        STIFFNESS,               // tangent stiffess matrix
        FORCE                    // force vector
    };
    
    
    class StructuralFluidInteractionAssembly:
    public MAST::AssemblyBase,
    public MAST::AssemblyElemOperations {
        
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
         *   Reattaches to the same system that was attached earlier.
         *
         *   This cannot be called if the clear_discipline_and_system() method
         *   has been called.
         */
        inline virtual void
        reattach_to_system() {
            // nothing to be done here.
        }

        /*!
         *   clears association with a system to this discipline, and vice-a-versa
         */
        virtual void
        clear_discipline_and_system( );

        
        /*!
         *   tells the object about the quantity to be assembled in the 
         *   matrix
         */
        inline void
        set_quantity_to_assemble(MAST::StructuralQuantityType qty) {
            
            _qty_type = qty;
        }
        
                
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

        
        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const libMesh::Elem& elem);

        /*!
         *   some simulations frequently deal with 1D/2D elements in 3D space,
         *   which requires use of MAST::LocalElemFE.
         */
        virtual bool
        if_use_local_elem() const {
            
            return true;
        }
        
        /*!
         *   sets additional data for local elem FE.
         */
        virtual void
        set_local_fe_data(MAST::LocalElemFE& fe,
                          const libMesh::Elem& e) const;

    protected:
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(bool if_jac,
                                        RealVectorX& vec,
                                        RealMatrixX& mat);

        
        /*!
         *   performs the frequency domain element aerodynamic force
         *   calculations over \par elem, and returns
         *   the element vector quantity in
         *   \par vec, respectively.
         */
        virtual void _elem_aerodynamic_force_calculations(ComplexVectorX& vec);

        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(bool if_jac,
                                                    RealVectorX& vec,
                                                    RealMatrixX& mat);

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        _elem_second_derivative_dot_solution_assembly(RealMatrixX& mat);

        /*!
         *   this defines the quantity to be assembled
         */
        MAST::StructuralQuantityType _qty_type;
        
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
