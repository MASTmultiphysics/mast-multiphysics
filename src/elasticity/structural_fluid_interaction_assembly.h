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

#ifndef __mast__structural_fluid_interaction_assembly_h__
#define __mast__structural_fluid_interaction_assembly_h__


// MAST includes
#include "base/nonlinear_implicit_assembly.h"

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
    public MAST::NonlinearImplicitAssembly {
        
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
         *   tells the object about the quantity to be assembled in the 
         *   matrix
         */
        inline void
        set_quantity_to_assemble(MAST::StructuralQuantityType qty) {
            
            _qty_type = qty;
        }
        
        
        /*!
         *    function that assembles the matrices and vectors quantities for
         *    nonlinear solution
         */
        virtual void
        residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>* R,
                               libMesh::SparseMatrix<Real>*  J,
                               libMesh::NonlinearImplicitSystem& S);
        
        
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
        (const libMesh::ParameterVector& parameters,
         const unsigned int i,
         std::vector<libMesh::NumericVector<Real>*>& basis,
         std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map);

        
    protected:
        
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::auto_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& vec,
                                        RealMatrixX& mat);

        
        /*!
         *   performs the frequency domain element aerodynamic force
         *   calculations over \par elem, and returns
         *   the element vector quantity in
         *   \par vec, respectively.
         */
        virtual void _elem_aerodynamic_force_calculations(MAST::ElementBase& elem,
                                                          ComplexVectorX& vec);

        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    bool if_jac,
                                                    RealVectorX& vec,
                                                    RealMatrixX& mat);

        
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
