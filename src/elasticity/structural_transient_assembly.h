/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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

#ifndef __mast__structural_transient_assembly_elem_operations_h__
#define __mast__structural_transient_assembly_elem_operations_h__


// MAST includes
#include "base/transient_assembly_elem_operations.h"



namespace MAST {
    
    
    class StructuralTransientAssemblyElemOperations:
    public MAST::TransientAssemblyElemOperations {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        StructuralTransientAssemblyElemOperations();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~StructuralTransientAssemblyElemOperations();
        
        //**************************************************************
        //these methods are provided for use by the solvers
        //**************************************************************
        
        
        /*!
         *   This call for first order ode should not be used for this
         *   transient assembly
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_x_dot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac);
        
        /*!
         *   performs the element calculations over \p elem, and returns
         *   the element vector and matrix quantities in \p mat and
         *   \p vec, respectively. \p if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_xddot,
                                       RealMatrixX& f_m_jac_xdot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac_xdot,
                                       RealMatrixX& f_x_jac);
        
        /*!
         *   Calculates the product of Jacobian-solution, and Jacobian-velocity
         *   over the element for a system of the form
         *   \f[ f_m(x,\dot{x}) + f_x(x) = 0 \f].
         *   \param f    = \f$ df(x)/dx \cdot dx +
         *    df_m(x,\dot{x})/dx \cdot dx +
         *    df_m(x,\dot{x})/d\dot{x} \cdot d{\dot x} \f$
         */
        virtual void
        linearized_jacobian_solution_product(RealVectorX& f);

        
        /*!
         *   performs the element sensitivity calculations over \p elem,
         *   and returns the component of  element residual sensitivity in
         *   \p f_m and \p f_x.
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   RealVectorX& f_m,
                                                   RealVectorX& f_x);
        
        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \p elem,
         *   and returns the matrix in \p vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat);

        /*!
         *   initializes the object for the geometric element \p elem. This
         *   expects the object to be in a cleared state, so the user should
         *   call \p clear_elem() between successive initializations.
         */
        virtual void
        init(const libMesh::Elem& elem);

    protected:
        
        
    };
    
    
}


#endif // __mast__structural_transient_assembly_elem_operations_h__

