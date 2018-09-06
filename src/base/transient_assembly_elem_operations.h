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


#ifndef __mast_transient_assembly_elem_operations_h__
#define __mast_transient_assembly_elem_operations_h__

// MAST includes
#include "base/assembly_elem_operation.h"
#include "base/mast_data_types.h"

namespace MAST {
    
    class TransientAssemblyElemOperations:
    public MAST::AssemblyElemOperations {
        
    public:
        TransientAssemblyElemOperations();
        
        virtual ~TransientAssemblyElemOperations();
        
        /*!
         *   performs the element calculations over \par elem, for a
         *   system of the form
         *   \f[ f_m(x,\dot{x}) + f_x(x) = 0 \f].
         *   \param f_x     = \f$ f(x) \f$
         *   \param f_m    =  \f$ f_m(x,\dot{x}) \f$
         *   \param f_m_jac_xdot = \f$ \frac{\partial (f_m(x)}{\partial \dot{x}} \f$
         *   \param f_m_jac = \f$ \frac{\partial (f_m(x)}{\partial x} \f$
         *   \param f_x_jac = \f$ \frac{\partial (f(x)}{\partial x} \f$
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_xdot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac) = 0;
        
        
        
        /*!
         *   performs the element calculations over \par elem, for a
         *   system of the form
         *   \f[ f_m(x,\ddot{x}, \dot{x}) + f_x(x, \dot{x}) = 0 \f].
         *   \param f_x     = \f$ f(x,\dot{x})  \f$
         *   \param f_m    =  \f$ f_m(x,\dot{x}) \f$
         *   \param f_m_jac_xddot = \f$ \frac{\partial (f_m(x)}{\partial \ddot{x}} \f$
         *   \param f_m_jac_xdot = \f$ \frac{\partial (f_m(x)}{\partial \dot{x}} \f$
         *   \param f_m_jac = \f$ \frac{\partial (f_m(x)}{\partial x} \f$
         *   \param f_x_jac_xdot = \f$ \frac{\partial (f(x)}{\partial \dot{x}} \f$
         *   \param f_x_jac = \f$ \frac{\partial (f(x)}{\partial x} \f$
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_xddot,
                                       RealMatrixX& f_m_jac_xdot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac_xdot,
                                       RealMatrixX& f_x_jac) = 0;
        
        /*!
         *   Calculates the product of Jacobian-solution, and Jacobian-velocity
         *   over the element for a system of the form
         *   \f[ f_m(x,\dot{x}) + f_x(x) = 0 \f].
         *   \param f    = \f$ df(x)/dx \cdot dx +
         *    df_m(x,\dot{x})/dx \cdot dx +
         *    df_m(x,\dot{x})/d\dot{x} \cdot d{\dot x} \f$
         */
        virtual void
        linearized_jacobian_solution_product(RealVectorX& f) = 0;
        
        
        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the component of  element residual sensitivity in
         *   \par f_m and \par f_x.
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   RealVectorX& f_m,
                                                   RealVectorX& f_x) = 0;


    protected:
        
    };
}


#endif // __mast_transient_assembly_elem_operations_h__

