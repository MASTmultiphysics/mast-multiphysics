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

#ifndef __mast__level_set_transient_assembly__
#define __mast__level_set_transient_assembly__

// MAST includes
#include "base/transient_assembly_elem_operations.h"



namespace MAST {
    
    
    class LevelSetTransientAssemblyElemOperations:
    public MAST::TransientAssemblyElemOperations {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        LevelSetTransientAssemblyElemOperations();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~LevelSetTransientAssemblyElemOperations();
        
        
        /*!
         *   set element reference solution for reinitialization
         */
        void set_elem_reference_solution(const RealVectorX& sol);
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_x_dot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac);
        
        /*!
         *   This call for second order ode should not be used for this
         *   transient assembly
         */
        virtual void elem_calculations(bool if_jac,
                                       RealVectorX& f_m,
                                       RealVectorX& f_x,
                                       RealMatrixX& f_m_jac_xddot,
                                       RealMatrixX& f_m_jac_xdot,
                                       RealMatrixX& f_m_jac,
                                       RealMatrixX& f_x_jac_xdot,
                                       RealMatrixX& f_x_jac) {
            libmesh_error();
        }
        
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
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                                   RealVectorX& vec);
        
        
        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        elem_second_derivative_dot_solution_assembly(RealMatrixX& mat);
        
        virtual void
        elem_sensitivity_calculations(const MAST::FunctionBase& f,
                                      RealVectorX& f_m,
                                      RealVectorX& f_x);

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
        
    };
    
    
}

#endif // __mast__level_set_transient_assembly__
