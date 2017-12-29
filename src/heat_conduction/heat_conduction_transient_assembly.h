/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__heat_conduction_transient_assembly__
#define __mast__heat_conduction_transient_assembly__

// MAST includes
#include "base/transient_assembly.h"



namespace MAST {
    
    
    class HeatConductionTransientAssembly:
    public MAST::TransientAssembly {
    public:
        
        /*!
         *   constructor associates this assembly object with the system
         */
        HeatConductionTransientAssembly();
        
        
        /*!
         *   destructor resets the association of this assembly object with
         *   the system
         */
        virtual ~HeatConductionTransientAssembly();
        
        //**************************************************************
        //these methods are provided for use by the solvers
        //**************************************************************
        
        
        /*!
         *   performs the element calculations over \par elem, and returns
         *   the element vector and matrix quantities in \par mat and
         *   \par vec, respectively. \par if_jac tells the method to also
         *   assemble the Jacobian, in addition to the residual vector.
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
                                        RealVectorX& f_m,
                                        RealVectorX& f_x,
                                        RealMatrixX& f_m_jac_x_dot,
                                        RealMatrixX& f_m_jac,
                                        RealMatrixX& f_x_jac);
        
        /*!
         *   This call for second order ode should not be used for this 
         *   transient assembly
         */
        virtual void _elem_calculations(MAST::ElementBase& elem,
                                        bool if_jac,
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
        _linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                              RealVectorX& f);

        /*!
         *   performs the element sensitivity calculations over \par elem,
         *   and returns the element residual sensitivity in \par vec .
         */
        virtual void _elem_sensitivity_calculations(MAST::ElementBase& elem,
                                                    RealVectorX& vec);
        

        /*!
         *   calculates \f$ d ([J] \{\Delta X\})/ dX  \f$ over \par elem,
         *   and returns the matrix in \par vec .
         */
        virtual void
        _elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                                      RealMatrixX& mat);

    protected:
        
        
        /*!
         *   @returns a smart-pointer to a newly created element for
         *   calculation of element quantities.
         */
        virtual std::unique_ptr<MAST::ElementBase>
        _build_elem(const libMesh::Elem& elem);

    };
    
    
}

#endif // __mast__heat_conduction_transient_assembly__
