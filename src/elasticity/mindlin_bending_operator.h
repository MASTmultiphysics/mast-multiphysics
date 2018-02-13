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

#ifndef __mast_mindlin_bending_operator_h__
#define __mast_mindlin_bending_operator_h__

// MAST includes
#include "elasticity/bending_operator.h"


namespace MAST {
    
    
    class MindlinBendingOperator:
    public MAST::BendingOperator2D {
        
    public:
        
        MindlinBendingOperator(MAST::StructuralElementBase& elem);
        
        virtual ~MindlinBendingOperator();
        
        /*!
         *   returns true if this bending operator supports a transverse shear component
         */
        virtual bool include_transverse_shear_energy() const {
            return true;
        }
        
        /*!
         *   initialze the bending strain operator for Mindlin element, withouth
         *   the z-location. This is useful for use with element stiffness matrix
         *   integration where the D matrix is calculated by section integration by
         *   the ElementPropertyCard2D.
         */
        virtual void
        initialize_bending_strain_operator (const MAST::FEBase& fe,
                                            const unsigned int qp,
                                            MAST::FEMOperatorMatrix& Bmat);
        
        /*!
         *    initializes the bending strain operator for the specified quadrature
         * point and z-location.
         */
        void
        initialize_bending_strain_operator_for_z(const MAST::FEBase& fe,
                                                 const unsigned int qp,
                                                 const Real z,
                                                 MAST::FEMOperatorMatrix& Bmat_bend);
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual(bool request_jacobian,
                                            RealVectorX& local_f,
                                            RealMatrixX& local_jac);

        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX& local_f,
                                                        RealMatrixX& local_jac);
        
        /*!
         *   calculate the transverse shear component for the element
         */
        virtual void
        calculate_transverse_shear_residual_boundary_velocity
        (const MAST::FunctionBase& p,
         const unsigned int s,
         const MAST::FieldFunction<RealVectorX>& vel_f,
         RealVectorX& local_f);

    protected:
        
        void
        _transverse_shear_operations(const std::vector<std::vector<Real> >& phi,
                                     const std::vector<std::vector<libMesh::RealVectorValue> >& dphi,
                                     const std::vector<Real>& JxW,
                                     const unsigned int     qp,
                                     const RealMatrixX&     material,
                                     FEMOperatorMatrix&     Bmat_trans,
                                     RealVectorX&           phi_vec,
                                     RealVectorX&           vec_n2,
                                     RealVectorX&           vec_2,
                                     RealMatrixX&           mat_n2n2,
                                     RealMatrixX&           mat_2n2,
                                     bool                   request_jacobian,
                                     RealVectorX&           local_f,
                                     RealMatrixX&           local_jac);
        
        /*!
         *   reduction in quadrature for shear energy
         */
        unsigned int _shear_quadrature_reduction;
    };
}



#endif // __mast_mindlin_bending_operator_h__
