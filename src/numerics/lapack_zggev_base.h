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

#ifndef __mast__lapack_zggev_interface_base_h__
#define __mast__lapack_zggev_interface_base_h__


// MAST includes
#include "base/mast_data_types.h"




namespace MAST {
    
    class LAPACK_ZGGEV_Base{
        
    public:
        
        LAPACK_ZGGEV_Base():
        info_val(-1)
        { }
        
        /*!
         *    computes the eigensolution for A x = \lambda B x. A & B will be
         *    overwritten
         */
        virtual void compute(const ComplexMatrixX& A,
                             const ComplexMatrixX& B,
                             bool computeEigenvectors = true) = 0;
        
        const ComplexMatrixX& A() const {
            libmesh_assert(info_val == 0);
            return this->_A;
        }
        
        
        const ComplexMatrixX& B() const {
            libmesh_assert(info_val == 0);
            return this->_B;
        }
        
        
        const ComplexVectorX& alphas() const {
            libmesh_assert(info_val == 0);
            return this->alpha;
        }
        
        const ComplexVectorX& betas() const {
            libmesh_assert(info_val == 0);
            return this->beta;
        }
        
        const ComplexMatrixX& left_eigenvectors() const {
            libmesh_assert(info_val == 0);
            return this->VL;
        }
        
        const ComplexMatrixX& right_eigenvectors() const {
            libmesh_assert(info_val == 0);
            return this->VR;
        }
        
        /*!
         *    Scales the right eigenvector so that the inner product with respect
         *    to the B matrix is equal to an Identity matrix, i.e.
         *    VL* B * VR = I
         */
        void scale_eigenvectors_to_identity_innerproduct() {
            libmesh_assert(info_val == 0);
            
            // this product should be an identity matrix
            ComplexMatrixX r = this->VL.conjugate().transpose() * _B * this->VR;
            
            // scale the right eigenvectors by the inverse of the inner-product
            // diagonal
            Complex val;
            for (unsigned int i=0; i<_B.cols(); i++) {
                val = r(i,i);
                if (std::abs(val) > 0.)
                    this->VR.col(i) *= (1./val);
            }
        }
        
        void print_inner_product(std::ostream& out) const {
            libmesh_assert(info_val == 0);
            ComplexMatrixX r;
            r = this->VL.conjugate().transpose() * _A * this->VR;
            out << "conj(VL)' * A * VR" << std::endl
            << r << std::endl;
            
            r = this->VL.conjugate().transpose() * _B * this->VR;
            out << "conj(VL)' * B * VR" << std::endl
            << r << std::endl;
            
        }
        
    protected:
        
        ComplexMatrixX _A;
        
        ComplexMatrixX _B;
        
        ComplexMatrixX VL;
        
        ComplexMatrixX VR;
        
        ComplexVectorX alpha;
        
        ComplexVectorX beta;
        
        int info_val;
    };
}



#endif // __mast__lapack_zggev_interface_base_h__
