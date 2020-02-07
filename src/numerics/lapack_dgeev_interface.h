/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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

#ifndef __mast__lapack_dgeev_interface_h__
#define __mast__lapack_dgeev_interface_h__


// MAST includes
#include "base/mast_data_types.h"


extern "C" {
    
    /*
     *  =====================================================================
     *  Purpose
     *  =======
     *
     *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
     *  eigenvalues and, optionally, the left and/or right eigenvectors.
     *
     *  The right eigenvector v(j) of A satisfies
     *                   A * v(j) = lambda(j) * v(j)
     *  where lambda(j) is its eigenvalue.
     *  The left eigenvector u(j) of A satisfies
     *                u(j)**H * A = lambda(j) * u(j)**H
     *  where u(j)**H denotes the conjugate transpose of u(j).
     *
     *  The computed eigenvectors are normalized to have Euclidean norm
     *  equal to 1 and largest component real.
     *
     *  Arguments
     *  =========
     *
     *  JOBVL   (input) CHARACTER*1
     *          = 'N': left eigenvectors of A are not computed;
     *          = 'V': left eigenvectors of A are computed.
     *
     *  JOBVR   (input) CHARACTER*1
     *          = 'N': right eigenvectors of A are not computed;
     *          = 'V': right eigenvectors of A are computed.
     *
     *  N       (input) INTEGER
     *          The order of the matrix A. N >= 0.
     *
     *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
     *          On entry, the N-by-N matrix A.
     *          On exit, A has been overwritten.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of the array A.  LDA >= max(1,N).
     *
     *  WR      (output) DOUBLE PRECISION array, dimension (N)
     *  WI      (output) DOUBLE PRECISION array, dimension (N)
     *          WR and WI contain the real and imaginary parts,
     *          respectively, of the computed eigenvalues.  Complex
     *          conjugate pairs of eigenvalues appear consecutively
     *          with the eigenvalue having the positive imaginary part
     *          first.
     *
     *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
     *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
     *          after another in the columns of VL, in the same order
     *          as their eigenvalues.
     *          If JOBVL = 'N', VL is not referenced.
     *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
     *          the j-th column of VL.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
     *          u(j+1) = VL(:,j) - i*VL(:,j+1).
     *
     *  LDVL    (input) INTEGER
     *          The leading dimension of the array VL.  LDVL >= 1; if
     *          JOBVL = 'V', LDVL >= N.
     *
     *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
     *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
     *          after another in the columns of VR, in the same order
     *          as their eigenvalues.
     *          If JOBVR = 'N', VR is not referenced.
     *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
     *          the j-th column of VR.
     *          If the j-th and (j+1)-st eigenvalues form a complex
     *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
     *          v(j+1) = VR(:,j) - i*VR(:,j+1).
     *
     *  LDVR    (input) INTEGER
     *          The leading dimension of the array VR.  LDVR >= 1; if
     *          JOBVR = 'V', LDVR >= N.
     *
     *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
     *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     *
     *  LWORK   (input) INTEGER
     *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
     *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
     *          performance, LWORK must generally be larger.
     *
     *          If LWORK = -1, then a workspace query is assumed; the routine
     *          only calculates the optimal size of the WORK array, returns
     *          this value as the first entry of the WORK array, and no error
     *          message related to LWORK is issued by XERBLA.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value.
     *          > 0:  if INFO = i, the QR algorithm failed to compute all the
     *                eigenvalues, and no eigenvectors have been computed;
     *                elements i+1:N of WR and WI contain eigenvalues which
     *                have converged.
     *  =====================================================================
     */
    extern void dgeev_(const char*    jobvl,
                       const char*    jobvr,
                       int*     n,
                       double*  a,
                       int*     lda,
                       double*  w_r,
                       double*  w_i,
                       double*  vl,
                       int*     ldvl,
                       double*  vr,
                       int*     ldvr,
                       double*  work,
                       int*     lwork,
                       int*     info);
    
}


namespace MAST {
    
    
    class LAPACK_DGEEV{
        
    public:
        
        LAPACK_DGEEV():
        info_val(-1)
        { }
        
        /*!
         *    computes the eigensolution for \f$ A x = \lambda I x\f$. A & B will be
         *    overwritten
         */
        void compute(const RealMatrixX& A,
                     bool computeEigenvectors = true);
        
        ComputationInfo info() const;
        
        const RealMatrixX& A() const {
            libmesh_assert(info_val == 0);
            return this->_A;
        }
        
        
        const ComplexVectorX& eig_vals() const {
            libmesh_assert(info_val == 0);
            return this->W;
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
            ComplexMatrixX r = this->VL.conjugate().transpose() * this->VR;
            
            // scale the right eigenvectors by the inverse of the inner-product
            // diagonal
            Complex val;
            for (unsigned int i=0; i<_A.cols(); i++) {
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
            
            r = this->VL.conjugate().transpose() * this->VR;
            out << "conj(VL)' * B * VR" << std::endl
            << r << std::endl;
            
        }
        
    protected:
        
        RealMatrixX    _A;
        
        ComplexMatrixX VL;
        
        ComplexMatrixX VR;
        
        ComplexVectorX W;
        
        int info_val;
    };
    
}

#endif // __mast__lapack_dgeev_interface_h__
