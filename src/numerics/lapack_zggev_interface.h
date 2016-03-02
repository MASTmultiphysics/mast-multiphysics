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

#ifndef __mast__lapack_zggev_interface_h__
#define __mast__lapack_zggev_interface_h__


// MAST includes
#include "base/mast_data_types.h"


extern "C" {
    /*
     *  =====================================================================
     *  Purpose
     *  =======
     *
     *  ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices
     *  (A,B), the generalized eigenvalues, and optionally, the left and/or
     *  right generalized eigenvectors.
     *
     *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     *  singular. It is usually represented as the pair (alpha,beta), as
     *  there is a reasonable interpretation for beta=0, and even for both
     *  being zero.
     *
     *  The right generalized eigenvector v(j) corresponding to the
     *  generalized eigenvalue lambda(j) of (A,B) satisfies
     *
     *               A * v(j) = lambda(j) * B * v(j).
     *
     *  The left generalized eigenvector u(j) corresponding to the
     *  generalized eigenvalues lambda(j) of (A,B) satisfies
     *
     *               u(j)**H * A = lambda(j) * u(j)**H * B
     *
     *  where u(j)**H is the conjugate-transpose of u(j).
     *
     *  Arguments
     *  =========
     *
     *  JOBVL   (input) CHARACTER*1
     *          = 'N':  do not compute the left generalized eigenvectors;
     *          = 'V':  compute the left generalized eigenvectors.
     *
     *  JOBVR   (input) CHARACTER*1
     *          = 'N':  do not compute the right generalized eigenvectors;
     *          = 'V':  compute the right generalized eigenvectors.
     *
     *  N       (input) INTEGER
     *          The order of the matrices A, B, VL, and VR.  N >= 0.
     *
     *  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
     *          On entry, the matrix A in the pair (A,B).
     *          On exit, A has been overwritten.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of A.  LDA >= max(1,N).
     *
     *  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
     *          On entry, the matrix B in the pair (A,B).
     *          On exit, B has been overwritten.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of B.  LDB >= max(1,N).
     *
     *  ALPHA   (output) COMPLEX*16 array, dimension (N)
     *  BETA    (output) COMPLEX*16 array, dimension (N)
     *          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
     *          generalized eigenvalues.
     *
     *          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
     *          underflow, and BETA(j) may even be zero.  Thus, the user
     *          should avoid naively computing the ratio alpha/beta.
     *          However, ALPHA will be always less than and usually
     *          comparable with norm(A) in magnitude, and BETA always less
     *          than and usually comparable with norm(B).
     *
     *  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
     *          If JOBVL = 'V', the left generalized eigenvectors u(j) are
     *          stored one after another in the columns of VL, in the same
     *          order as their eigenvalues.
     *          Each eigenvector is scaled so the largest component has
     *          abs(real part) + abs(imag. part) = 1.
     *          Not referenced if JOBVL = 'N'.
     *
     *  LDVL    (input) INTEGER
     *          The leading dimension of the matrix VL. LDVL >= 1, and
     *          if JOBVL = 'V', LDVL >= N.
     *
     *  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
     *          If JOBVR = 'V', the right generalized eigenvectors v(j) are
     *          stored one after another in the columns of VR, in the same
     *          order as their eigenvalues.
     *          Each eigenvector is scaled so the largest component has
     *          abs(real part) + abs(imag. part) = 1.
     *          Not referenced if JOBVR = 'N'.
     *
     *  LDVR    (input) INTEGER
     *          The leading dimension of the matrix VR. LDVR >= 1, and
     *          if JOBVR = 'V', LDVR >= N.
     *
     *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
     *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     *
     *  LWORK   (input) INTEGER
     *          The dimension of the array WORK.  LWORK >= max(1,2*N).
     *          For good performance, LWORK must generally be larger.
     *
     *          If LWORK = -1, then a workspace query is assumed; the routine
     *          only calculates the optimal size of the WORK array, returns
     *          this value as the first entry of the WORK array, and no error
     *          message related to LWORK is issued by XERBLA.
     *
     *  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value.
     *          =1,...,N:
     *                The QZ iteration failed.  No eigenvectors have been
     *                calculated, but ALPHA(j) and BETA(j) should be
     *                correct for j=INFO+1,...,N.
     *          > N:  =N+1: other then QZ iteration failed in DHGEQZ,
     *                =N+2: error return from DTGEVC.
     *
     *  =====================================================================
     */
    extern int zggev_(char*                 jobvl,
                      char*                 jobvr,
                      int*                  n,
                      std::complex<double>* a,
                      int*                  lda,
                      std::complex<double>* b,
                      int*                  ldb,
                      std::complex<double>* alpha,
                      std::complex<double>* beta,
                      std::complex<double>* vl,
                      int*                  ldvl,
                      std::complex<double>* vr,
                      int*                  ldvr,
                      std::complex<double>* work,
                      int*                  lwork,
                      double*               rwork,
                      int*                  info);
    
}


namespace MAST {
    
    class LAPACK_ZGGEV{
        
    public:
        
        LAPACK_ZGGEV():
        info_val(-1)
        { }
        
        /*!
         *    computes the eigensolution for A x = \lambda B x. A & B will be
         *    overwritten
         */
        void compute(ComplexMatrixX& A,
                     ComplexMatrixX& B,
                     bool computeEigenvectors = true);
        
        ComputationInfo info() const;
        
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



#endif // __mast__lapack_zggev_interface_h__
