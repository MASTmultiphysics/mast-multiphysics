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

#ifndef __mast__lapack_zggevx_interface_h__
#define __mast__lapack_zggevx_interface_h__


// MAST includes
#include "base/mast_data_types.h"
#include "numerics/lapack_zggev_base.h"


extern "C" {
    /*
     *  =====================================================================
     *  Purpose
     *  =======
     *
     *  ZGGEVX computes for a pair of N-by-N complex nonsymmetric matrices
     *  (A,B) the generalized eigenvalues, and optionally, the left and/or
     *  right generalized eigenvectors.
     *
     *  Optionally, it also computes a balancing transformation to improve
     *  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     *  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     *  the eigenvalues (RCONDE), and reciprocal condition numbers for the
     *  right eigenvectors (RCONDV).
     *
     *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     *  singular. It is usually represented as the pair (alpha,beta), as
     *  there is a reasonable interpretation for beta=0, and even for both
     *  being zero.
     *
     *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     *  of (A,B) satisfies
     *                   A * v(j) = lambda(j) * B * v(j) .
     *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     *  of (A,B) satisfies
     *                   u(j)**H * A  = lambda(j) * u(j)**H * B.
     *  where u(j)**H is the conjugate-transpose of u(j).
     *
     *
     *  Arguments
     *  =========
     *
     *  BALANC  (input) CHARACTER*1
     *          Specifies the balance option to be performed:
     *          = 'N':  do not diagonally scale or permute;
     *          = 'P':  permute only;
     *          = 'S':  scale only;
     *          = 'B':  both permute and scale.
     *          Computed reciprocal condition numbers will be for the
     *          matrices after permuting and/or balancing. Permuting does
     *          not change condition numbers (in exact arithmetic), but
     *          balancing does.
     *
     *  JOBVL   (input) CHARACTER*1
     *          = 'N':  do not compute the left generalized eigenvectors;
     *          = 'V':  compute the left generalized eigenvectors.
     *
     *  JOBVR   (input) CHARACTER*1
     *          = 'N':  do not compute the right generalized eigenvectors;
     *          = 'V':  compute the right generalized eigenvectors.
     *
     *  SENSE   (input) CHARACTER*1
     *          Determines which reciprocal condition numbers are computed.
     *          = 'N': none are computed;
     *          = 'E': computed for eigenvalues only;
     *          = 'V': computed for eigenvectors only;
     *          = 'B': computed for eigenvalues and eigenvectors.
     *
     *  N       (input) INTEGER
     *          The order of the matrices A, B, VL, and VR.  N >= 0.
     *
     *  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
     *          On entry, the matrix A in the pair (A,B).
     *          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
     *          or both, then A contains the first part of the complex Schur
     *          form of the "balanced" versions of the input A and B.
     *
     *  LDA     (input) INTEGER
     *          The leading dimension of A.  LDA >= max(1,N).
     *
     *  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
     *          On entry, the matrix B in the pair (A,B).
     *          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
     *          or both, then B contains the second part of the complex
     *          Schur form of the "balanced" versions of the input A and B.
     *
     *  LDB     (input) INTEGER
     *          The leading dimension of B.  LDB >= max(1,N).
     *
     *  ALPHA   (output) COMPLEX*16 array, dimension (N)
     *  BETA    (output) COMPLEX*16 array, dimension (N)
     *          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized
     *          eigenvalues.
     *
     *          Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or
     *          underflow, and BETA(j) may even be zero.  Thus, the user
     *          should avoid naively computing the ratio ALPHA/BETA.
     *          However, ALPHA will be always less than and usually
     *          comparable with norm(A) in magnitude, and BETA always less
     *          than and usually comparable with norm(B).
     *
     *  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
     *          If JOBVL = 'V', the left generalized eigenvectors u(j) are
     *          stored one after another in the columns of VL, in the same
     *          order as their eigenvalues.
     *          Each eigenvector will be scaled so the largest component
     *          will have abs(real part) + abs(imag. part) = 1.
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
     *          Each eigenvector will be scaled so the largest component
     *          will have abs(real part) + abs(imag. part) = 1.
     *          Not referenced if JOBVR = 'N'.
     *
     *  LDVR    (input) INTEGER
     *          The leading dimension of the matrix VR. LDVR >= 1, and
     *          if JOBVR = 'V', LDVR >= N.
     *
     *  ILO     (output) INTEGER
     *  IHI     (output) INTEGER
     *          ILO and IHI are integer values such that on exit
     *          A(i,j) = 0 and B(i,j) = 0 if i > j and
     *          j = 1,...,ILO-1 or i = IHI+1,...,N.
     *          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
     *
     *  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
     *          Details of the permutations and scaling factors applied
     *          to the left side of A and B.  If PL(j) is the index of the
     *          row interchanged with row j, and DL(j) is the scaling
     *          factor applied to row j, then
     *            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
     *                      = DL(j)  for j = ILO,...,IHI
     *                      = PL(j)  for j = IHI+1,...,N.
     *          The order in which the interchanges are made is N to IHI+1,
     *          then 1 to ILO-1.
     *
     *  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
     *          Details of the permutations and scaling factors applied
     *          to the right side of A and B.  If PR(j) is the index of the
     *          column interchanged with column j, and DR(j) is the scaling
     *          factor applied to column j, then
     *            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
     *                      = DR(j)  for j = ILO,...,IHI
     *                      = PR(j)  for j = IHI+1,...,N
     *          The order in which the interchanges are made is N to IHI+1,
     *          then 1 to ILO-1.
     *
     *  ABNRM   (output) DOUBLE PRECISION
     *          The one-norm of the balanced matrix A.
     *
     *  BBNRM   (output) DOUBLE PRECISION
     *          The one-norm of the balanced matrix B.
     *
     *  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
     *          If SENSE = 'E' or 'B', the reciprocal condition numbers of
     *          the eigenvalues, stored in consecutive elements of the array.
     *          If SENSE = 'N' or 'V', RCONDE is not referenced.
     *
     *  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
     *          If JOB = 'V' or 'B', the estimated reciprocal condition
     *          numbers of the eigenvectors, stored in consecutive elements
     *          of the array. If the eigenvalues cannot be reordered to
     *          compute RCONDV(j), RCONDV(j) is set to 0; this can only occur
     *          when the true value would be very small anyway.
     *          If SENSE = 'N' or 'E', RCONDV is not referenced.
     *
     *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
     *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     *
     *  LWORK   (input) INTEGER
     *          The dimension of the array WORK. LWORK >= max(1,2*N).
     *          If SENSE = 'E', LWORK >= max(1,4*N).
     *          If SENSE = 'V' or 'B', LWORK >= max(1,2*N*N+2*N).
     *
     *          If LWORK = -1, then a workspace query is assumed; the routine
     *          only calculates the optimal size of the WORK array, returns
     *          this value as the first entry of the WORK array, and no error
     *          message related to LWORK is issued by XERBLA.
     *
     *  RWORK   (workspace) REAL array, dimension (lrwork)
     *          lrwork must be at least max(1,6*N) if BALANC = 'S' or 'B',
     *          and at least max(1,2*N) otherwise.
     *          Real workspace.
     *
     *  IWORK   (workspace) INTEGER array, dimension (N+2)
     *          If SENSE = 'E', IWORK is not referenced.
     *
     *  BWORK   (workspace) LOGICAL array, dimension (N)
     *          If SENSE = 'N', BWORK is not referenced.
     *
     *  INFO    (output) INTEGER
     *          = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value.
     *          = 1,...,N:
     *                The QZ iteration failed.  No eigenvectors have been
     *                calculated, but ALPHA(j) and BETA(j) should be correct
     *                for j=INFO+1,...,N.
     *          > N:  =N+1: other than QZ iteration failed in ZHGEQZ.
     *                =N+2: error return from ZTGEVC.
     *
     *  Further Details
     *  ===============
     *
     *  Balancing a matrix pair (A,B) includes, first, permuting rows and
     *  columns to isolate eigenvalues, second, applying diagonal similarity
     *  transformation to the rows and columns to make the rows and columns
     *  as close in norm as possible. The computed reciprocal condition
     *  numbers correspond to the balanced matrix. Permuting rows and columns
     *  will not change the condition numbers (in exact arithmetic) but
     *  diagonal scaling will.  For further explanation of balancing, see
     *  section 4.11.1.2 of LAPACK Users' Guide.
     *
     *  An approximate error bound on the chordal distance between the i-th
     *  computed generalized eigenvalue w and the corresponding exact
     *  eigenvalue lambda is
     *
     *       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
     *
     *  An approximate error bound for the angle between the i-th computed
     *  eigenvector VL(i) or VR(i) is given by
     *
     *       EPS * norm(ABNRM, BBNRM) / DIF(i).
     *
     *  For further explanation of the reciprocal condition numbers RCONDE
     *  and RCONDV, see section 4.11 of LAPACK User's Guide.
     *
     *  =====================================================================
     */
    extern int zggevx_(char*                 balanc,
                       char*                 jobvl,
                       char*                 jobvr,
                       char*                 sens,
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
                       int*                  ilo,
                       int*                  ihi,
                       double*               lscale,
                       double*               rscale,
                       double*               abnrm,
                       double*               bbnrm,
                       double*               rconde,
                       double*               rcondv,
                       std::complex<double>* work,
                       int*                  lwork,
                       double*               rwork,
                       int*                  iwork,
                       bool*                 bwork,
                       int*                  info);
    
}


namespace MAST {
    
    class LAPACK_ZGGEVX:
    public MAST::LAPACK_ZGGEV_Base {
        
    public:
        
        LAPACK_ZGGEVX():
        MAST::LAPACK_ZGGEV_Base()
        { }
        
        /*!
         *    computes the eigensolution for A x = \lambda B x. A & B will be
         *    overwritten
         */
        virtual void compute(const ComplexMatrixX& A,
                             const ComplexMatrixX& B,
                             bool computeEigenvectors = true);
        
    };
}



#endif // __mast__lapack_zggevx_interface_h__
