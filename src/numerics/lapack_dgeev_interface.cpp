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


// MAST includes
#include "numerics/lapack_dgeev_interface.h"


void
MAST::LAPACK_DGEEV::compute(const RealMatrixX &A,
                            bool computeEigenvectors) {
    
    libmesh_assert(A.cols() == A.rows());
    
    _A = A;
    
    RealMatrixX
    Amat = A;
    
    int n = (int)A.cols();
    
    char L='N',R='N';
    
    if (computeEigenvectors) {
        
        L = 'V'; R = 'V';
        VL.setZero(n, n);
        VR.setZero(n, n);
    }
    
    int
    lwork=16*n;
    
    info_val=-1;
    
    RealVectorX
    work,
    w_r,
    w_i;
    
    RealMatrixX
    vecl,
    vecr;
    
    W.setZero(n);
    work.setZero(lwork);
    w_r.setZero(n);
    w_i.setZero(n);
    vecl.setZero(n,n);
    vecr.setZero(n,n);
    
    Real
    *a_vals    = Amat.data(),
    *w_r_v     = w_r.data(),
    *w_i_v     = w_i.data(),
    *vecl_v    = vecl.data(),
    *vecr_v    = vecr.data(),
    *work_v    = work.data();
    
    
    dgeev_(&L, &R, &n,
           &(a_vals[0]), &n,
           &(w_r_v[0]), &(w_i_v[0]),
           &(vecl_v[0]), &n, &(vecr_v[0]), &n,
           &(work_v[0]), &lwork,
           &info_val);
    
    // now sort the eigenvalues for complex conjugates
    unsigned int n_located = 0;
    while (n_located < n) {
        
        // if the imaginary part of the eigenvalue is non-zero, it is a
        // complex conjugate
        if (w_i(n_located) != 0.) { // complex conjugate
            
            W(  n_located) = std::complex<double>(w_r(n_located),  w_i(n_located));
            W(1+n_located) = std::complex<double>(w_r(n_located), -w_i(n_located));
            
            // copy the eigenvectors if they were requested
            if (computeEigenvectors) {
                
                std::complex<double> iota = std::complex<double>(0, 1.);
                
                VL.col(  n_located) = (vecl.col(  n_located).cast<Complex>() +
                                       vecl.col(1+n_located).cast<Complex>() * iota);
                VL.col(1+n_located) = (vecl.col(  n_located).cast<Complex>() -
                                       vecl.col(1+n_located).cast<Complex>() * iota);
                VR.col(  n_located) = (vecr.col(  n_located).cast<Complex>() +
                                       vecr.col(1+n_located).cast<Complex>() * iota);
                VR.col(1+n_located) = (vecr.col(  n_located).cast<Complex>() -
                                       vecr.col(1+n_located).cast<Complex>() * iota);
            }
            
            // two complex conjugate roots were found
            n_located +=2;
        }
        else {
            
            W(  n_located) = std::complex<double>(w_r(n_located),  0.);
            
            // copy the eigenvectors if they were requested
            if (computeEigenvectors) {
                
                VL.col(n_located) = vecl.col(n_located).cast<Complex>();
                VR.col(n_located) = vecr.col(n_located).cast<Complex>();
            }
            
            // only one real root was found
            n_located++;
        }
    }
    
    if (info_val  != 0)
        libMesh::out
        << "Warning!!  DGEEV returned with nonzero info = "
        << info_val << std::endl;
}


