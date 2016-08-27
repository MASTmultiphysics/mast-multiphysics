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


// MAST includes
#include "numerics/lapack_zggevx_interface.h"


void
MAST::LAPACK_ZGGEVX::compute(const ComplexMatrixX &A,
                             const ComplexMatrixX &B,
                             bool computeEigenvectors) {
    
    libmesh_assert(A.cols() == A.rows() &&
                   B.cols() == A.rows() &&
                   B.cols() == B.rows());
    
    _A = A;
    _B = B;
    
    ComplexMatrixX
    Amat = _A,
    Bmat = _B;
    
    int n = (int)A.cols();
    
    char BAL='B', L='N',R='N', S='E';
    
    if (computeEigenvectors) {
        
        L = 'V'; R = 'V'; S='B';
        VL.setZero(n, n);
        VR.setZero(n, n);
    }
    
    int
    lwork    =4*(n*n+n),
    l_rwork  =8*n,
    ilo      = 0,
    ihi      = 0;
    info_val =-1,
    
    alpha.setZero(n);
    beta.setZero(n);
    ComplexVectorX
    work   = ComplexVectorX::Zero(lwork);
    RealVectorX
    rwork  = RealVectorX::Zero(l_rwork),
    lscale = RealVectorX::Zero(n),
    rscale = RealVectorX::Zero(n),
    rconde = RealVectorX::Zero(n),
    rcondv = RealVectorX::Zero(n);
    
    Complex
    *a_vals  = Amat.data(),
    *b_vals  = Bmat.data(),
    *alpha_v = alpha.data(),
    *beta_v  = beta.data(),
    *VL_v    = VL.data(),
    *VR_v    = VR.data(),
    *work_v  = work.data();
    
    Real
    *rwork_v  = rwork.data(),
    *lscale_v = lscale.data(),
    *rscale_v = rscale.data(),
    *rconde_v = rconde.data(),
    *rcondv_v = rcondv.data(),
    abnrm     = 0.,
    bbnrm     = 0.;
    
    std::vector<int>
    iwork(n+2, 0);
    
    bool
    bwork[n];
    
    zggevx_(&BAL, &L, &R, &S, &n,
            &(a_vals[0]), &n,
            &(b_vals[0]), &n,
            &(alpha_v[0]), &(beta_v[0]),
            &(VL_v[0]), &n,
            &(VR_v[0]), &n,
            &ilo, &ihi,
            &(lscale_v[0]), &(rscale_v[0]),
            &abnrm, &bbnrm,
            &(rconde_v[0]), &(rcondv_v[0]),
            &(work_v[0]), &lwork,
            &(rwork_v[0]),
            &(iwork[0]),
            &(bwork[0]),
            &info_val);
    
    if (info_val  != 0)
        libMesh::out
        << "Warning!!  ZGGEVX returned with nonzero info = "
        << info_val << std::endl;
}


