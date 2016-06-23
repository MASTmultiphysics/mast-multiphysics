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

/*
// MAST includes
#include "aeroelasticity/pk_flutter_solution.h"
#include "aeroelasticity/pk_flutter_root.h"
#include "numerics/lapack_zggev_interface.h"


void
MAST::PKFlutterSolution::init (const MAST::PKFlutterSolver& solver,
                               const Real k_red,
                               const Real v_ref,
                               const Real bref,
                               const RealMatrixX& kmat,
                               const LAPACK_ZGGEV& eig_sol) {
    
    // make sure that it hasn't already been initialized
    libmesh_assert(!_roots.size());
    
    _ref_val     = v_ref;
    _k_red       = k_red;
    _stiff_mat   = kmat;
    _Bmat        = eig_sol.B();
    
    const ComplexMatrixX
    &VR          = eig_sol.right_eigenvectors(),
    &VL          = eig_sol.left_eigenvectors();
    const ComplexVectorX
    &num         = eig_sol.alphas(),
    &den         = eig_sol.betas();
    
    unsigned int
    nvals        = (unsigned int)_Bmat.rows();
    
    _roots.resize(nvals);
    
    // iterate over the roots and initialize the vector
    for (unsigned int i=0; i<nvals; i++) {
        
        MAST::PKFlutterRoot* root = new MAST::PKFlutterRoot;
        root->init(k_red,
                   v_ref,
                   bref,
                   num(i),
                   den(i),
                   kmat,
                   VR.col(i),
                   VL.col(i));
        
        _roots[i] = root;
    }
}

*/