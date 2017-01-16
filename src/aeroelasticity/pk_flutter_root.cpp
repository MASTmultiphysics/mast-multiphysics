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


// MAST includes
#include "aeroelasticity/pk_flutter_root.h"


MAST::PKFlutterRoot::PKFlutterRoot():
MAST::FlutterRootBase()
{ }



void
MAST::PKFlutterRoot::init(const Real k_red_val,
                          const Real v_ref_val,
                          const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const RealMatrixX& Kmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left) {
    
    kr      = k_red_val;
    V       = v_ref_val;
    
    if (std::abs(den) > 0.) {
        
        root = num/den;
        g     = std::real(root);
        omega = std::imag(root);
    }

    // calculate the modal participation vector
    // Since PK flutter solver works with complex matrices, we cannot look at
    // roots only in the upper or lower half of the complex domain. Here,
    // however, we only are interested in the modal participation vector.
    const unsigned int
    nvals         = (int)Kmat.rows();
    eig_vec_right = evec_right.topRows(nvals);
    eig_vec_left  = evec_left.topRows(nvals);
    
    // use the stiffness matrix and the first half of the eigenvectors to
    // calculate the modal participation based on strain energy
    
    ComplexVectorX
    k_q           = Kmat * eig_vec_right;
    modal_participation.resize(nvals, 1);
    
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(eig_vec_right(i)) * k_q(i));
    
    modal_participation *= (1./modal_participation.sum());
}



//bool
//MAST::PKFlutterRoot::is_similar(MAST::FlutterRootBase &r) const {
//    
//    const Real tol = 1.0e-6;
//    bool similar = false;
//    
//    // currently only the modal participation is used to check for
//    // similarity
//    RealVectorX
//    diff = (modal_participation - r.modal_participation) * (omega-r.omega);
//    
//    if (diff.norm() <= tol)
//        similar = true;
//    
//    return similar;
//}




