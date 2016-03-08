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
#include "aeroelasticity/ug_flutter_root.h"


MAST::UGFlutterRoot::UGFlutterRoot():
MAST::FlutterRootBase ()
{ }



MAST::UGFlutterRoot::UGFlutterRoot(const MAST::UGFlutterRoot& f):
MAST::FlutterRootBase (f)
{ }




void
MAST::UGFlutterRoot::init(const Real kr_ref_val,
                          const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left) {
    
    kr = kr_ref_val;
    
    if (std::abs(den) > 0.) {
        
        root = num/den;
        if (std::real(root) > 0.) {
            
            V     = sqrt(1./std::real(root));
            g     = std::imag(root)/std::real(root);
            omega = kr*V/b_ref;
            if_nonphysical_root = false;
        }
        else {
            
            V     = 0.;
            g     = 0.;
            omega = 0.;
            if_nonphysical_root = true;
        }
    }
    
    // calculate the modal participation vector
    const unsigned int
    nvals          = (int)Bmat.rows();
    eig_vec_right  = evec_right;
    eig_vec_left   = evec_left;
    ComplexVectorX
    k_q            = Bmat * evec_right;
    modal_participation.resize(nvals, 1);
    
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(evec_right(i)) * k_q(i));
    
    modal_participation *= (1./modal_participation.sum());
}


