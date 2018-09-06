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


// MAST includes
#include "aeroelasticity/time_domain_flutter_root.h"


MAST::TimeDomainFlutterRoot::TimeDomainFlutterRoot():
MAST::FlutterRootBase()
{ }



MAST::TimeDomainFlutterRoot::TimeDomainFlutterRoot(const TimeDomainFlutterRoot& f):
MAST::FlutterRootBase(f)
{ }



void
MAST::TimeDomainFlutterRoot::init(const Real v_ref_val,
                                  const Complex num,
                                  const Complex den,
                                  const RealMatrixX& Bmat,
                                  const ComplexVectorX& evec_right,
                                  const ComplexVectorX& evec_left) {
    
    if (std::abs(den) > 0.) {
        
        root                = num/den;
        V                   = v_ref_val;
        omega               = std::imag(root);
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

