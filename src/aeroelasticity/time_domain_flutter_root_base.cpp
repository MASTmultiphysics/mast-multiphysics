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
#include "aeroelasticity/time_domain_flutter_root_base.h"


MAST::TimeDomainFlutterRootBase::TimeDomainFlutterRootBase():
has_sensitivity_data  (false),
V                     (0.),
g                     (0.),
omega                 (0.),
V_sens                (0.),
root                  (0.),
root_sens             (0.)
{ }



MAST::TimeDomainFlutterRootBase::TimeDomainFlutterRootBase(const TimeDomainFlutterRootBase& f):
has_sensitivity_data   (f.has_sensitivity_data),
V                      (f.V),
g                      (f.g),
omega                  (f.omega),
V_sens                 (f.V_sens),
root                   (f.root),
root_sens              (f.root_sens),
eig_vec_right          (f.eig_vec_right),
eig_vec_left           (f.eig_vec_left),
modal_participation    (f.modal_participation)
{ }



void
MAST::TimeDomainFlutterRootBase::copy_root(const MAST::TimeDomainFlutterRootBase& f) {
    
    has_sensitivity_data   = f.has_sensitivity_data;
    V                      = f.V;
    g                      = f.g;
    omega                  = f.omega;
    V_sens                 = f.V_sens;
    root                   = f.root;
    root_sens              = f.root_sens;
    eig_vec_right          = f.eig_vec_right;
    eig_vec_left           = f.eig_vec_left;
    modal_participation    = f.modal_participation;
}


void
MAST::TimeDomainFlutterRootBase::init(const Real v_ref_val,
                            const Complex num,
                            const Complex den,
                            const RealMatrixX& Bmat,
                            const ComplexVectorX& evec_right,
                            const ComplexVectorX& evec_left) {
    
    if (std::abs(den) > 0.)
    {
        root                = num/den;
        V                   = v_ref_val;
        g                   = std::real(root);
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

