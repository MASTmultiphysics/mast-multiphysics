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
#include "examples/fsi/beam_flag_flutter_solution/beam_flag_frequency_domain_displacement.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"


MAST::BeamFlagFrequencyDomainDisplacement::
BeamFlagFrequencyDomainDisplacement(MAST::SystemInitialization& sys,
                                    const std::string& nm):
MAST::ComplexMeshFieldFunction(sys, nm)
{ }




MAST::BeamFlagFrequencyDomainDisplacement::~BeamFlagFrequencyDomainDisplacement() {
    
}




void
MAST::BeamFlagFrequencyDomainDisplacement::operator() (const libMesh::Point& p,
                                                       const Real t,
                                                       ComplexVectorX& v) const {
    
    // make sure that the object was initialized
    libmesh_assert(_function_re);
    
    libMesh::Point
    pt    = p;
    pt(1) = 0;   // since the beam is assumed to be at y=0
    
    DenseRealVector v_re, v_im;
    (*_function_re)(pt, t, v_re);
    (*_function_im)(pt, t, v_im);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v_re.size());
    
    // now copy this to the output vector
    v = ComplexVectorX::Zero(v_re.size());
    for (unsigned int i=0; i<v_re.size(); i++)
        v(i) = std::complex<Real>(v_re(i), v_im(i));
}





void
MAST::BeamFlagFrequencyDomainDisplacement::perturbation(const libMesh::Point& p,
                                                        const Real t,
                                                        ComplexVectorX& v) const {
    
    // make sure that the object was initialized
    libmesh_assert(_perturbed_function_re);
    
    libMesh::Point
    pt    = p;
    pt(1) = 0;   // since the beam is assumed to be at y=0

    DenseRealVector v_re, v_im;
    (*_perturbed_function_re)(pt, t, v_re);
    (*_perturbed_function_im)(pt, t, v_im);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v_re.size());
    
    // now copy this to the output vector
    v = ComplexVectorX::Zero(v_re.size());
    for (unsigned int i=0; i<v_re.size(); i++)
        v(i) = std::complex<Real>(v_re(i), v_im(i));
}



