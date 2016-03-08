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
#include "boundary_condition/flexible_surface_motion.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"


MAST::FlexibleSurfaceMotion::FlexibleSurfaceMotion(MAST::SystemInitialization& sys):
MAST::SurfaceMotionBase(),
_system(sys),
_freq(NULL) {
    
}




MAST::FlexibleSurfaceMotion::~FlexibleSurfaceMotion() {
    
}







void
MAST::FlexibleSurfaceMotion::init(MAST::FrequencyFunction& freq,
                                  libMesh::NumericVector<Real>& sol) {
    
    // make sure that the object has been cleared before setting the data again
    //libmesh_assert(!_freq);
    
    libMesh::System& sys = _system.system();

    _freq             = &freq;

    // first initialize the solution to the given vector
    _sol.reset(libMesh::NumericVector<Real>::build(sys.comm()).release());
    _sol->init(sol.size(), true, libMesh::SERIAL);
    
    // now localize the give solution to this objects's vector
    sol.localize(*_sol);
    
    // if the mesh function has not been created so far, initialize it
    // we use only the u, v, w variables here.
    const std::vector<unsigned int>& sys_vars = _system.vars();
    std::vector<unsigned int> vars(3);
    for (unsigned int i=0; i<3; i++) vars[i] = sys_vars[i];
    
    _function.reset(new libMesh::MeshFunction(sys.get_equation_systems(),
                                              *_sol,
                                              sys.get_dof_map(),
                                              vars));
    _function->init();
}





void
MAST::FlexibleSurfaceMotion::freq_domain_motion(const libMesh::Point& p,
                                                const libMesh::Point& n,
                                                ComplexVectorX& w,
                                                ComplexVectorX& dn_rot) {
    
    
    w.setZero();
    dn_rot.setZero();
    
    libmesh_assert(_function.get()); // should be initialized before this call
    
    // translation is obtained by direct interpolation of the u,v,w vars
    
    DenseRealVector v(3);
    (*_function)(p, 0., v);
    
    
    // now copy the values to u_trans
    Complex iota(0., 1.);
    for (unsigned int i=0; i<3; i++) w(i) = v(i);
    
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient> gradients;
    _function->gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    ComplexVectorX
    rot = ComplexVectorX::Zero(3);
    
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}





