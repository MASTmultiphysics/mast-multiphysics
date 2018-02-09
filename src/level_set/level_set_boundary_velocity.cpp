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
#include "level_set/level_set_boundary_velocity.h"


MAST::LevelSetBoundaryVelocity::LevelSetBoundaryVelocity(const unsigned int dim):
MAST::FieldFunction<RealVectorX>("phi"),
_phi(nullptr) {
    
}

MAST::LevelSetBoundaryVelocity::~LevelSetBoundaryVelocity() {
    
    if (_phi)
        delete _phi;
}

void
MAST::LevelSetBoundaryVelocity::init(MAST::SystemInitialization& sys,
                                     const libMesh::NumericVector<Real>& sol,
                                     const libMesh::NumericVector<Real>& dsol) {
    
    if (!_phi)
        _phi = new MAST::MeshFieldFunction(sys, "phi");
    else
        _phi->clear();
    _phi->init(sol, &dsol);
}


void
MAST::LevelSetBoundaryVelocity::operator() (const libMesh::Point& p,
                                            const Real t,
                                            RealVectorX& v) const {
    
    libmesh_assert(_phi);
    
    Real
    tol    = 1.e-6;
    
    RealVectorX
    val      = RealVectorX::Zero(1),
    grad     = RealVectorX::Zero(_dim),
    dval     = RealVectorX::Zero(1);
    
    RealMatrixX
    gradmat  = RealMatrixX::Zero(1, _dim);
    
    // the velocity is identified using the level set function gradient
    // and its sensitivity
    (*_phi)(p, t, val);
    
    // since the function provides the velocity at the boundary, the
    // level set value should be close to zero.
    libmesh_assert_less_equal(std::fabs(val(0)), tol);
    
    _phi->gradient(p, t, gradmat);
    _phi->perturbation(p, t, dval);
    // copy the matrix to the vector for further calculations
    grad = gradmat.row(0);
    
    // at boundary, phi(x) = 0
    // so,  dphi/dp + grad(phi) . V = 0
    //      dphi/dp + grad(phi) . (-grad(phi)  Vn) = 0   [since V = -grad(phi) Vn]
    //      dphi/dp -(grad(phi) .   grad(phi)) Vn  = 0
    //      Vn  =  (dphi/dp) / |grad(phi)|^2
    //      V   =  -grad(phi) Vn = -grad(phi) * (dphi/dp) / |grad(phi)|^2
    v = -dval(0)/grad.squaredNorm()*grad;
}

