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
#include "boundary_condition/function_surface_motion.h"
#include "base/system_initialization.h"
#include "aeroelasticity/frequency_function.h"
#include "aeroelasticity/displacement_function_base.h"



MAST::FunctionSurfaceMotion::FunctionSurfaceMotion(const std::string& nm):
MAST::SurfaceMotionBase(nm),
_freq(nullptr),
_w   (nullptr) {
    
}




MAST::FunctionSurfaceMotion::~FunctionSurfaceMotion() {
    
}







void
MAST::FunctionSurfaceMotion::init(MAST::FrequencyFunction&           freq,
                                  MAST::DisplacementFunctionBase&    w) {
    
    _freq             = &freq;
    _w                = &w;
}





void
MAST::FunctionSurfaceMotion::time_domain_motion(const Real t,
                                                const libMesh::Point& p,
                                                const libMesh::Point& n,
                                                RealVectorX& wdot,
                                                RealVectorX& dn_rot) {
    
    
    wdot.setZero();
    dn_rot.setZero();
    
    RealVectorX
    w      = RealVectorX::Zero(3),
    dwdx   = RealVectorX::Zero(3),
    dwdy   = RealVectorX::Zero(3),
    dwdz   = RealVectorX::Zero(3),
    rot    = RealVectorX::Zero(3);
    
    Real
    freq   = 0.;
    
    (*_freq)(freq);
    (*_w)   (p, 0.,    w);
    _w->dwdx(p, 0., dwdx);
    _w->dwdy(p, 0., dwdy);
    _w->dwdz(p, 0., dwdz);
    
    
    // now copy the values to u_trans
    wdot  = w * freq * cos(freq*t);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    rot(0) = dwdy(2) - dwdz(1); // dwz/dy - dwy/dz
    rot(1) = dwdz(0) - dwdx(2); // dwx/dz - dwz/dx
    rot(2) = dwdx(1) - dwdy(0); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}




void
MAST::FunctionSurfaceMotion::freq_domain_motion(const libMesh::Point& p,
                                                const libMesh::Point& n,
                                                ComplexVectorX& w,
                                                ComplexVectorX& dn_rot) {
    
    
    w.setZero();
    dn_rot.setZero();
    
    RealVectorX
    w_R    = RealVectorX::Zero(3),
    dwdx   = RealVectorX::Zero(3),
    dwdy   = RealVectorX::Zero(3),
    dwdz   = RealVectorX::Zero(3),
    rot    = RealVectorX::Zero(3);
    
    Real
    freq   = 0.;
    
    (*_freq)(freq);
    (*_w)   (p, 0.,  w_R);
    _w->dwdx(p, 0., dwdx);
    _w->dwdy(p, 0., dwdy);
    _w->dwdz(p, 0., dwdz);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    w    =  w_R.cast<Complex>();
    
    // now prepare the rotation vector
    rot(0) = dwdy(2) - dwdz(1); // dwz/dy - dwy/dz
    rot(1) = dwdz(0) - dwdx(2); // dwx/dz - dwz/dx
    rot(2) = dwdx(1) - dwdy(0); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}





