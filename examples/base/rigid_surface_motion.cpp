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
#include "examples/base/rigid_surface_motion.h"
#include "aeroelasticity/frequency_function.h"



MAST::RigidSurfaceMotion::RigidSurfaceMotion():
_freq(nullptr),
_plunge_amplitude(0.),
_pitch_amplitude(0.),
_pitch_phase(0.) {
    
}


MAST::RigidSurfaceMotion::~RigidSurfaceMotion() {
    
}




void
MAST::RigidSurfaceMotion::init(MAST::FrequencyFunction& freq,
                               const RealVectorX& plunge_vector,
                               const RealVectorX& pitch_axis,
                               const RealVectorX& hinge_location,
                               const Real plunge_amplitude,
                               const Real pitch_amplitude,
                               const Real pitch_phase) {

    // make sure that the object has been cleared before setting the data again
    libmesh_assert(!_freq);
    
    _freq             = &freq;
    _plunge_vector    = plunge_vector;
    _pitch_axis       = pitch_axis;
    _hinge_location   = hinge_location;
    _plunge_amplitude = plunge_amplitude;
    _pitch_amplitude  = pitch_amplitude;
    _pitch_phase      = pitch_phase;
}



void
MAST::RigidSurfaceMotion::time_domain_motion(const Real t,
                                             const libMesh::Point& p,
                                             const libMesh::Point& n,
                                             RealVectorX& wdot,
                                             RealVectorX& dn_rot) {
    
    wdot.setZero();
    dn_rot.setZero();
    
    RealVector3
    r     = RealVector3::Zero(3),
    rot   = RealVector3::Zero(3);
    
    Real
    freq  = 0.;
    
    // get the value of frequency
    (*_freq)(freq);
    
    // include pitch-based displacement
    if (fabs(_pitch_amplitude) > 0.) {
        
        // normal distance from pitching axis to the given point
        for (unsigned int i=0; i<3; i++) r(i) = p(i);
        r      -= _hinge_location;
        
        rot     = _pitch_axis.cross(r);  // omega x r
        wdot    = rot * _pitch_amplitude * freq * cos(freq*t + _pitch_phase);
        
        for (unsigned int i=0; i<3; i++) r(i) = n(i);
        rot     = _pitch_axis.cross(r);  // omega x r
        dn_rot  = rot * _pitch_amplitude * sin(freq*t + _pitch_phase);
    }
    
    // add the plunge-based displacement
    if (fabs(_plunge_amplitude) > 0.)
        wdot.topRows(3)   += _plunge_vector * _plunge_amplitude * freq * cos(freq*t);
}




void
MAST::RigidSurfaceMotion::freq_domain_motion(const libMesh::Point& p,
                                             const libMesh::Point& n,
                                             ComplexVectorX& w,
                                             ComplexVectorX& dn_rot) {
    
    w.setZero();
    dn_rot.setZero();
    
    RealVector3
    r     = RealVector3::Zero(3),
    rot   = RealVector3::Zero(3);
    
    
    // include pitch-based displacement
    if (fabs(_pitch_amplitude) > 0.) {
        
        // normal distance from pitching axis to the given point
        Complex
        pitch_scale =
        Complex(cos(_pitch_phase), sin(_pitch_phase)) * _pitch_amplitude;
        
        for (unsigned int i=0; i<3; i++) r(i) = p(i);
        r      -= _hinge_location;
        
        rot     = _pitch_axis.cross(r);  // omega x r
        w       = rot.cast<Complex>() * pitch_scale;
        
        for (unsigned int i=0; i<3; i++) r(i) = n(i);
        rot     = _pitch_axis.cross(r);  // omega x r
        dn_rot  = rot.cast<Complex>() * pitch_scale;
    }
    
    // add the plunge-based displacement
    if (fabs(_plunge_amplitude) > 0.)
        w.topRows(3)    += _plunge_vector.cast<Complex>() * _plunge_amplitude;
}


