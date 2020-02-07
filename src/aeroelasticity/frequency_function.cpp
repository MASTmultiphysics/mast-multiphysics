/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "aeroelasticity/frequency_function.h"
#include "base/parameter.h"


MAST::FrequencyFunction::FrequencyFunction(const std::string& nm,
                                           MAST::FieldFunction<Real>& omega,
                                           MAST::FieldFunction<Real>& velocity,
                                           MAST::FieldFunction<Real>& b_ref):
MAST::FieldFunction<Real>(nm),
_if_red_freq(false),
_omega(omega),
_velocity(velocity),
_b_ref(b_ref) {
    
    _functions.insert(&_omega);
    _functions.insert(&_velocity);
    _functions.insert(&_b_ref);
}



MAST::FrequencyFunction::~FrequencyFunction() {
    
}



void
MAST::FrequencyFunction::operator() (Real& v) const {
    
    Real b, V;
    
    _omega(v);
    
    /*if (_if_red_freq) {
        
        _b_ref   (b);
        _velocity(V);
        
        v *= b/V;
    }*/
}




void
MAST::FrequencyFunction::derivative (    const MAST::FunctionBase& f,
                                     Real& v) const {
    
    Real w, b, vel, dw, db, dvel;
    
    _omega(w);  _omega.derivative( f, dw);
    v     = dw;
    
    /*if (_if_red_freq) {
        
        _b_ref     (b);  _b_ref.derivative   (d, f,   db);
        _velocity(vel);  _velocity.derivative( f, dvel);
        
        v *= b/vel;
        v += w * (db/vel - b/vel/vel*dvel);
    }*/
}




void
MAST::FrequencyFunction::nondimensionalizing_factor(Real& v) {
    
    if (!_if_red_freq)
        v = 1.;
    else {
        
        Real vel;
        _velocity(vel);
        _b_ref(v);
        v /= vel;
    }
}

