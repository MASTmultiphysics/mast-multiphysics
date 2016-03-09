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
    
    _functions.insert(_omega.master());
    _functions.insert(_velocity.master());
    _functions.insert(_b_ref.master());
}


MAST::FrequencyFunction::FrequencyFunction(const MAST::FrequencyFunction& f):
MAST::FieldFunction<Real>(f),
_if_red_freq(f._if_red_freq),
_omega(f._omega),
_velocity(f._velocity),
_b_ref(f._b_ref) {
    
    _functions.insert(_omega.master());
    _functions.insert(_velocity.master());
    _functions.insert(_b_ref.master());
}


MAST::FrequencyFunction::~FrequencyFunction() {
    
}





std::auto_ptr<MAST::FieldFunction<Real> >
MAST::FrequencyFunction::clone() const {
    
    MAST::FrequencyFunction* rval = new MAST::FrequencyFunction(*this);
    
    return std::auto_ptr<MAST::FieldFunction<Real> >(rval);
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
MAST::FrequencyFunction::derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     Real& v) const {
    
    Real w, b, vel, dw, db, dvel;
    
    _omega(w);  _omega.derivative(d, f, dw);
    v     = dw;
    
    /*if (_if_red_freq) {
        
        _b_ref     (b);  _b_ref.derivative   (d, f,   db);
        _velocity(vel);  _velocity.derivative(d, f, dvel);
        
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

