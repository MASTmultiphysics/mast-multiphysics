/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "base/real_output_function.h"


MAST::RealOutputFunction::RealOutputFunction(MAST::OutputQuantityType t):
MAST::OutputFunctionBase(t),
_value(0.) {
    
}



MAST::RealOutputFunction::~RealOutputFunction() {
    
}



void
MAST::RealOutputFunction::clear() {
    
    _value = 0.;
    _derivative.setZero();
    _sensitivity.clear();
}



void
MAST::RealOutputFunction::set_value(Real val) {
    
    _value = val;
}



void
MAST::RealOutputFunction::add_value(Real val) {
    
    _value += val;
}


Real
MAST::RealOutputFunction::get_value() const {
    
    return _value;
}



void
MAST::RealOutputFunction::set_derivative(const RealVectorX& dvdX) {
    
    _derivative = dvdX;
}



const RealVectorX&
MAST::RealOutputFunction::get_derivative() const {
    
    return _derivative;
}



void
MAST::RealOutputFunction::set_sensitivity(const MAST::FunctionBase* f,
                                          const Real dvdp) {
    
    // make sure that the function does not alreadey exists
    
    std::map<const MAST::FunctionBase*, Real>::const_iterator
    it   =  _sensitivity.find(f);
    
    libmesh_assert(it == _sensitivity.end());
    
    _sensitivity.insert(std::pair<const MAST::FunctionBase*, Real>(f, dvdp));
}



void
MAST::RealOutputFunction::add_sensitivity(const MAST::FunctionBase* f,
                                          const Real dvdp) {
    
    // make sure that the function exists
    
    std::map<const MAST::FunctionBase*, Real>::iterator
    it   =  _sensitivity.find(f);
    
    if(it == _sensitivity.end())
        _sensitivity.insert(std::pair<const MAST::FunctionBase*, Real>(f, dvdp));
    else
        it->second += dvdp;
}



Real
MAST::RealOutputFunction::get_sensitivity(const MAST::FunctionBase* f) const {

    // make sure that the function exists
    
    std::map<const MAST::FunctionBase*, Real>::const_iterator
    it   =  _sensitivity.find(f);
    
    libmesh_assert(it != _sensitivity.end());
    
    return it->second;
}



