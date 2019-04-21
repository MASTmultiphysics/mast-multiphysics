/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "examples/base/multilinear_interpolation.h"


MAST::MultilinearInterpolation::MultilinearInterpolation
(const std::string& nm,
 std::map<Real, MAST::FieldFunction<Real>*>& values):
MAST::FieldFunction<Real>(nm),
_values(values) {
    
    // make sure that the size of the provided values is finite
    libmesh_assert(values.size() > 0);
    
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it = values.begin(), end = values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        _functions.insert(it->second);
}





MAST::MultilinearInterpolation::~MultilinearInterpolation() {
    
}



void
MAST::MultilinearInterpolation::operator() (const libMesh::Point& p,
                                            Real t,
                                            Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        (*it1->second)(p, t, f0);
        (*it2->second)(p, t, f1);
        // now interpolate
        v =  (f0 +
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}



void
MAST::MultilinearInterpolation::derivative(          const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           Real t,
                                           Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        it1->second->derivative( f, p, t, f0);
        it2->second->derivative( f, p, t, f1);
        // now interpolate
        v =  (f0 +
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}




MAST::SectionOffset::SectionOffset(const std::string& nm,
                                   const MAST::FieldFunction<Real> &thickness,
                                   const Real scale):
MAST::FieldFunction<Real>(nm),
_dim(thickness),
_scale(scale) {
    
    _functions.insert(&thickness);
}



MAST::SectionOffset::~SectionOffset() { }



void
MAST::SectionOffset::operator() (const libMesh::Point& p,
                                 Real t,
                                 Real& v) const {
    
    _dim(p, t, v);
    v *= 0.5*_scale;
}


void
MAST::SectionOffset::derivative(   const MAST::FunctionBase& f,
                                const libMesh::Point& p,
                                Real t,
                                Real& v) const {
    _dim.derivative( f, p, t, v);
    v *= 0.5*_scale;
}

