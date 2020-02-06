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
#include "coordinates/polar_coordinate.h"


MAST::PolarCoordinate::PolarCoordinate(const std::string& nm,
                                       MAST::FieldFunction<Real>& theta):
MAST::CoordinateBase(nm),
_theta(theta) {
    
    _functions.insert(&theta);
}


MAST::PolarCoordinate::~PolarCoordinate() { }



void
MAST::PolarCoordinate::operator() (const libMesh::Point& p,
                                   const Real t,
                                   RealMatrixX& v) const {
    v.setZero(3,3);
    Real theta;
    _theta(p,t,theta);
    v(0,0) = cos(theta);
    v(1,0) = sin(theta);
    
    v(0,1) = -sin(theta);
    v(1,1) =  cos(theta);
    
    v(2,2) = 1.;
}



void
MAST::PolarCoordinate::derivative (  const MAST::FunctionBase& f,
                                   const libMesh::Point& p,
                                   const Real t,
                                   RealMatrixX& v) const {
    v.setZero(3,3);
    Real theta, dtheta;
    _theta(p,t,theta);
    _theta.derivative(f,p,t,dtheta);
    
    v(0,0) = -sin(theta)*dtheta;
    v(1,0) =  cos(theta)*dtheta;
    
    v(0,1) = -cos(theta)*dtheta;
    v(1,1) = -sin(theta)*dtheta;
}
