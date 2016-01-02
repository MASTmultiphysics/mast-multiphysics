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
#include "coordinates/basis_matrix_coordinate.h"


MAST::BasisMatrixCoordinate::BasisMatrixCoordinate(const std::string& nm,
                                                   MAST::FieldFunction<RealMatrixX>* basis):
MAST::CoordinateBase(nm),
_basis(basis) {
    
    _functions.insert(basis->master());
}




MAST::BasisMatrixCoordinate::~BasisMatrixCoordinate() {
    
    delete _basis;
}




MAST::BasisMatrixCoordinate::BasisMatrixCoordinate(const MAST::BasisMatrixCoordinate& c):
MAST::CoordinateBase(c),
_basis(c._basis->clone().release()) {
    
    _functions.insert(_basis->master());
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::BasisMatrixCoordinate::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    (new MAST::BasisMatrixCoordinate(*this));
}




void
MAST::BasisMatrixCoordinate::operator() (const libMesh::Point& p,
                                         const Real t,
                                         RealMatrixX& v) const {
    
    (*_basis)(p, t, v);
}




void
MAST::BasisMatrixCoordinate::derivative (const MAST::DerivativeType d,
                                         const MAST::FunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         RealMatrixX& v) const {
    
    _basis->derivative(d, f, p, t, v);
}



