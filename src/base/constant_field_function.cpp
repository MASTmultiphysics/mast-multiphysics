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
#include "base/constant_field_function.h"
#include "base/parameter.h"



MAST::ConstantFieldFunction::
ConstantFieldFunction(const std::string& nm,
                      const MAST::Parameter& p):
MAST::FieldFunction<Real>(nm),
_p(p) {
    
    _functions.insert(p.master());
}


MAST::ConstantFieldFunction::~ConstantFieldFunction() {
    
}



MAST::ConstantFieldFunction::ConstantFieldFunction(const MAST::ConstantFieldFunction& f):
MAST::FieldFunction<Real>(f),
_p(f._p) {
    
    _functions.insert(f.master());
}



std::auto_ptr<MAST::FieldFunction<Real> >
MAST::ConstantFieldFunction::clone() const {
    MAST::FieldFunction<Real>* rval =
    new MAST::ConstantFieldFunction(*this);
    
    return std::auto_ptr<MAST::FieldFunction<Real> >(rval);
}



void
MAST::ConstantFieldFunction::operator() (const libMesh::Point& p,
                                         const Real t,
                                         Real& v) const {
    v = _p();
}




void
MAST::ConstantFieldFunction::derivative (const MAST::DerivativeType d,
                                         const MAST::FunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         Real& v) const {
    
    v = _p.depends_on(f)?1:0;
}

