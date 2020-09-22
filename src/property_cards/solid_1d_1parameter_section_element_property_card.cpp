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
#include "property_cards/solid_1d_1parameter_section_element_property_card.h"
//#include "property_cards/material_property_card_base.h"
//#include "base/field_function_base.h"
//#include "base/elem_base.h"


MAST::Solid1D1ParameterSectionProperty::Area::Area(
    const std::function<void(Real&, Real&)> func,
    const std::function<void(Real&, Real&, Real&)> dfunc,
    const MAST::FieldFunction<Real>& DIM1):
    MAST::FieldFunction<Real>("Area"),
    _func(func), _dfunc(dfunc), _DIM1(DIM1)
    {
        _functions.insert(&DIM1);
    }
    
void MAST::Solid1D1ParameterSectionProperty::Area::operator() 
    (const libMesh::Point& p, const Real t, Real& m) const {
        Real DIM1;
        _DIM1(p, t, DIM1);           
        _func(DIM1, m);
    }
    
void MAST::Solid1D1ParameterSectionProperty::Area::derivative(
    const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
    Real& m) const {
        Real DIM1, dDIM1;
        _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
        _dfunc(DIM1, dDIM1, m);
    }
    

MAST::Solid1D1ParameterSectionProperty::TorsionalConstant::
TorsionalConstant(const std::function<void(Real&, Real&)> func,
                  const std::function<void(Real&, Real&, Real&)> dfunc,
                  const MAST::FieldFunction<Real>& DIM1):
                  MAST::FieldFunction<Real>("TorsionalConstant"),
                  _func(func), _dfunc(dfunc), _DIM1(DIM1){
    _functions.insert(&DIM1);
}
    

void MAST::Solid1D1ParameterSectionProperty::TorsionalConstant::
operator() (const libMesh::Point& p, const Real t, Real& m) const {
    Real DIM1;
    _DIM1(p, t, DIM1);
    _func(DIM1, m);
}

    
void MAST::Solid1D1ParameterSectionProperty::TorsionalConstant::
derivative (MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m) {
    Real DIM1, dDIM1;
    _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
    _dfunc(DIM1, dDIM1, m);
}
    

    
MAST::Solid1D1ParameterSectionProperty::PolarInertia::
PolarInertia(const std::function<void(Real&, Real&)> func,
                         const std::function<void(Real&, Real&, Real&)> dfunc,
                         const std::function<void(Real&, Real&)> Afunc,
                         const std::function<void(Real&, Real&, Real&)> dAfunc,
                         const MAST::FieldFunction<Real>& DIM1,
                         const MAST::FieldFunction<Real>&  hy_offset,
                         const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("PolarInertia"),
            _func(func),
            _dfunc(dfunc),
            _Afunc(Afunc),
            _dAfunc(dAfunc),
            _DIM1(DIM1),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&DIM1);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }

void MAST::Solid1D1ParameterSectionProperty::PolarInertia::
operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, offy, offz, A;
                _DIM1(p, t, DIM1);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                _Afunc(DIM1, A);
                
                _func(DIM1, m);
                m += A*pow(offy,2) + A*pow(offz,2); // Account for offsets
            }
    

void MAST::Solid1D1ParameterSectionProperty::PolarInertia::
derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, dDIM1, offy, offz, doffy, doffz, A, dA;
                _DIM1        (p, t, DIM1);           _DIM1.derivative( f, p, t, dDIM1);
                _hy_offset (p, t, offy);  _hy_offset.derivative( f, p, t, doffy);
                _hz_offset (p, t, offz);  _hz_offset.derivative( f, p, t, doffz);
                
                
                _Afunc(DIM1, A);
                _dAfunc(DIM1, dDIM1, dA);
                                
                _dfunc(DIM1, dDIM1, m);
                m += dA*pow(offy,2) + 2.*A*offy*doffy + dA*pow(offz,2) + 2.*A*offz*doffz;
            }
    
    
MAST::Solid1D1ParameterSectionProperty::AreaYMoment::
AreaYMoment(const std::function<void(Real&, Real&)> Afunc,
                        const std::function<void(Real&, Real&, Real&)> dAfunc,
                        const MAST::FieldFunction<Real>& DIM1,
                        const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _Afunc(Afunc),
            _dAfunc(dAfunc),
            _DIM1(DIM1),
            _hz_offset(hz_offset) {
                _functions.insert(&DIM1);
                _functions.insert(&hz_offset);
            }

void MAST::Solid1D1ParameterSectionProperty::AreaYMoment::
operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, off, A;
                _DIM1(p, t, DIM1);
                _hz_offset(p, t, off);
                
                _Afunc(DIM1, A);
                m = A*off;
            }

void MAST::Solid1D1ParameterSectionProperty::AreaYMoment::
derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, off, dDIM1, doff, A, dA;
                _DIM1        (p, t, DIM1);         _DIM1.derivative( f, p, t, dDIM1);
                _hz_offset (p, t, off); _hz_offset.derivative( f, p, t, doff);
                
                _Afunc(DIM1, A);
                _dAfunc(DIM1, dDIM1, dA);
                
                m = dA*off + A*doff;
            }
    

MAST::Solid1D1ParameterSectionProperty::AreaZMoment::
AreaZMoment(const std::function<void(Real&, Real&)> Afunc,
            const std::function<void(Real&, Real&, Real&)> dAfunc,
            const MAST::FieldFunction<Real>& DIM1,
            const MAST::FieldFunction<Real>&  hy_offset):
MAST::FieldFunction<Real>("AreaZMoment"),
_Afunc(Afunc),
_dAfunc(dAfunc),
_DIM1(DIM1),
_hy_offset(hy_offset) {
    _functions.insert(&DIM1);
    _functions.insert(&hy_offset);
}
    
void MAST::Solid1D1ParameterSectionProperty::AreaZMoment::
operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, off, A;
                _DIM1(p, t, DIM1);
                _hy_offset(p, t, off);
                
                _Afunc(DIM1, A);
                m = A*off;
            }
    
void MAST::Solid1D1ParameterSectionProperty::AreaZMoment::
derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real DIM1, off, dDIM1, doff, A, dA;
                _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
                _hy_offset(p, t, off); _hy_offset.derivative( f, p, t, doff);
                
                _Afunc(DIM1, A);
                _dAfunc(DIM1, dDIM1, dA);
                
                m = dA*off + A*doff;
            }

MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix::
AreaInertiaMatrix(const std::function<void(Real&, Real&)> Izfunc,
                          const std::function<void(Real&, Real&, Real&)> dIzfunc,
                          const std::function<void(Real&, Real&)> Iyfunc,
                          const std::function<void(Real&, Real&, Real&)> dIyfunc,
                          const std::function<void(Real&, Real&)> Afunc,
                          const std::function<void(Real&, Real&, Real&)> dAfunc,
                          const MAST::FieldFunction<Real>& DIM1,
                          const MAST::FieldFunction<Real>&  hy_offset,
                          const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
            _Izfunc(Izfunc),
            _dIzfunc(dIzfunc),
            _Iyfunc(Iyfunc),
            _dIyfunc(dIyfunc),
            _Afunc(Afunc),
            _dAfunc(dAfunc),
            _DIM1(DIM1),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&DIM1);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }

void MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix::
operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real DIM1, offy, offz, A;
                m = RealMatrixX::Zero(2,2);
                _DIM1(p, t, DIM1);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                _Afunc(DIM1, A);
                _Izfunc(DIM1, m(0,0));
                _Iyfunc(DIM1, m(1,1));
                
                m(0,0) += A*pow(offy,2) ; // Account for offset
                m(0,1) = A*offy*offz;
                m(1,0) = m(0,1);
                m(1,1) += A*pow(offz,2) ; // Account for offset
            }

void MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix::
derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real DIM1, offy, offz, A, dDIM1, doffy, doffz, dA;
                m = RealMatrixX::Zero(2,2);
                _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
                _hy_offset(p, t, offy); _hy_offset.derivative( f, p, t, doffy);
                _hz_offset(p, t, offz); _hz_offset.derivative( f, p, t, doffz);
                
                _Afunc(DIM1, A);
                _dAfunc(DIM1, dDIM1, dA);
                _dIzfunc(DIM1, dDIM1, m(0,0));
                _dIyfunc(DIM1, dDIM1, m(1,1));
                
                m(0,0) += dA*pow(offy,2) + A*2.*offy*doffy;
                m(0,1) = dA*offy*offz + A*doffy*offz + A*offy*doffz;
                m(1,0) = m(0,1);
                m(1,1) += dA*pow(offz,2) + A*2.*offz*doffz;
            }
