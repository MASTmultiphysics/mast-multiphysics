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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

// MAST includes
#include "property_cards/solid_1d_8parameter_section_element_property_card.h"


MAST::Solid1D8ParameterSectionProperty::Area::
Area(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                              Real&, Real&, Real&, Real&)> func,
    const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                            Real&, Real&, Real&, Real&, Real&, 
                            Real&, Real&, Real&, Real&, Real&,
                            Real&, Real&)> dfunc,
     const MAST::FieldFunction<Real>& DIM1,
     const MAST::FieldFunction<Real>& DIM2,
     const MAST::FieldFunction<Real>& DIM3,
     const MAST::FieldFunction<Real>& DIM4,
     const MAST::FieldFunction<Real>& DIM5,
     const MAST::FieldFunction<Real>& DIM6,
     const MAST::FieldFunction<Real>& DIM7,
     const MAST::FieldFunction<Real>& DIM8):
MAST::FieldFunction<Real>("Area"),
_func(func), _dfunc(dfunc), _DIM1(DIM1), _DIM2(DIM2), _DIM3(DIM3), 
_DIM4(DIM4), _DIM5(DIM5), _DIM6(DIM6), _DIM7(DIM7), _DIM8(DIM8){
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
}
    
void MAST::Solid1D8ParameterSectionProperty::Area::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _func(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, m);
}
    
void MAST::Solid1D8ParameterSectionProperty::Area::
derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
           Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8;
    _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _dfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, 
           dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, m);
}


MAST::Solid1D8ParameterSectionProperty::TorsionalConstant::
TorsionalConstant(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&)> func,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&,
                                           Real&, Real&)> dfunc,
                  const MAST::FieldFunction<Real>& DIM1,
                  const MAST::FieldFunction<Real>& DIM2,
                  const MAST::FieldFunction<Real>& DIM3,
                  const MAST::FieldFunction<Real>& DIM4,
                  const MAST::FieldFunction<Real>& DIM5,
                  const MAST::FieldFunction<Real>& DIM6,
                  const MAST::FieldFunction<Real>& DIM7,
                  const MAST::FieldFunction<Real>& DIM8):
MAST::FieldFunction<Real>("TorsionalConstant"),
_func(func), _dfunc(dfunc), _DIM1(DIM1), _DIM2(DIM2), 
_DIM3(DIM3), _DIM4(DIM4), _DIM5(DIM5), _DIM6(DIM6), 
_DIM7(DIM7), _DIM8(DIM8)
{
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
}
    

void MAST::Solid1D8ParameterSectionProperty::TorsionalConstant::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _func(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, m);
}

    
void MAST::Solid1D8ParameterSectionProperty::TorsionalConstant::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8;
    _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _dfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, dDIM3, 
           dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, m);
}
    

    
MAST::Solid1D8ParameterSectionProperty::PolarInertia::
PolarInertia(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&)> func,
             const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&, Real&,
                                      Real&, Real&)> dfunc,
             const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&)> Afunc,
             const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&, Real&, 
                                      Real&, Real&, Real&, Real&, Real&,
                                      Real&, Real&)> dAfunc,
            const MAST::FieldFunction<Real>& DIM1,
            const MAST::FieldFunction<Real>& DIM2,
            const MAST::FieldFunction<Real>& DIM3,
            const MAST::FieldFunction<Real>& DIM4,
            const MAST::FieldFunction<Real>& DIM5,
            const MAST::FieldFunction<Real>& DIM6,
            const MAST::FieldFunction<Real>& DIM7,
            const MAST::FieldFunction<Real>& DIM8,
            const MAST::FieldFunction<Real>&  hy_offset,
            const MAST::FieldFunction<Real>&  hz_offset):
MAST::FieldFunction<Real>("PolarInertia"),
_func(func),
_dfunc(dfunc),
_Afunc(Afunc),
_dAfunc(dAfunc),
_DIM1(DIM1), _DIM2(DIM2), _DIM3(DIM3), _DIM4(DIM4), _DIM5(DIM5), 
_DIM6(DIM6), _DIM7(DIM7), _DIM8(DIM8),
_hy_offset(hy_offset),
_hz_offset(hz_offset) 
{
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
    _functions.insert(&hy_offset);
    _functions.insert(&hz_offset);
}

void MAST::Solid1D8ParameterSectionProperty::PolarInertia::
operator() (const libMesh::Point& p,const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, offy, offz, A;
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _hy_offset(p, t, offy);
    _hz_offset(p, t, offz);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    
    _func(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, m);
    m += A*pow(offy,2) + A*pow(offz,2); // Account for offsets
}
    

void MAST::Solid1D8ParameterSectionProperty::PolarInertia::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p,
            const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8; 
    Real offy, offz, doffy, doffz, A, dA;
    _DIM1        (p, t, DIM1);           _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _hy_offset (p, t, offy);  _hy_offset.derivative( f, p, t, doffy);
    _hz_offset (p, t, offz);  _hz_offset.derivative( f, p, t, doffz);
    
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    _dAfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
            dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, dA);
                    
    _dfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, dDIM3, 
           dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, m);
    m += dA*pow(offy,2) + 2.*A*offy*doffy + dA*pow(offz,2) + 2.*A*offz*doffz;
}
    
    
MAST::Solid1D8ParameterSectionProperty::AreaYMoment::
AreaYMoment(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&)> Afunc,
            const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&, Real&,
                                     Real&, Real&)> dAfunc,
            const MAST::FieldFunction<Real>& DIM1,
            const MAST::FieldFunction<Real>& DIM2,
            const MAST::FieldFunction<Real>& DIM3,
            const MAST::FieldFunction<Real>& DIM4,
            const MAST::FieldFunction<Real>& DIM5,
            const MAST::FieldFunction<Real>& DIM6,
            const MAST::FieldFunction<Real>& DIM7,
            const MAST::FieldFunction<Real>& DIM8,
            const MAST::FieldFunction<Real>& hz_offset):
MAST::FieldFunction<Real>("AreaYMoment"),
_Afunc(Afunc),
_dAfunc(dAfunc),
_DIM1(DIM1), _DIM2(DIM2), _DIM3(DIM3), _DIM4(DIM4), _DIM5(DIM5), 
_DIM6(DIM6), _DIM7(DIM7), _DIM8(DIM8),
_hz_offset(hz_offset) 
{
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
    _functions.insert(&hz_offset);
}

void MAST::Solid1D8ParameterSectionProperty::AreaYMoment::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, off, A;
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _hz_offset(p, t, off);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    m = A*off;
}

void MAST::Solid1D8ParameterSectionProperty::AreaYMoment::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8; 
    Real off, doff, A, dA;
    _DIM1        (p, t, DIM1);         _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _hz_offset (p, t, off); _hz_offset.derivative( f, p, t, doff);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    _dAfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
            dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, dA);
    
    m = dA*off + A*doff;
}
    


MAST::Solid1D8ParameterSectionProperty::AreaZMoment::
AreaZMoment(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&)> Afunc,
            const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&, Real&, 
                                     Real&, Real&, Real&, Real&, Real&,
                                     Real&, Real&)> dAfunc,
            const MAST::FieldFunction<Real>& DIM1,
            const MAST::FieldFunction<Real>& DIM2,
            const MAST::FieldFunction<Real>& DIM3,
            const MAST::FieldFunction<Real>& DIM4,
            const MAST::FieldFunction<Real>& DIM5,
            const MAST::FieldFunction<Real>& DIM6,
            const MAST::FieldFunction<Real>& DIM7,
            const MAST::FieldFunction<Real>& DIM8,
            const MAST::FieldFunction<Real>& hy_offset):
MAST::FieldFunction<Real>("AreaZMoment"),
_Afunc(Afunc),
_dAfunc(dAfunc),
_DIM1(DIM1), _DIM2(DIM2), _DIM3(DIM3), _DIM4(DIM4), _DIM5(DIM5), _DIM6(DIM6),
_DIM7(DIM7), _DIM8(DIM8),
_hy_offset(hy_offset) 
{
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
    _functions.insert(&hy_offset);
}
    
void MAST::Solid1D8ParameterSectionProperty::AreaZMoment::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, off, A;
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _hy_offset(p, t, off);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    m = A*off;
}
    
void MAST::Solid1D8ParameterSectionProperty::AreaZMoment::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            Real& m) const
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8;
    Real off, doff, A, dA;
    _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _hy_offset(p, t, off); _hy_offset.derivative( f, p, t, doff);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    _dAfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
            dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, dA);
    
    m = dA*off + A*doff;
}

MAST::Solid1D8ParameterSectionProperty::AreaInertiaMatrix::
AreaInertiaMatrix(const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&)> Izfunc,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&,
                                           Real&, Real&)> dIzfunc,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&)> Iyfunc,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&,
                                           Real&, Real&)> dIyfunc,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&)> Afunc,
                  const std::function<void(Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&, 
                                           Real&, Real&, Real&, Real&, Real&,
                                           Real&, Real&)> dAfunc,
                  const MAST::FieldFunction<Real>& DIM1,
                  const MAST::FieldFunction<Real>& DIM2,
                  const MAST::FieldFunction<Real>& DIM3,
                  const MAST::FieldFunction<Real>& DIM4,
                  const MAST::FieldFunction<Real>& DIM5,
                  const MAST::FieldFunction<Real>& DIM6,
                  const MAST::FieldFunction<Real>& DIM7,
                  const MAST::FieldFunction<Real>& DIM8,
                  const MAST::FieldFunction<Real>&  hy_offset,
                  const MAST::FieldFunction<Real>&  hz_offset):
MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
_Izfunc(Izfunc),
_dIzfunc(dIzfunc),
_Iyfunc(Iyfunc),
_dIyfunc(dIyfunc),
_Afunc(Afunc),
_dAfunc(dAfunc),
_DIM1(DIM1), _DIM2(DIM2), _DIM3(DIM3), _DIM4(DIM4), _DIM5(DIM5),
_DIM6(DIM6), _DIM7(DIM7), _DIM8(DIM8),
_hy_offset(hy_offset),
_hz_offset(hz_offset) 
{
    _functions.insert(&DIM1);
    _functions.insert(&DIM2);
    _functions.insert(&DIM3);
    _functions.insert(&DIM4);
    _functions.insert(&DIM5);
    _functions.insert(&DIM6);
    _functions.insert(&DIM7);
    _functions.insert(&DIM8);
    _functions.insert(&hy_offset);
    _functions.insert(&hz_offset);
}

void MAST::Solid1D8ParameterSectionProperty::AreaInertiaMatrix::
operator() (const libMesh::Point& p, const Real t, RealMatrixX& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, offy, offz, A;
    m = RealMatrixX::Zero(2,2);
    _DIM1(p, t, DIM1);
    _DIM2(p, t, DIM2);
    _DIM3(p, t, DIM3);
    _DIM4(p, t, DIM4);
    _DIM5(p, t, DIM5);
    _DIM6(p, t, DIM6);
    _DIM7(p, t, DIM7);
    _DIM8(p, t, DIM8);
    _hy_offset(p, t, offy);
    _hz_offset(p, t, offz);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    _Izfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, m(0,0));
    _Iyfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, m(1,1));
    
    m(0,0) += A*pow(offy,2) ; // Account for offset
    m(0,1) = A*offy*offz;
    m(1,0) = m(0,1);
    m(1,1) += A*pow(offz,2) ; // Account for offset
}

void MAST::Solid1D8ParameterSectionProperty::AreaInertiaMatrix::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            RealMatrixX& m) const 
{
    Real DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8;
    Real dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8;
    Real offy, offz, A, doffy, doffz, dA;
    m = RealMatrixX::Zero(2,2);
    _DIM1(p, t, DIM1); _DIM1.derivative( f, p, t, dDIM1);
    _DIM2(p, t, DIM2); _DIM2.derivative( f, p, t, dDIM2);
    _DIM3(p, t, DIM3); _DIM3.derivative( f, p, t, dDIM3);
    _DIM4(p, t, DIM4); _DIM4.derivative( f, p, t, dDIM4);
    _DIM5(p, t, DIM5); _DIM5.derivative( f, p, t, dDIM5);
    _DIM6(p, t, DIM6); _DIM6.derivative( f, p, t, dDIM6);
    _DIM7(p, t, DIM7); _DIM7.derivative( f, p, t, dDIM7);
    _DIM8(p, t, DIM8); _DIM8.derivative( f, p, t, dDIM8);
    _hy_offset(p, t, offy); _hy_offset.derivative( f, p, t, doffy);
    _hz_offset(p, t, offz); _hz_offset.derivative( f, p, t, doffz);
    
    _Afunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, A);
    _dAfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
            dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, dA);
    _dIzfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
             dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, m(0,0));
    _dIyfunc(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, dDIM1, dDIM2, 
             dDIM3, dDIM4, dDIM5, dDIM6, dDIM7, dDIM8, m(1,1));
    
    m(0,0) += dA*pow(offy,2) + A*2.*offy*doffy;
    m(0,1) = dA*offy*offz + A*doffy*offz + A*offy*doffz;
    m(1,0) = m(0,1);
    m(1,1) += dA*pow(offz,2) + A*2.*offz*doffz;
}
