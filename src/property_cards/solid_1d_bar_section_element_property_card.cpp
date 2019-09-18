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
#include "property_cards/solid_1d_2parameter_section_element_property_card.h"
#include "property_cards/solid_1d_bar_section_element_property_card.h"


// Rectangle (BAR in Siemens NX Nastran and Astros 21.2)
void MAST::Solid1DBarSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = DIM1*DIM2;
}

void MAST::Solid1DBarSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dA){
    dA = DIM1*dDIM2+DIM2*dDIM1;
}

void MAST::Solid1DBarSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (DIM1*(DIM2*DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIz){
    dIz = ((DIM2*DIM2*DIM2)*dDIM1)/1.2E+1+(DIM1*(DIM2*DIM2)*dDIM2)/4.0;
}

void MAST::Solid1DBarSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM2)/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIy){
    dIy = ((DIM1*DIM1*DIM1)*dDIM2)/1.2E+1+((DIM1*DIM1)*DIM2*dDIM1)/4.0;
}

void MAST::Solid1DBarSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (DIM1*DIM2*(DIM1*DIM1+DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIp){
    dIp = (DIM1*dDIM2*(DIM1*DIM1+(DIM2*DIM2)*3.0))/1.2E+1+(DIM2*dDIM1*((DIM1*DIM1)*3.0+DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcJ1(Real& DIM1, Real& DIM2, Real& J1){
    J1 = (1.0/(DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2*DIM2)*((DIM1*DIM1*DIM1*DIM1)*DIM2*-2.52E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*4.0E+2+(DIM2*DIM2*DIM2*DIM2*DIM2)*2.1E+1))/1.2E+3;
}

void MAST::Solid1DBarSectionProperty::calcdJ1(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ1){
    dJ1 = (1.0/(DIM1*DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2*DIM2)*dDIM1*((DIM1*DIM1*DIM1*DIM1*DIM1)*1.0E+2-(DIM2*DIM2*DIM2*DIM2*DIM2)*2.1E+1))/3.0E+2+(1.0/(DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2)*dDIM2*((DIM1*DIM1*DIM1*DIM1)*DIM2*-4.2E+1+(DIM1*DIM1*DIM1*DIM1*DIM1)*5.0E+1+(DIM2*DIM2*DIM2*DIM2*DIM2)*7.0))/5.0E+1;
}

void MAST::Solid1DBarSectionProperty::calcJ2(Real& DIM1, Real& DIM2, Real& J2){
    J2 = ((DIM1*DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2)*(DIM1*(DIM2*DIM2*DIM2*DIM2)*-2.52E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*2.1E+1+(DIM2*DIM2*DIM2*DIM2*DIM2)*4.0E+2))/1.2E+3;
}

void MAST::Solid1DBarSectionProperty::calcdJ2(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ2){
    dJ2 = (DIM1*DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2*DIM2)*dDIM2*((DIM1*DIM1*DIM1*DIM1*DIM1)*2.1E+1-(DIM2*DIM2*DIM2*DIM2*DIM2)*1.0E+2)*(-1.0/3.0E+2)+((DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2)*dDIM1*(DIM1*(DIM2*DIM2*DIM2*DIM2)*-4.2E+1+(DIM1*DIM1*DIM1*DIM1*DIM1)*7.0+(DIM2*DIM2*DIM2*DIM2*DIM2)*5.0E+1))/5.0E+1;
}


void MAST::Solid1DBarSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    if (DIM1>DIM2)
    {
        MAST::Solid1DBarSectionProperty::calcJ1(DIM1, DIM2, J);
    }
    else
    {
        MAST::Solid1DBarSectionProperty::calcJ2(DIM1, DIM2, J);
    }
}


void MAST::Solid1DBarSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ){
    if (DIM1>DIM2)
    {
        MAST::Solid1DBarSectionProperty::calcdJ1(DIM1, DIM2, dDIM1, dDIM2, dJ);
    }
    else
    {
        MAST::Solid1DBarSectionProperty::calcdJ2(DIM1, DIM2, dDIM1, dDIM2, dJ);
    }
}

void MAST::Solid1DBarSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v;
    DIM1(DIM1v); DIM2(DIM2v);
    if (DIM1v<=0){
        libmesh_error_msg("DIM1<=0");
    }
    else if (DIM2v<=0){
        libmesh_error_msg("DIM2<=0");
    }
    
    _A.reset(new MAST::Solid1D2ParameterSectionProperty::Area(MAST::Solid1DBarSectionProperty::calcA,
                                                              MAST::Solid1DBarSectionProperty::calcdA,
                                                              DIM1, DIM2));
    
    _Ay.reset(new MAST::Solid1D2ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DBarSectionProperty::calcA,
                                                                MAST::Solid1DBarSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D2ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DBarSectionProperty::calcA,
                                                                MAST::Solid1DBarSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D2ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DBarSectionProperty::calcJ,
                                                                MAST::Solid1DBarSectionProperty::calcdJ,
                                                                DIM1, DIM2));
    
    _Ip.reset(new MAST::Solid1D2ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DBarSectionProperty::calcIp,
                                                                MAST::Solid1DBarSectionProperty::calcdIp,
                                                                MAST::Solid1DBarSectionProperty::calcA,
                                                                MAST::Solid1DBarSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D2ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DBarSectionProperty::calcIz,
                                                                MAST::Solid1DBarSectionProperty::calcdIz,
                                                                MAST::Solid1DBarSectionProperty::calcIy,
                                                                MAST::Solid1DBarSectionProperty::calcdIy,
                                                                MAST::Solid1DBarSectionProperty::calcA,
                                                                MAST::Solid1DBarSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
