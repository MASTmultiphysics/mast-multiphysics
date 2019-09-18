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
#include "property_cards/solid_1d_ellipse_section_element_property_card.h"


#define PI 3.1415926535897932


// Ellipse
void MAST::Solid1DEllipseSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = (DIM1*DIM2*PI)/4.0;
}

void MAST::Solid1DEllipseSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dA){
    dA = (PI*(DIM1*dDIM2+DIM2*dDIM1))/4.0;
}

void MAST::Solid1DEllipseSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (DIM1*(DIM2*DIM2*DIM2)*PI)/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIz){
    dIz = ((DIM2*DIM2)*PI*(DIM1*dDIM2*3.0+DIM2*dDIM1))/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM2*PI)/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIy){
    dIy = ((DIM1*DIM1)*PI*(DIM1*dDIM2+DIM2*dDIM1*3.0))/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (DIM1*DIM2*PI*(DIM1*DIM1+DIM2*DIM2))/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIp){
    dIp = (DIM1*dDIM2*PI*(DIM1*DIM1+(DIM2*DIM2)*3.0))/6.4E+1+(DIM2*dDIM1*PI*((DIM1*DIM1)*3.0+DIM2*DIM2))/6.4E+1;
}

void MAST::Solid1DEllipseSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    J = ((DIM1*DIM1*DIM1)*(DIM2*DIM2*DIM2)*PI)/((DIM1*DIM1)*1.6E+1+(DIM2*DIM2)*1.6E+1);
}

void MAST::Solid1DEllipseSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ){
    dJ = ((DIM1*DIM1)*(DIM2*DIM2)*PI*1.0/pow(DIM1*DIM1+DIM2*DIM2,2.0)*((DIM1*DIM1*DIM1)*dDIM2*3.0+(DIM2*DIM2*DIM2)*dDIM1*3.0+(DIM1*DIM1)*DIM2*dDIM1+DIM1*(DIM2*DIM2)*dDIM2))/1.6E+1;
}


void MAST::Solid1DEllipseSectionElementPropertyCard::init() {
    
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
    
    _A.reset(new MAST::Solid1D2ParameterSectionProperty::Area(MAST::Solid1DEllipseSectionProperty::calcA,
                                                              MAST::Solid1DEllipseSectionProperty::calcdA,
                                                              DIM1, DIM2));
    
    _Ay.reset(new MAST::Solid1D2ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DEllipseSectionProperty::calcA,
                                                                MAST::Solid1DEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D2ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DEllipseSectionProperty::calcA,
                                                                MAST::Solid1DEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D2ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DEllipseSectionProperty::calcJ,
                                                                MAST::Solid1DEllipseSectionProperty::calcdJ,
                                                                DIM1, DIM2));
    
    _Ip.reset(new MAST::Solid1D2ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DEllipseSectionProperty::calcIp,
                                                                MAST::Solid1DEllipseSectionProperty::calcdIp,
                                                                MAST::Solid1DEllipseSectionProperty::calcA,
                                                                MAST::Solid1DEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D2ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DEllipseSectionProperty::calcIz,
                                                                MAST::Solid1DEllipseSectionProperty::calcdIz,
                                                                MAST::Solid1DEllipseSectionProperty::calcIy,
                                                                MAST::Solid1DEllipseSectionProperty::calcdIy,
                                                                MAST::Solid1DEllipseSectionProperty::calcA,
                                                                MAST::Solid1DEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
