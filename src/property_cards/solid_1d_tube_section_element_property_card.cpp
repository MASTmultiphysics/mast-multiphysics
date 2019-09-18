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
#include "property_cards/solid_1d_tube_section_element_property_card.h"


#define PI 3.1415926535897932


// Annulus (TUBE in Siemens NX Nastran and Astros 21.2)
void MAST::Solid1DTubeSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = PI*(DIM1*DIM1-DIM2*DIM2);
}

void MAST::Solid1DTubeSectionProperty::calcdA(Real& DIM1, Real& DIM2, 
                                              Real& dDIM1, Real& dDIM2,
                                              Real& dA){
    dA = PI*(DIM1*dDIM1-DIM2*dDIM2)*2.0;
}

void MAST::Solid1DTubeSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/4.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIz(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2, 
                                               Real& dIz){
    dIz = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2);
}

void MAST::Solid1DTubeSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/4.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIy(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2, 
                                               Real& dIy){
    dIy = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2);
}

void MAST::Solid1DTubeSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/2.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIp(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2,
                                               Real& dIp){
    dIp = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2)*2.0;
}

void MAST::Solid1DTubeSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    J = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/2.0;
}

void MAST::Solid1DTubeSectionProperty::calcdJ(Real& DIM1, Real& DIM2, 
                                              Real& dDIM1, Real& dDIM2, 
                                              Real& dJ){
    dJ = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2)*2.0;
}


void MAST::Solid1DTubeSectionElementPropertyCard::init() {
    
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
    else if (DIM2v>=DIM1v){
        libmesh_error_msg("DIM2>=DIM1");
    }
    
    _A.reset(new MAST::Solid1D2ParameterSectionProperty::Area(MAST::Solid1DTubeSectionProperty::calcA,
                                                              MAST::Solid1DTubeSectionProperty::calcdA,
                                                              DIM1, DIM2));
    
    _Ay.reset(new MAST::Solid1D2ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D2ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D2ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DTubeSectionProperty::calcJ,
                                                                MAST::Solid1DTubeSectionProperty::calcdJ,
                                                                DIM1, DIM2));
    
    _Ip.reset(new MAST::Solid1D2ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DTubeSectionProperty::calcIp,
                                                                MAST::Solid1DTubeSectionProperty::calcdIp,
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D2ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DTubeSectionProperty::calcIz,
                                                                MAST::Solid1DTubeSectionProperty::calcdIz,
                                                                MAST::Solid1DTubeSectionProperty::calcIy,
                                                                MAST::Solid1DTubeSectionProperty::calcdIy,
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
