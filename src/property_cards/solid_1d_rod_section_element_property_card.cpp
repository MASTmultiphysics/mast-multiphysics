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
#include "property_cards/solid_1d_rod_section_element_property_card.h"


#define PI 3.1415926535897932


void MAST::Solid1DRodSectionProperty::calcA(Real& DIM1, Real& A){
    A = PI*DIM1*DIM1;
}

void MAST::Solid1DRodSectionProperty::calcdA(Real& DIM1, Real& dDIM1, Real& dA){
    dA = 2.*PI*DIM1*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIz(Real& DIM1, Real& Iz){
    Iz = PI*pow(DIM1,4.)/4.0;
}

void MAST::Solid1DRodSectionProperty::calcdIz(Real& DIM1, Real& dDIM1, Real& dIz){
    dIz = PI*pow(DIM1,3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIy(Real& DIM1, Real& Iy){
    Iy = PI*pow(DIM1,4.)/4.0; 
}

void MAST::Solid1DRodSectionProperty::calcdIy(Real& DIM1, Real& dDIM1, Real& dIy){
    dIy = PI*pow(DIM1,3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIp(Real& DIM1, Real& Ip){
    Ip = PI*pow(DIM1, 4.)/2.0;
}

void MAST::Solid1DRodSectionProperty::calcdIp(Real& DIM1, Real& dDIM1, Real& dIp){
    dIp = 2.*PI*pow(DIM1, 3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcJ(Real& DIM1, Real& J){
    J = PI*pow(DIM1, 4.)/2.0;
}

void MAST::Solid1DRodSectionProperty::calcdJ(Real& DIM1, Real& dDIM1, Real& dJ){
    dJ = 2.*PI*pow(DIM1, 3.)*dDIM1;
}


void MAST::Solid1DRodSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v; DIM1(DIM1v);
    if (DIM1v<=0){
        libmesh_error_msg("DIM1<=0");
    }
    
    _A.reset(new MAST::Solid1D1ParameterSectionProperty::Area(MAST::Solid1DRodSectionProperty::calcA,
                                                              MAST::Solid1DRodSectionProperty::calcdA,
                                                              DIM1));
    
    _Ay.reset(new MAST::Solid1D1ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1, 
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D1ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1, 
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D1ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DRodSectionProperty::calcJ,
                                                                MAST::Solid1DRodSectionProperty::calcdJ,
                                                                DIM1));
    
    _Ip.reset(new MAST::Solid1D1ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DRodSectionProperty::calcIp,
                                                                MAST::Solid1DRodSectionProperty::calcdIp,
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DRodSectionProperty::calcIz,
                                                                MAST::Solid1DRodSectionProperty::calcdIz,
                                                                MAST::Solid1DRodSectionProperty::calcIy,
                                                                MAST::Solid1DRodSectionProperty::calcdIy,
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
