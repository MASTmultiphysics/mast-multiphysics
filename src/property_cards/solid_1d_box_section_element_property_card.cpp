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
#include "property_cards/solid_1d_box_section_element_property_card.h"


// Hollow Rectangle (BOX in Siemens NX Nastran and Astros 21.2)
void MAST::Solid1DBoxSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A){
    A = DIM1*DIM2-(DIM1-DIM4*2.0)*(DIM2-DIM3*2.0);
}

void MAST::Solid1DBoxSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA){
    dA = DIM3*dDIM1*2.0+DIM4*dDIM2*2.0+dDIM3*(DIM1*2.0-DIM4*4.0)+dDIM4*(DIM2*2.0-DIM3*4.0);
}

void MAST::Solid1DBoxSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz){
    Iz = (DIM1*(DIM2*DIM2*DIM2))/1.2E+1-((DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,3.0))/1.2E+1;
}

void MAST::Solid1DBoxSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz){
    dIz = dDIM1*((DIM2*DIM2*DIM2)/1.2E+1-pow(DIM2-DIM3*2.0,3.0)/1.2E+1)+(dDIM4*pow(DIM2-DIM3*2.0,3.0))/6.0+dDIM2*((DIM1*(DIM2*DIM2))/4.0-((DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,2.0))/4.0)+(dDIM3*(DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,2.0))/2.0;
}

void MAST::Solid1DBoxSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM2)/1.2E+1-(pow(DIM1-DIM4*2.0,3.0)*(DIM2-DIM3*2.0))/1.2E+1;
}

void MAST::Solid1DBoxSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy){
    dIy = dDIM2*((DIM1*DIM1*DIM1)/1.2E+1-pow(DIM1-DIM4*2.0,3.0)/1.2E+1)+(dDIM3*pow(DIM1-DIM4*2.0,3.0))/6.0+dDIM1*(((DIM1*DIM1)*DIM2)/4.0-(pow(DIM1-DIM4*2.0,2.0)*(DIM2-DIM3*2.0))/4.0)+(dDIM4*pow(DIM1-DIM4*2.0,2.0)*(DIM2-DIM3*2.0))/2.0;
}

void MAST::Solid1DBoxSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip){
    Ip = (DIM1*(DIM2*DIM2*DIM2))/1.2E+1+((DIM1*DIM1*DIM1)*DIM2)/1.2E+1-((DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,3.0))/1.2E+1-(pow(DIM1-DIM4*2.0,3.0)*(DIM2-DIM3*2.0))/1.2E+1;
}

void MAST::Solid1DBoxSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp){
    dIp = dDIM3*(((DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,2.0))/2.0+pow(DIM1-DIM4*2.0,3.0)/6.0)+dDIM4*((pow(DIM1-DIM4*2.0,2.0)*(DIM2-DIM3*2.0))/2.0+pow(DIM2-DIM3*2.0,3.0)/6.0)+dDIM1*(((DIM1*DIM1)*DIM2)/4.0-(pow(DIM1-DIM4*2.0,2.0)*(DIM2-DIM3*2.0))/4.0+(DIM2*DIM2*DIM2)/1.2E+1-pow(DIM2-DIM3*2.0,3.0)/1.2E+1)+dDIM2*((DIM1*(DIM2*DIM2))/4.0-((DIM1-DIM4*2.0)*pow(DIM2-DIM3*2.0,2.0))/4.0+(DIM1*DIM1*DIM1)/1.2E+1-pow(DIM1-DIM4*2.0,3.0)/1.2E+1);
}

void MAST::Solid1DBoxSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J){
    J = (DIM3*DIM4*pow(DIM1-DIM4,2.0)*pow(DIM2-DIM3,2.0)*-2.0)/(DIM3*DIM3+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3);
}

void MAST::Solid1DBoxSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ){
    dJ = DIM3*dDIM4*(DIM1-DIM4)*pow(DIM2-DIM3,2.0)*1.0/pow(DIM3*DIM3+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3,2.0)*(DIM1*(DIM3*DIM3)+DIM1*(DIM4*DIM4)-(DIM3*DIM3)*DIM4*3.0-DIM4*DIM4*DIM4-DIM1*DIM2*DIM3+DIM2*DIM3*DIM4*3.0)*-2.0-DIM4*dDIM3*pow(DIM1-DIM4,2.0)*(DIM2-DIM3)*1.0/pow(DIM3*DIM3+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3,2.0)*(DIM2*(DIM3*DIM3)+DIM2*(DIM4*DIM4)-DIM3*(DIM4*DIM4)*3.0-DIM3*DIM3*DIM3-DIM1*DIM2*DIM4+DIM1*DIM3*DIM4*3.0)*2.0-DIM3*DIM4*dDIM1*(DIM1-DIM4)*pow(DIM2-DIM3,2.0)*1.0/pow(DIM3*DIM3+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3,2.0)*((DIM3*DIM3)*2.0+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3*2.0)*2.0-DIM3*DIM4*dDIM2*pow(DIM1-DIM4,2.0)*(DIM2-DIM3)*1.0/pow(DIM3*DIM3+DIM4*DIM4-DIM1*DIM4-DIM2*DIM3,2.0)*(DIM3*DIM3+(DIM4*DIM4)*2.0-DIM1*DIM4*2.0-DIM2*DIM3)*2.0;
}


void MAST::Solid1DBoxSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v, DIM3v, DIM4v, DIM5v;
    DIM1(DIM1v); DIM2(DIM2v); DIM3(DIM3v); DIM4(DIM4v);
    if (DIM1v<=0){
        libmesh_error_msg("DIM1<=0");
    }
    else if (DIM2v<=0){
        libmesh_error_msg("DIM2<=0");
    }
    else if (DIM3v<=0){
        libmesh_error_msg("DIM3<=0");
    }
    else if (DIM4v<=0){
        libmesh_error_msg("DIM4<=0");
    }
    else if (DIM4v>=(0.5*DIM1v)){
        libmesh_error_msg("DIM4>=0.5*DIM1");
    }
    else if (DIM3v>=(0.5*DIM2v)){
        libmesh_error_msg("DIM3>=0.5*DIM2");
    }
    
    
    _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DBoxSectionProperty::calcA,
                                                              MAST::Solid1DBoxSectionProperty::calcdA,
                                                              DIM1, DIM2, DIM3, 
                                                              DIM4));
    
    _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DBoxSectionProperty::calcA,
                                                                MAST::Solid1DBoxSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DBoxSectionProperty::calcA,
                                                                MAST::Solid1DBoxSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DBoxSectionProperty::calcJ,
                                                                MAST::Solid1DBoxSectionProperty::calcdJ,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4));
    
    _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DBoxSectionProperty::calcIp,
                                                                MAST::Solid1DBoxSectionProperty::calcdIp,
                                                                MAST::Solid1DBoxSectionProperty::calcA,
                                                                MAST::Solid1DBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DBoxSectionProperty::calcIz,
                                                                MAST::Solid1DBoxSectionProperty::calcdIz,
                                                                MAST::Solid1DBoxSectionProperty::calcIy,
                                                                MAST::Solid1DBoxSectionProperty::calcdIy,
                                                                MAST::Solid1DBoxSectionProperty::calcA,
                                                                MAST::Solid1DBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
