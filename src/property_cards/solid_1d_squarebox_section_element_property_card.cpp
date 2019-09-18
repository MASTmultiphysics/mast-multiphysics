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
#include "property_cards/solid_1d_squarebox_section_element_property_card.h"


// SquareBox
void MAST::Solid1DSquareBoxSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = DIM2*(DIM1-DIM2)*4.0;
}

void MAST::Solid1DSquareBoxSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dA){
    dA = DIM2*dDIM1*4.0+dDIM2*(DIM1*4.0-DIM2*8.0);
}

void MAST::Solid1DSquareBoxSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (DIM1*DIM1*DIM1*DIM1)/1.2E+1-pow(DIM1-DIM2*2.0,4.0)/1.2E+1;
}

void MAST::Solid1DSquareBoxSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIz){
    dIz = dDIM1*((DIM1*DIM1*DIM1)/3.0-pow(DIM1-DIM2*2.0,3.0)/3.0)+dDIM2*pow(DIM1-DIM2*2.0,3.0)*(2.0/3.0);
}

void MAST::Solid1DSquareBoxSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = (DIM1*DIM1*DIM1*DIM1)/1.2E+1-pow(DIM1-DIM2*2.0,4.0)/1.2E+1;
}

void MAST::Solid1DSquareBoxSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIy){
    dIy = dDIM1*((DIM1*DIM1*DIM1)/3.0-pow(DIM1-DIM2*2.0,3.0)/3.0)+dDIM2*pow(DIM1-DIM2*2.0,3.0)*(2.0/3.0);
}

void MAST::Solid1DSquareBoxSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (DIM1*DIM1*DIM1*DIM1)/6.0-pow(DIM1-DIM2*2.0,4.0)/6.0;
}

void MAST::Solid1DSquareBoxSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIp){
    dIp = dDIM1*((DIM1*DIM1*DIM1)*(2.0/3.0)-pow(DIM1-DIM2*2.0,3.0)*(2.0/3.0))+dDIM2*pow(DIM1-DIM2*2.0,3.0)*(4.0/3.0);
}

void MAST::Solid1DSquareBoxSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    J = (DIM1*(DIM2*DIM2*DIM2))/2.4E+1-((DIM1*DIM1*DIM1)*DIM2)/2.4E+1+(DIM1*DIM1*DIM1*DIM1)*(7.0/4.8E+1)-(DIM2*DIM2*DIM2*DIM2)*(7.0/4.8E+1);
}

void MAST::Solid1DSquareBoxSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ){
    dJ = dDIM1*((DIM1*DIM1)*DIM2*(-1.0/8.0)+(DIM1*DIM1*DIM1)*(7.0/1.2E+1)+(DIM2*DIM2*DIM2)/2.4E+1)-dDIM2*(DIM1*(DIM2*DIM2)*(-1.0/8.0)+(DIM1*DIM1*DIM1)/2.4E+1+(DIM2*DIM2*DIM2)*(7.0/1.2E+1));
}



void MAST::Solid1DSquareBoxSectionElementPropertyCard::init() {
    
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
    
    _A.reset(new MAST::Solid1D2ParameterSectionProperty::Area(MAST::Solid1DSquareBoxSectionProperty::calcA,
                                                              MAST::Solid1DSquareBoxSectionProperty::calcdA,
                                                              DIM1, DIM2));
    
    _Ay.reset(new MAST::Solid1D2ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DSquareBoxSectionProperty::calcA,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D2ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DSquareBoxSectionProperty::calcA,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D2ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DSquareBoxSectionProperty::calcJ,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdJ,
                                                                DIM1, DIM2));
    
    _Ip.reset(new MAST::Solid1D2ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DSquareBoxSectionProperty::calcIp,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdIp,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcA,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D2ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DSquareBoxSectionProperty::calcIz,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdIz,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcIy,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdIy,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcA,
                                                                MAST::Solid1DSquareBoxSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
