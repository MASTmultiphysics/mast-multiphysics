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
#include "property_cards/solid_1d_square_section_element_property_card.h"


// Square (Speical case of BAR in Siemens NX Nastran and Astros 21.2
void MAST::Solid1DSquareSectionProperty::calcA(Real& DIM1, Real& A){
    A = DIM1*DIM1;
}

void MAST::Solid1DSquareSectionProperty::calcdA(Real& DIM1, Real& dDIM1, Real& dA){
    dA = DIM1*dDIM1*2.0;
}

void MAST::Solid1DSquareSectionProperty::calcIz(Real& DIM1, Real& Iz){
    Iz = (DIM1*DIM1*DIM1*DIM1)/1.2E+1;
}

void MAST::Solid1DSquareSectionProperty::calcdIz(Real& DIM1, Real& dDIM1, Real& dIz){
    dIz = ((DIM1*DIM1*DIM1)*dDIM1)/3.0;
}

void MAST::Solid1DSquareSectionProperty::calcIy(Real& DIM1, Real& Iy){
    Iy = (DIM1*DIM1*DIM1*DIM1)/1.2E+1;
}

void MAST::Solid1DSquareSectionProperty::calcdIy(Real& DIM1, Real& dDIM1, Real& dIy){
    dIy = ((DIM1*DIM1*DIM1)*dDIM1)/3.0;
}

void MAST::Solid1DSquareSectionProperty::calcIp(Real& DIM1, Real& Ip){
    Ip = (DIM1*DIM1*DIM1*DIM1)/6.0;
}

void MAST::Solid1DSquareSectionProperty::calcdIp(Real& DIM1, Real& dDIM1, Real& dIp){
    dIp = (DIM1*DIM1*DIM1)*dDIM1*(2.0/3.0);
}

void MAST::Solid1DSquareSectionProperty::calcJ(Real& DIM1, Real& J){
    J = (DIM1*DIM1*DIM1*DIM1)*(9.0/6.4E+1);
}

void MAST::Solid1DSquareSectionProperty::calcdJ(Real& DIM1, Real& dDIM1, Real& dJ){
    dJ = (DIM1*DIM1*DIM1)*dDIM1*(9.0/1.6E+1);
}


void MAST::Solid1DSquareSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v;
    DIM1(DIM1v);
    if (DIM1v<=0){
        libmesh_error_msg("DIM1<=0");
    }
    
    _A.reset(new MAST::Solid1D1ParameterSectionProperty::Area(MAST::Solid1DSquareSectionProperty::calcA,
                                                              MAST::Solid1DSquareSectionProperty::calcdA,
                                                              DIM1));
    
    _Ay.reset(new MAST::Solid1D1ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DSquareSectionProperty::calcA,
                                                                MAST::Solid1DSquareSectionProperty::calcdA,
                                                                DIM1, 
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D1ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DSquareSectionProperty::calcA,
                                                                MAST::Solid1DSquareSectionProperty::calcdA,
                                                                DIM1, 
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D1ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DSquareSectionProperty::calcJ,
                                                                MAST::Solid1DSquareSectionProperty::calcdJ,
                                                                DIM1));
    
    _Ip.reset(new MAST::Solid1D1ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DSquareSectionProperty::calcIp,
                                                                MAST::Solid1DSquareSectionProperty::calcdIp,
                                                                MAST::Solid1DSquareSectionProperty::calcA,
                                                                MAST::Solid1DSquareSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DSquareSectionProperty::calcIz,
                                                                MAST::Solid1DSquareSectionProperty::calcdIz,
                                                                MAST::Solid1DSquareSectionProperty::calcIy,
                                                                MAST::Solid1DSquareSectionProperty::calcdIy,
                                                                MAST::Solid1DSquareSectionProperty::calcA,
                                                                MAST::Solid1DSquareSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
