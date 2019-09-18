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
#include "property_cards/solid_1d_reghexagon_section_element_property_card.h"


// RegHexagon (Speical case of HEXA in Nastran
void MAST::Solid1DRegHexagonSectionProperty::calcA(Real& DIM1, Real& A){
    A = (sqrt(3.0)*(DIM1*DIM1))/2.0;
}

void MAST::Solid1DRegHexagonSectionProperty::calcdA(Real& DIM1, Real& dDIM1, Real& dA){
    dA = sqrt(3.0)*DIM1*dDIM1;
}

void MAST::Solid1DRegHexagonSectionProperty::calcIz(Real& DIM1, Real& Iz){
    Iz = sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1)*(5.0/1.44E+2);
}

void MAST::Solid1DRegHexagonSectionProperty::calcdIz(Real& DIM1, Real& dDIM1, Real& dIz){
    dIz = sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1*(5.0/3.6E+1);
}

void MAST::Solid1DRegHexagonSectionProperty::calcIy(Real& DIM1, Real& Iy){
    Iy = sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1)*(5.0/1.44E+2);
}

void MAST::Solid1DRegHexagonSectionProperty::calcdIy(Real& DIM1, Real& dDIM1, Real& dIy){
    dIy = sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1*(5.0/3.6E+1);
}

void MAST::Solid1DRegHexagonSectionProperty::calcIp(Real& DIM1, Real& Ip){
    Ip = sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1)*(5.0/7.2E+1);
}

void MAST::Solid1DRegHexagonSectionProperty::calcdIp(Real& DIM1, Real& dDIM1, Real& dIp){
    dIp = sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1*(5.0/1.8E+1);
}

void MAST::Solid1DRegHexagonSectionProperty::calcJ(Real& DIM1, Real& J){
    J = (DIM1*DIM1*DIM1*DIM1)*1.154E-1;
}

void MAST::Solid1DRegHexagonSectionProperty::calcdJ(Real& DIM1, Real& dDIM1, Real& dJ){
    dJ = (DIM1*DIM1*DIM1)*dDIM1*4.616E-1;
}




void MAST::Solid1DRegHexagonSectionElementPropertyCard::init() {
    
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
    
    _A.reset(new MAST::Solid1D1ParameterSectionProperty::Area(MAST::Solid1DRegHexagonSectionProperty::calcA,
                                                              MAST::Solid1DRegHexagonSectionProperty::calcdA,
                                                              DIM1));
    
    _Ay.reset(new MAST::Solid1D1ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DRegHexagonSectionProperty::calcA,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdA,
                                                                DIM1, 
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D1ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DRegHexagonSectionProperty::calcA,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdA,
                                                                DIM1, 
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D1ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DRegHexagonSectionProperty::calcJ,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdJ,
                                                                DIM1));
    
    _Ip.reset(new MAST::Solid1D1ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DRegHexagonSectionProperty::calcIp,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdIp,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcA,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DRegHexagonSectionProperty::calcIz,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdIz,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcIy,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdIy,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcA,
                                                                MAST::Solid1DRegHexagonSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
