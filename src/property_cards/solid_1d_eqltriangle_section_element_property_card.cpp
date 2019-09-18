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
#include "property_cards/solid_1d_eqltriangle_section_element_property_card.h"


// EqlTriangle
void MAST::Solid1DEqlTriangleSectionProperty::calcA(Real& DIM1, Real& A){
    A = (sqrt(3.0)*(DIM1*DIM1))/4.0;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcdA(Real& DIM1, Real& dDIM1, Real& dA){
    dA = (sqrt(3.0)*DIM1*dDIM1)/2.0;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcIz(Real& DIM1, Real& Iz){
    Iz = (sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1))/9.6E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcdIz(Real& DIM1, Real& dDIM1, Real& dIz){
    dIz = (sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1)/2.4E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcIy(Real& DIM1, Real& Iy){
    Iy = (sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1))/9.6E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcdIy(Real& DIM1, Real& dDIM1, Real& dIy){
    dIy = (sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1)/2.4E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcIp(Real& DIM1, Real& Ip){
    Ip = (sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1))/4.8E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcdIp(Real& DIM1, Real& dDIM1, Real& dIp){
    dIp = (sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1)/1.2E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcJ(Real& DIM1, Real& J){
    J = (sqrt(3.0)*(DIM1*DIM1*DIM1*DIM1))/8.0E+1;
}

void MAST::Solid1DEqlTriangleSectionProperty::calcdJ(Real& DIM1, Real& dDIM1, Real& dJ){
    dJ = (sqrt(3.0)*(DIM1*DIM1*DIM1)*dDIM1)/2.0E+1;
}



void MAST::Solid1DEqlTriangleSectionElementPropertyCard::init() {
    
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
    
    _A.reset(new MAST::Solid1D1ParameterSectionProperty::Area(MAST::Solid1DEqlTriangleSectionProperty::calcA,
                                                              MAST::Solid1DEqlTriangleSectionProperty::calcdA,
                                                              DIM1));
    
    _Ay.reset(new MAST::Solid1D1ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcA,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdA,
                                                                DIM1, 
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D1ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcA,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdA,
                                                                DIM1, 
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D1ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcJ,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdJ,
                                                                DIM1));
    
    _Ip.reset(new MAST::Solid1D1ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcIp,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdIp,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcA,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcIz,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdIz,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcIy,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdIy,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcA,
                                                                MAST::Solid1DEqlTriangleSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
