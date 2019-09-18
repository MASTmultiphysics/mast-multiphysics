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
#include "property_cards/solid_1d_hollowellipse_section_element_property_card.h"

#define PI 3.141592653589793

// Hollow Ellipse
void MAST::Solid1DHollowEllipseSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A){
    A = (PI*(DIM1*DIM2-DIM3*DIM4))/4.0;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA){
    dA = (PI*(DIM1*dDIM2+DIM2*dDIM1-DIM3*dDIM4-DIM4*dDIM3))/4.0;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz){
    Iz = (PI*(DIM1*(DIM2*DIM2*DIM2)-DIM3*(DIM4*DIM4*DIM4)))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz){
    dIz = (PI*((DIM2*DIM2*DIM2)*dDIM1-(DIM4*DIM4*DIM4)*dDIM3+DIM1*(DIM2*DIM2)*dDIM2*3.0-DIM3*(DIM4*DIM4)*dDIM4*3.0))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy){
    Iy = (PI*((DIM1*DIM1*DIM1)*DIM2-(DIM3*DIM3*DIM3)*DIM4))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy){
    dIy = (PI*((DIM1*DIM1*DIM1)*dDIM2-(DIM3*DIM3*DIM3)*dDIM4+(DIM1*DIM1)*DIM2*dDIM1*3.0-(DIM3*DIM3)*DIM4*dDIM3*3.0))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip){
    Ip = (PI*(DIM1*(DIM2*DIM2*DIM2)-DIM3*(DIM4*DIM4*DIM4)))/6.4E+1+(PI*((DIM1*DIM1*DIM1)*DIM2-(DIM3*DIM3*DIM3)*DIM4))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp){
    dIp = (DIM1*dDIM2*PI*(DIM1*DIM1+(DIM2*DIM2)*3.0))/6.4E+1+(DIM2*dDIM1*PI*((DIM1*DIM1)*3.0+DIM2*DIM2))/6.4E+1-(DIM3*dDIM4*PI*(DIM3*DIM3+(DIM4*DIM4)*3.0))/6.4E+1-(DIM4*dDIM3*PI*((DIM3*DIM3)*3.0+DIM4*DIM4))/6.4E+1;
}

void MAST::Solid1DHollowEllipseSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J){
    J = ((DIM1*DIM1*DIM1)*PI*(DIM2*DIM2*DIM2*DIM2-DIM4*DIM4*DIM4*DIM4))/(DIM2*(DIM1*DIM1+DIM2*DIM2)*1.6E+1);
}

void MAST::Solid1DHollowEllipseSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ){
    dJ = ((DIM1*DIM1*DIM1)*(DIM4*DIM4*DIM4)*dDIM4*PI*(-1.0/4.0))/(DIM2*(DIM1*DIM1+DIM2*DIM2))+((DIM1*DIM1*DIM1)*1.0/(DIM2*DIM2)*dDIM2*PI*1.0/pow(DIM1*DIM1+DIM2*DIM2,2.0)*(DIM2*DIM2*DIM2*DIM2*DIM2*DIM2+(DIM1*DIM1)*(DIM2*DIM2*DIM2*DIM2)*3.0+(DIM1*DIM1)*(DIM4*DIM4*DIM4*DIM4)+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)*3.0))/1.6E+1+((DIM1*DIM1)*dDIM1*PI*(DIM1*DIM1+(DIM2*DIM2)*3.0)*(DIM2*DIM2*DIM2*DIM2-DIM4*DIM4*DIM4*DIM4)*1.0/pow(DIM1*DIM1+DIM2*DIM2,2.0))/(DIM2*1.6E+1);
}



void MAST::Solid1DHollowEllipseSectionElementPropertyCard::init() {
    
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
    else if (DIM3v>=DIM1v){
        libmesh_error_msg("DIM3>=DIM1");
    }
    else if (DIM4v>=DIM2v){
        libmesh_error_msg("DIM4>=DIM2");
    }
    
    
    _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DHollowEllipseSectionProperty::calcA,
                                                              MAST::Solid1DHollowEllipseSectionProperty::calcdA,
                                                              DIM1, DIM2, DIM3, 
                                                              DIM4));
    
    _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcA,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcA,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcJ,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdJ,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4));
    
    _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcIp,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdIp,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcA,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcIz,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdIz,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcIy,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdIy,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcA,
                                                                MAST::Solid1DHollowEllipseSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
