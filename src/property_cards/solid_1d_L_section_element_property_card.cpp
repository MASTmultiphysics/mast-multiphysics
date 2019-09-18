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
#include "property_cards/solid_1d_L_section_element_property_card.h"


// Unequal Angle (T in Siemens NX Nastran)
void MAST::Solid1DLSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A){
    A = DIM1*DIM3+DIM4*(DIM2-DIM3);
}

void MAST::Solid1DLSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA){
    dA = DIM3*dDIM1+DIM4*dDIM2+dDIM3*(DIM1-DIM4)+dDIM4*(DIM2-DIM3);
}

void MAST::Solid1DLSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz){
    Iz = (DIM1*(DIM3*DIM3*DIM3))/1.2E+1+(DIM4*pow(DIM2-DIM3,3.0))/1.2E+1+DIM4*(DIM2-DIM3)*pow(DIM2/2.0-(DIM2*DIM4*(DIM2-DIM3))/(DIM1*DIM3*2.0+DIM4*(DIM2-DIM3)*2.0),2.0)+DIM1*(DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.0/pow(DIM1*DIM3*2.0+DIM4*(DIM2-DIM3)*2.0,2.0)*pow(DIM2-DIM3,2.0);
}

void MAST::Solid1DLSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz){
    dIz = (dDIM3*(DIM1-DIM4)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*pow(DIM1*(DIM3*DIM3)-(DIM2*DIM2)*DIM4-(DIM3*DIM3)*DIM4+DIM2*DIM3*DIM4*2.0,2.0))/4.0+(DIM3*dDIM1*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)*3.0+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*4.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*2.0-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*2.0-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*6.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*2.0))/1.2E+1+(DIM4*dDIM2*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*pow(-DIM1*(DIM3*DIM3)+(DIM2*DIM2)*DIM4+(DIM3*DIM3)*DIM4+DIM1*DIM2*DIM3*2.0-DIM2*DIM3*DIM4*2.0,2.0))/4.0+(dDIM4*(DIM2-DIM3)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*4.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*6.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*2.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*2.0-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*4.0-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*4.0-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*6.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*6.0+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*2.0))/1.2E+1;
}

void MAST::Solid1DLSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM3)/1.2E+1+((DIM4*DIM4*DIM4)*(DIM2-DIM3))/1.2E+1+DIM1*DIM3*pow(DIM1*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM3*(DIM1-DIM4))/(DIM1*DIM3*2.0+DIM4*(DIM2-DIM3)*2.0),2.0)+(DIM1*DIM1)*(DIM3*DIM3)*DIM4*1.0/pow(DIM1*DIM3*2.0+DIM4*(DIM2-DIM3)*2.0,2.0)*pow(DIM1-DIM4,2.0)*(DIM2-DIM3);
}

void MAST::Solid1DLSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy){
    dIy = (dDIM4*(DIM2-DIM3)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*pow((DIM1*DIM1)*DIM3-DIM2*(DIM4*DIM4)+DIM3*(DIM4*DIM4)-DIM1*DIM3*DIM4*2.0,2.0))/4.0+(DIM4*dDIM2*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)*3.0+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*4.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*2.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*2.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*6.0+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*2.0))/1.2E+1+(DIM3*dDIM1*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*pow((DIM1*DIM1)*DIM3-DIM2*(DIM4*DIM4)+DIM3*(DIM4*DIM4)+DIM1*DIM2*DIM4*2.0-DIM1*DIM3*DIM4*2.0,2.0))/4.0+(dDIM3*(DIM1-DIM4)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)+(DIM1*DIM1)*(DIM2*DIM2)*(DIM4*DIM4)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*6.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*2.0-DIM1*(DIM2*DIM2)*(DIM4*DIM4*DIM4)*2.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*4.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*4.0-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*6.0+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*6.0+(DIM1*DIM1*DIM1)*DIM2*DIM3*DIM4*2.0))/1.2E+1;
}

void MAST::Solid1DLSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip){
    Ip = ((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*6.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*6.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*2.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*2.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*4.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*4.0-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*4.0-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*4.0-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*6.0-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*6.0+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*4.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*4.0+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*4.0+(DIM1*DIM1*DIM1)*DIM2*DIM3*DIM4*4.0)/(DIM1*DIM3*1.2E+1+DIM2*DIM4*1.2E+1-DIM3*DIM4*1.2E+1);
}

void MAST::Solid1DLSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp){
    dIp = (dDIM4*(DIM2-DIM3)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)*3.0+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)*3.0+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)*3.0+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*1.8E+1+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*6.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*2.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*6.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*2.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*1.2E+1-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*1.2E+1-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*4.0-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*4.0-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*6.0-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*6.0+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*1.2E+1+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*6.0+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*2.0))/1.2E+1+(dDIM3*(DIM1-DIM4)*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)*3.0+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)*3.0+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)*3.0+(DIM1*DIM1)*(DIM2*DIM2)*(DIM4*DIM4)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*6.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*1.8E+1-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*6.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*2.0-DIM1*(DIM2*DIM2)*(DIM4*DIM4*DIM4)*2.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*4.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*4.0-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*1.2E+1-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.2E+1-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*6.0-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*6.0+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*6.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*1.2E+1+(DIM1*DIM1*DIM1)*DIM2*DIM3*DIM4*2.0))/1.2E+1+(DIM3*dDIM1*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)*3.0+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)*3.0+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)*3.0+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)*3.0+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)+(DIM1*DIM1)*(DIM2*DIM2)*(DIM4*DIM4)*1.2E+1+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*1.8E+1+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*4.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*2.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*6.0-DIM1*(DIM2*DIM2)*(DIM4*DIM4*DIM4)*1.2E+1-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*1.2E+1-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*1.2E+1-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*2.0-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*6.0-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*3.0E+1+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*2.4E+1+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*2.0+(DIM1*DIM1*DIM1)*DIM2*DIM3*DIM4*1.2E+1))/1.2E+1+(DIM4*dDIM2*1.0/pow(DIM1*DIM3+DIM2*DIM4-DIM3*DIM4,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)*3.0+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)*3.0+(DIM2*DIM2)*(DIM4*DIM4*DIM4*DIM4)+(DIM2*DIM2*DIM2*DIM2)*(DIM4*DIM4)*3.0+(DIM3*DIM3)*(DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3*DIM3*DIM3)*(DIM4*DIM4)*3.0+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*1.2E+1+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*4.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM4*DIM4)*1.8E+1-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*6.0-DIM2*DIM3*(DIM4*DIM4*DIM4*DIM4)*2.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*1.2E+1-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*2.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*6.0-DIM2*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*1.2E+1-(DIM2*DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.2E+1-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*3.0E+1+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*2.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*2.4E+1+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*1.2E+1))/1.2E+1;
}

void MAST::Solid1DLSectionProperty::calcJ1_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_h){
    J1_h = (DIM3*DIM3*DIM3)*(DIM1-DIM4/2.0)*((DIM3*((DIM3*DIM3*DIM3*DIM3)*1.0/pow(DIM1*2.0-DIM4,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM1*2.0-DIM4)+1.0/3.0);
}

void MAST::Solid1DLSectionProperty::calcdJ1_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_h){
    dJ1_h = dDIM3*((DIM3*DIM3*DIM3)*((DIM3*DIM3*DIM3*DIM3)*1.0/pow(DIM1*2.0-DIM4,5.0)*(5.6E+1/2.5E+1)+(((DIM3*DIM3*DIM3*DIM3)*1.0/pow(DIM1*2.0-DIM4,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM1*2.0-DIM4))*(DIM1-DIM4/2.0)+(DIM3*DIM3)*(DIM1-DIM4/2.0)*((DIM3*((DIM3*DIM3*DIM3*DIM3)*1.0/pow(DIM1*2.0-DIM4,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM1*2.0-DIM4)+1.0/3.0)*3.0)-((DIM3*DIM3*DIM3)*dDIM1*1.0/pow(DIM1*2.0-DIM4,5.0)*(DIM1*(DIM4*DIM4*DIM4*DIM4)*-2.5E+2+(DIM1*DIM1*DIM1*DIM1)*DIM4*2.0E+3-(DIM1*DIM1*DIM1*DIM1*DIM1)*8.0E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*1.68E+2+(DIM4*DIM4*DIM4*DIM4*DIM4)*2.5E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4)*1.0E+3-(DIM1*DIM1*DIM1)*(DIM4*DIM4)*2.0E+3))/7.5E+1+((DIM3*DIM3*DIM3)*dDIM4*1.0/pow(DIM1*2.0-DIM4,5.0)*(DIM1*(DIM4*DIM4*DIM4*DIM4)*-2.5E+2+(DIM1*DIM1*DIM1*DIM1)*DIM4*2.0E+3-(DIM1*DIM1*DIM1*DIM1*DIM1)*8.0E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*1.68E+2+(DIM4*DIM4*DIM4*DIM4*DIM4)*2.5E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4)*1.0E+3-(DIM1*DIM1*DIM1)*(DIM4*DIM4)*2.0E+3))/1.5E+2;
}

void MAST::Solid1DLSectionProperty::calcJ2_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_h){
    J2_h = DIM3*((((1.0/(DIM3*DIM3*DIM3*DIM3)*pow(DIM1-DIM4/2.0,4.0))/1.2E+1-1.0)*(DIM1*(2.1E+1/1.0E+2)-DIM4*(2.1E+1/2.0E+2)))/DIM3+1.0/3.0)*pow(DIM1-DIM4/2.0,3.0);
}

void MAST::Solid1DLSectionProperty::calcdJ2_h(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_h){
    dJ2_h = (1.0/(DIM3*DIM3*DIM3*DIM3)*dDIM1*pow(DIM1*2.0-DIM4,2.0)*(DIM1*(DIM3*DIM3*DIM3*DIM3)*-1.344E+3+DIM1*(DIM4*DIM4*DIM4*DIM4)*7.0E+1-(DIM1*DIM1*DIM1*DIM1)*DIM4*5.6E+2+(DIM3*DIM3*DIM3*DIM3)*DIM4*6.72E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*2.24E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*1.6E+3-(DIM4*DIM4*DIM4*DIM4*DIM4)*7.0-(DIM1*DIM1)*(DIM4*DIM4*DIM4)*2.8E+2+(DIM1*DIM1*DIM1)*(DIM4*DIM4)*5.6E+2))/6.4E+3-(1.0/(DIM3*DIM3*DIM3*DIM3)*dDIM4*pow(DIM1*2.0-DIM4,2.0)*(DIM1*(DIM3*DIM3*DIM3*DIM3)*-1.344E+3+DIM1*(DIM4*DIM4*DIM4*DIM4)*7.0E+1-(DIM1*DIM1*DIM1*DIM1)*DIM4*5.6E+2+(DIM3*DIM3*DIM3*DIM3)*DIM4*6.72E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*2.24E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*1.6E+3-(DIM4*DIM4*DIM4*DIM4*DIM4)*7.0-(DIM1*DIM1)*(DIM4*DIM4*DIM4)*2.8E+2+(DIM1*DIM1*DIM1)*(DIM4*DIM4)*5.6E+2))/1.28E+4+(1.0/(DIM3*DIM3*DIM3*DIM3*DIM3)*dDIM3*pow(DIM1*2.0-DIM4,3.0)*(DIM1*(DIM4*DIM4*DIM4*DIM4)*-2.1E+2+(DIM1*DIM1*DIM1*DIM1)*DIM4*1.68E+3-(DIM1*DIM1*DIM1*DIM1*DIM1)*6.72E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*3.2E+3+(DIM4*DIM4*DIM4*DIM4*DIM4)*2.1E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4)*8.4E+2-(DIM1*DIM1*DIM1)*(DIM4*DIM4)*1.68E+3))/7.68E+4;
}

void MAST::Solid1DLSectionProperty::calcJ1_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_v){
    J1_v = (DIM4*DIM4*DIM4)*(DIM2-DIM3/2.0)*((DIM4*((DIM4*DIM4*DIM4*DIM4)*1.0/pow(DIM2*2.0-DIM3,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM2*2.0-DIM3)+1.0/3.0);
}

void MAST::Solid1DLSectionProperty::calcdJ1_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_v){
    dJ1_v = dDIM4*((DIM4*DIM4*DIM4)*((DIM4*DIM4*DIM4*DIM4)*1.0/pow(DIM2*2.0-DIM3,5.0)*(5.6E+1/2.5E+1)+(((DIM4*DIM4*DIM4*DIM4)*1.0/pow(DIM2*2.0-DIM3,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM2*2.0-DIM3))*(DIM2-DIM3/2.0)+(DIM4*DIM4)*(DIM2-DIM3/2.0)*((DIM4*((DIM4*DIM4*DIM4*DIM4)*1.0/pow(DIM2*2.0-DIM3,4.0)*(4.0/3.0)-1.0)*(2.1E+1/5.0E+1))/(DIM2*2.0-DIM3)+1.0/3.0)*3.0)-((DIM4*DIM4*DIM4)*dDIM2*1.0/pow(DIM2*2.0-DIM3,5.0)*(DIM2*(DIM3*DIM3*DIM3*DIM3)*-2.5E+2+(DIM2*DIM2*DIM2*DIM2)*DIM3*2.0E+3-(DIM2*DIM2*DIM2*DIM2*DIM2)*8.0E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*2.5E+1+(DIM4*DIM4*DIM4*DIM4*DIM4)*1.68E+2+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*1.0E+3-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*2.0E+3))/7.5E+1+((DIM4*DIM4*DIM4)*dDIM3*1.0/pow(DIM2*2.0-DIM3,5.0)*(DIM2*(DIM3*DIM3*DIM3*DIM3)*-2.5E+2+(DIM2*DIM2*DIM2*DIM2)*DIM3*2.0E+3-(DIM2*DIM2*DIM2*DIM2*DIM2)*8.0E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*2.5E+1+(DIM4*DIM4*DIM4*DIM4*DIM4)*1.68E+2+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*1.0E+3-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*2.0E+3))/1.5E+2;
}

void MAST::Solid1DLSectionProperty::calcJ2_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_v){
    J2_v = DIM4*((((1.0/(DIM4*DIM4*DIM4*DIM4)*pow(DIM2-DIM3/2.0,4.0))/1.2E+1-1.0)*(DIM2*(2.1E+1/1.0E+2)-DIM3*(2.1E+1/2.0E+2)))/DIM4+1.0/3.0)*pow(DIM2-DIM3/2.0,3.0);
}

void MAST::Solid1DLSectionProperty::calcdJ2_v(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_v){
    dJ2_v = (1.0/(DIM4*DIM4*DIM4*DIM4)*dDIM2*pow(DIM2*2.0-DIM3,2.0)*(DIM2*(DIM3*DIM3*DIM3*DIM3)*7.0E+1-(DIM2*DIM2*DIM2*DIM2)*DIM3*5.6E+2-DIM2*(DIM4*DIM4*DIM4*DIM4)*1.344E+3+DIM3*(DIM4*DIM4*DIM4*DIM4)*6.72E+2+(DIM2*DIM2*DIM2*DIM2*DIM2)*2.24E+2-(DIM3*DIM3*DIM3*DIM3*DIM3)*7.0+(DIM4*DIM4*DIM4*DIM4*DIM4)*1.6E+3-(DIM2*DIM2)*(DIM3*DIM3*DIM3)*2.8E+2+(DIM2*DIM2*DIM2)*(DIM3*DIM3)*5.6E+2))/6.4E+3-(1.0/(DIM4*DIM4*DIM4*DIM4)*dDIM3*pow(DIM2*2.0-DIM3,2.0)*(DIM2*(DIM3*DIM3*DIM3*DIM3)*7.0E+1-(DIM2*DIM2*DIM2*DIM2)*DIM3*5.6E+2-DIM2*(DIM4*DIM4*DIM4*DIM4)*1.344E+3+DIM3*(DIM4*DIM4*DIM4*DIM4)*6.72E+2+(DIM2*DIM2*DIM2*DIM2*DIM2)*2.24E+2-(DIM3*DIM3*DIM3*DIM3*DIM3)*7.0+(DIM4*DIM4*DIM4*DIM4*DIM4)*1.6E+3-(DIM2*DIM2)*(DIM3*DIM3*DIM3)*2.8E+2+(DIM2*DIM2*DIM2)*(DIM3*DIM3)*5.6E+2))/1.28E+4+(1.0/(DIM4*DIM4*DIM4*DIM4*DIM4)*dDIM4*pow(DIM2*2.0-DIM3,3.0)*(DIM2*(DIM3*DIM3*DIM3*DIM3)*-2.1E+2+(DIM2*DIM2*DIM2*DIM2)*DIM3*1.68E+3-(DIM2*DIM2*DIM2*DIM2*DIM2)*6.72E+2+(DIM3*DIM3*DIM3*DIM3*DIM3)*2.1E+1+(DIM4*DIM4*DIM4*DIM4*DIM4)*3.2E+3+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*8.4E+2-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*1.68E+3))/7.68E+4;
}


void MAST::Solid1DLSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J)
{
    Real J_h, J_v;
    if((DIM1-0.5*DIM4)>DIM3)
    {
        calcJ1_h(DIM1, DIM2, DIM3, DIM4, J_h);
    }
    else
    {
        calcJ2_h(DIM1, DIM2, DIM3, DIM4, J_h);
    }
    
    if((DIM2-0.5*DIM3)>DIM4)
    {
        calcJ1_v(DIM1, DIM2, DIM3, DIM4, J_v);
    }
    else
    {
        calcJ2_v(DIM1, DIM2, DIM3, DIM4, J_v);
    }
    J = J_h+J_v;
}

void MAST::Solid1DLSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ)
{
    Real dJ_h, dJ_v;
    if((DIM1-0.5*DIM4)>DIM3)
    {
        calcdJ1_h(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_h);
    }
    else
    {
        calcdJ2_h(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_h);
    }
    
    if((DIM2-0.5*DIM3)>DIM4)
    {
        calcdJ1_v(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_v);
    }
    else
    {
        calcdJ2_v(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_v);
    }
    dJ = dJ_h+dJ_v;
}


void MAST::Solid1DLSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DLSectionProperty::calcA,
                                                              MAST::Solid1DLSectionProperty::calcdA,
                                                              DIM1, DIM2, DIM3, 
                                                              DIM4));
    
    _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DLSectionProperty::calcA,
                                                                MAST::Solid1DLSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DLSectionProperty::calcA,
                                                                MAST::Solid1DLSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DLSectionProperty::calcJ,
                                                                MAST::Solid1DLSectionProperty::calcdJ,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4));
    
    _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DLSectionProperty::calcIp,
                                                                MAST::Solid1DLSectionProperty::calcdIp,
                                                                MAST::Solid1DLSectionProperty::calcA,
                                                                MAST::Solid1DLSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DLSectionProperty::calcIz,
                                                                MAST::Solid1DLSectionProperty::calcdIz,
                                                                MAST::Solid1DLSectionProperty::calcIy,
                                                                MAST::Solid1DLSectionProperty::calcdIy,
                                                                MAST::Solid1DLSectionProperty::calcA,
                                                                MAST::Solid1DLSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
