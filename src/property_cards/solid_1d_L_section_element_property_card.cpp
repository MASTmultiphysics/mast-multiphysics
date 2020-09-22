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


const std::vector<libMesh::Point> 
MAST::Solid1DLSectionElementPropertyCard::get_geom_points(const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v, DIM3v, DIM4v;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v); DIM3(p, t, DIM3v); DIM4(p, t, DIM4v);
        
    Real offset_y, offset_z;
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    libMesh::Point shift(0.5*DIM4v, 0.5*DIM3v);
    
    std::vector<libMesh::Point> points = {
        libMesh::Point(0., 0.) - shift + offset,
        libMesh::Point(DIM1v, 0.) - shift + offset,
        libMesh::Point(DIM1v, DIM3v) - shift + offset,
        libMesh::Point(DIM4v, DIM3v) - shift + offset,
        libMesh::Point(DIM4v, DIM2v) - shift + offset,
        libMesh::Point(0., DIM2v) - shift + offset
    };
    
    return points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DLSectionElementPropertyCard::get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v, DIM3v, DIM4v;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v); DIM3(p, t, DIM3v); DIM4(p, t, DIM4v);
    
    Real dDIM1v, dDIM2v, dDIM3v, dDIM4v;
    DIM1.derivative(f, p, t, dDIM1v);   DIM2.derivative(f, p, t, dDIM2v);
    DIM3.derivative(f, p, t, dDIM3v);   DIM4.derivative(f, p, t, dDIM4v);
    
    Real offset_y, offset_z, doffset_y, doffset_z;
    hy_off(p, t, offset_y);     hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z);     hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    libMesh::Point dshift(0.5*dDIM4v, 0.5*dDIM3v);
    
    std::vector<libMesh::Point> points = {
        libMesh::Point(0., 0.) - dshift + doffset,
        libMesh::Point(dDIM1v, 0.) - dshift + doffset,
        libMesh::Point(dDIM1v, dDIM3v) - dshift + doffset,
        libMesh::Point(dDIM4v, dDIM3v) - dshift + doffset,
        libMesh::Point(dDIM4v, dDIM2v) - dshift + doffset,
        libMesh::Point(0., dDIM2v) - dshift + doffset
    };
    
    return points;
}


const libMesh::Point
MAST::Solid1DLSectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_centroid(p, t);
}


const libMesh::Point
MAST::Solid1DLSectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_shear_center(p, t);
}


const libMesh::Point
MAST::Solid1DLSectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_centroid_derivative(f, p, t);
}


const libMesh::Point 
MAST::Solid1DLSectionElementPropertyCard::get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t) 
{
    return cross_section->get_shear_center_derivative(f, p, t);
}


const std::vector<libMesh::Point> 
MAST::Solid1DLSectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    // Ordered C, D, E, F as defined in MSC Nastran manual (See PBARL or PBEAML)
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, DIM3v, DIM4v;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v); DIM3(p, t, DIM3v); DIM4(p, t, DIM4v);
    
    Real offset_y, offset_z;
    hy_off(p, t, offset_y);     hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    // Stress Points
    std::vector<libMesh::Point> points = {
        libMesh::Point(0.5*DIM4v, DIM2v-0.5*DIM3v) - ps,
        libMesh::Point(DIM1v-0.5*DIM4v, -0.5*DIM3v)  - ps,
        libMesh::Point(-0.5*DIM4v, -0.5*DIM3v) - ps,
        libMesh::Point(-0.5*DIM4v, DIM2v-0.5*DIM3v) - ps
    };
    
    return points;
};


const std::vector<libMesh::Point> 
MAST::Solid1DLSectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v, DIM3v, DIM4v;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v); DIM3(p, t, DIM3v); DIM4(p, t, DIM4v);
    
    Real dDIM1v, dDIM2v, dDIM3v, dDIM4v;
    DIM1.derivative(f, p, t, dDIM1v);   DIM2.derivative(f, p, t, dDIM2v);
    DIM3.derivative(f, p, t, dDIM3v);   DIM4.derivative(f, p, t, dDIM4v);
    
    Real offset_y, offset_z,    doffset_y, doffset_z;
    hy_off(p, t, offset_y);     hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z);     hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y); 
    
    // Stress Points
    std::vector<libMesh::Point> points = {
        libMesh::Point(0.5*dDIM4v, dDIM2v-0.5*dDIM3v) - dps,
        libMesh::Point(dDIM1v-0.5*dDIM4v, -0.5*dDIM3v) - dps,
        libMesh::Point(-0.5*dDIM4v, -0.5*dDIM3v) - dps,
        libMesh::Point(-0.5*dDIM4v, dDIM2v-0.5*dDIM3v) - dps
    };
}


void MAST::Solid1DLSectionElementPropertyCard::init(const libMesh::LibMeshInit& init,
                                                    const uint n_target_elems,
                                                    const libMesh::ElemType element_type) 
{
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Create a cross section model of this section
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, *this, element_type));
    
//     _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DLSectionProperty::calcA,
//                                                               MAST::Solid1DLSectionProperty::calcdA,
//                                                               DIM1, DIM2, DIM3, 
//                                                               DIM4));
//     
//     _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
//                                                                 MAST::Solid1DLSectionProperty::calcA,
//                                                                 MAST::Solid1DLSectionProperty::calcdA,
//                                                                 DIM1, DIM2, 
//                                                                 DIM3, DIM4,
//                                                                 hz_off));
//     
//     _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
//                                                                 MAST::Solid1DLSectionProperty::calcA,
//                                                                 MAST::Solid1DLSectionProperty::calcdA,
//                                                                 DIM1, DIM2, 
//                                                                 DIM3, DIM4,
//                                                                 hy_off));
//     
//     _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
//                                                                 MAST::Solid1DLSectionProperty::calcJ,
//                                                                 MAST::Solid1DLSectionProperty::calcdJ,
//                                                                 DIM1, DIM2,
//                                                                 DIM3, DIM4));
//     
//     _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
//                                                                 MAST::Solid1DLSectionProperty::calcIp,
//                                                                 MAST::Solid1DLSectionProperty::calcdIp,
//                                                                 MAST::Solid1DLSectionProperty::calcA,
//                                                                 MAST::Solid1DLSectionProperty::calcdA,
//                                                                 DIM1, DIM2,
//                                                                 DIM3, DIM4,
//                                                                 hy_off,
//                                                                 hz_off));
//     
//     _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
//                                                                 MAST::Solid1DLSectionProperty::calcIz,
//                                                                 MAST::Solid1DLSectionProperty::calcdIz,
//                                                                 MAST::Solid1DLSectionProperty::calcIy,
//                                                                 MAST::Solid1DLSectionProperty::calcdIy,
//                                                                 MAST::Solid1DLSectionProperty::calcA,
//                                                                 MAST::Solid1DLSectionProperty::calcdA,
//                                                                 DIM1, DIM2,
//                                                                 DIM3, DIM4,
//                                                                 hy_off,
//                                                                 hz_off));
    
    _A.reset(new MAST::Solid1DnParameterSectionProperty::Area(*cross_section));
    
    _Ay.reset(new MAST::Solid1DnParameterSectionProperty::AreaYMoment(*cross_section));
    
    _Az.reset(new MAST::Solid1DnParameterSectionProperty::AreaZMoment(*cross_section));
    
    _AI.reset(new MAST::Solid1DnParameterSectionProperty::AreaInertiaMatrix(*cross_section));
    
    _Ip.reset(new MAST::Solid1DnParameterSectionProperty::PolarInertia(*cross_section));
    
    _J.reset(new MAST::Solid1DnParameterSectionProperty::TorsionalConstant(*cross_section));
    
    _Kappa.reset(new MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix(*cross_section));
    
    _Gamma.reset(new MAST::Solid1DnParameterSectionProperty::WarpingConstant(*cross_section));
    
    _initialized = true;
}


void MAST::Solid1DLSectionElementPropertyCard::create_cross_section(
    const libMesh::LibMeshInit& init, const uint n_target_elems, 
    const libMesh::ElemType element_type)
{
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, 
                                               *this, element_type));
}


void MAST::Solid1DLSectionElementPropertyCard::calculate_properties_pilkey()
{
    cross_section->calculate_geometric_properties();
    cross_section->calculate_warping_properties();
}
