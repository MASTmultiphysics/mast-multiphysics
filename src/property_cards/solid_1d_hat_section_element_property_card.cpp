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
#include "property_cards/solid_1d_hat_section_element_property_card.h"


// Open HAT (HAT in Nastran)
void MAST::Solid1DHatSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A){
    A = DIM2*(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0);
}

void MAST::Solid1DHatSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA){
    dA = DIM2*dDIM1*2.0+DIM2*dDIM3+DIM2*dDIM4*2.0+dDIM2*(DIM1*2.0-DIM2*4.0+DIM3+DIM4*2.0);
}

void MAST::Solid1DHatSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz){
    Iz = (DIM2*(DIM1*(DIM2*DIM2*DIM2)*-1.6E+1-(DIM1*DIM1*DIM1)*DIM2*1.6E+1+(DIM1*DIM1*DIM1)*DIM3*8.0+(DIM1*DIM1*DIM1)*DIM4*1.6E+1-(DIM2*DIM2*DIM2)*DIM3*4.0-(DIM2*DIM2*DIM2)*DIM4*5.6E+1+(DIM1*DIM1*DIM1*DIM1)*4.0+(DIM2*DIM2*DIM2*DIM2)*4.0+(DIM1*DIM1)*(DIM2*DIM2)*2.4E+1+(DIM2*DIM2)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4)*4.0+DIM1*(DIM2*DIM2)*DIM3*8.0-(DIM1*DIM1)*DIM2*DIM3*1.2E+1+DIM1*(DIM2*DIM2)*DIM4*1.12E+2-(DIM1*DIM1)*DIM2*DIM4*7.2E+1+(DIM1*DIM1)*DIM3*DIM4*2.4E+1+(DIM2*DIM2)*DIM3*DIM4*2.8E+1-DIM1*DIM2*DIM3*DIM4*4.8E+1))/(DIM1*2.4E+1-DIM2*2.4E+1+DIM3*1.2E+1+DIM4*2.4E+1);
}

void MAST::Solid1DHatSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz){
    dIz = (dDIM2*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2*DIM2)*1.36E+2-(DIM1*DIM1*DIM1*DIM1)*DIM2*6.4E+1+(DIM1*DIM1*DIM1*DIM1)*DIM3*2.0E+1+(DIM1*DIM1*DIM1*DIM1)*DIM4*4.0E+1+(DIM2*DIM2*DIM2*DIM2)*DIM3*4.4E+1+(DIM2*DIM2*DIM2*DIM2)*DIM4*3.76E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*8.0-(DIM2*DIM2*DIM2*DIM2*DIM2)*3.2E+1-(DIM1*DIM1)*(DIM2*DIM2*DIM2)*2.24E+2+(DIM1*DIM1*DIM1)*(DIM2*DIM2)*1.76E+2+(DIM1*DIM1*DIM1)*(DIM3*DIM3)*8.0+(DIM1*DIM1*DIM1)*(DIM4*DIM4)*3.2E+1+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*3.0-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*2.0E+1+(DIM2*DIM2)*(DIM4*DIM4*DIM4)*2.4E+1-(DIM2*DIM2*DIM2)*(DIM4*DIM4)*4.64E+2-DIM1*(DIM2*DIM2*DIM2)*DIM3*1.28E+2-(DIM1*DIM1*DIM1)*DIM2*DIM3*8.0E+1-DIM1*(DIM2*DIM2*DIM2)*DIM4*1.024E+3-(DIM1*DIM1*DIM1)*DIM2*DIM4*3.52E+2+(DIM1*DIM1*DIM1)*DIM3*DIM4*8.0E+1-(DIM2*DIM2*DIM2)*DIM3*DIM4*3.68E+2+DIM1*(DIM2*DIM2)*(DIM3*DIM3)*3.0E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*2.4E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*1.44E+2+DIM1*(DIM2*DIM2)*(DIM4*DIM4)*6.96E+2-(DIM1*DIM1)*DIM2*(DIM4*DIM4)*2.88E+2+(DIM1*DIM1)*(DIM2*DIM2)*DIM4*9.6E+2+(DIM1*DIM1)*DIM3*(DIM4*DIM4)*4.8E+1+(DIM1*DIM1)*(DIM3*DIM3)*DIM4*2.4E+1+(DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.8E+2+(DIM2*DIM2)*(DIM3*DIM3)*DIM4*9.0E+1-DIM1*DIM2*DIM3*(DIM4*DIM4)*1.92E+2-DIM1*DIM2*(DIM3*DIM3)*DIM4*9.6E+1+DIM1*(DIM2*DIM2)*DIM3*DIM4*6.48E+2-(DIM1*DIM1)*DIM2*DIM3*DIM4*3.84E+2))/1.2E+1+(DIM2*dDIM3*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-8.0-(DIM1*DIM1*DIM1)*DIM2*2.4E+1+(DIM1*DIM1*DIM1)*DIM4*4.8E+1-(DIM2*DIM2*DIM2)*DIM3*4.0-(DIM2*DIM2*DIM2)*DIM4*8.0+(DIM1*DIM1*DIM1*DIM1)*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*4.0+(DIM1*DIM1)*(DIM2*DIM2)*1.6E+1+(DIM1*DIM1)*(DIM4*DIM4)*4.8E+1+(DIM2*DIM2)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4)*5.2E+1+DIM1*(DIM2*DIM2)*DIM3*4.0-DIM1*DIM2*(DIM4*DIM4)*9.6E+1+DIM1*(DIM2*DIM2)*DIM4*5.6E+1-(DIM1*DIM1)*DIM2*DIM4*9.6E+1+(DIM2*DIM2)*DIM3*DIM4*4.0))/1.2E+1+(DIM2*dDIM4*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-1.52E+2-(DIM1*DIM1*DIM1)*DIM2*7.2E+1+(DIM1*DIM1*DIM1)*DIM3*2.4E+1-(DIM2*DIM2*DIM2)*DIM3*5.2E+1-(DIM2*DIM2*DIM2)*DIM4*8.0+(DIM1*DIM1*DIM1*DIM1)*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*5.2E+1+(DIM1*DIM1)*(DIM2*DIM2)*1.6E+2+(DIM1*DIM1)*(DIM3*DIM3)*1.2E+1+(DIM2*DIM2)*(DIM3*DIM3)*1.3E+1+(DIM2*DIM2)*(DIM4*DIM4)*4.0-DIM1*DIM2*(DIM3*DIM3)*2.4E+1+DIM1*(DIM2*DIM2)*DIM3*1.24E+2-(DIM1*DIM1)*DIM2*DIM3*9.6E+1+DIM1*(DIM2*DIM2)*DIM4*8.0+(DIM2*DIM2)*DIM3*DIM4*4.0))/6.0+(DIM2*dDIM1*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-1.6E+1-(DIM1*DIM1*DIM1)*DIM2*1.6E+1+(DIM1*DIM1*DIM1)*DIM3*8.0+(DIM1*DIM1*DIM1)*DIM4*1.6E+1-(DIM2*DIM2*DIM2)*DIM3*4.0-(DIM2*DIM2*DIM2)*DIM4*2.4E+1+(DIM1*DIM1*DIM1*DIM1)*4.0+(DIM2*DIM2*DIM2*DIM2)*4.0+(DIM1*DIM1)*(DIM2*DIM2)*2.4E+1+(DIM1*DIM1)*(DIM3*DIM3)*4.0+(DIM1*DIM1)*(DIM4*DIM4)*1.6E+1+(DIM2*DIM2)*(DIM3*DIM3)+(DIM2*DIM2)*(DIM4*DIM4)*3.6E+1-DIM1*DIM2*(DIM3*DIM3)*4.0+DIM1*(DIM2*DIM2)*DIM3*1.6E+1-(DIM1*DIM1)*DIM2*DIM3*2.0E+1-DIM1*DIM2*(DIM4*DIM4)*4.8E+1+DIM1*(DIM2*DIM2)*DIM4*6.4E+1-(DIM1*DIM1)*DIM2*DIM4*5.6E+1+DIM1*DIM3*(DIM4*DIM4)*1.6E+1+DIM1*(DIM3*DIM3)*DIM4*8.0+(DIM1*DIM1)*DIM3*DIM4*2.4E+1-DIM2*DIM3*(DIM4*DIM4)*1.6E+1-DIM2*(DIM3*DIM3)*DIM4*8.0+(DIM2*DIM2)*DIM3*DIM4*2.8E+1-DIM1*DIM2*DIM3*DIM4*4.8E+1))/2.0;
}

void MAST::Solid1DHatSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy){
    Iy = (DIM2*(DIM3*DIM3*DIM3))/1.2E+1+(DIM2*pow(DIM2+DIM4,3.0))/6.0+DIM2*(DIM2+DIM4)*pow(DIM2*(-1.0/2.0)+DIM3/2.0+DIM4/2.0,2.0)*2.0+(DIM2*(DIM1-DIM2*2.0)*((DIM2*DIM2)*4.0+(DIM3*DIM3)*3.0-DIM2*DIM3*6.0))/6.0;
}

void MAST::Solid1DHatSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy){
    dIy = (dDIM2*(DIM1*(DIM2*DIM2)*2.4E+1+DIM1*(DIM3*DIM3)*6.0-DIM2*(DIM3*DIM3)*1.2E+1+(DIM2*DIM2)*DIM3*3.6E+1+DIM3*(DIM4*DIM4)*1.2E+1+(DIM3*DIM3)*DIM4*6.0-(DIM2*DIM2*DIM2)*3.2E+1+DIM3*DIM3*DIM3+(DIM4*DIM4*DIM4)*8.0-DIM1*DIM2*DIM3*2.4E+1))/1.2E+1+dDIM3*((DIM2*(DIM3*DIM3))/4.0-(DIM2*(DIM1-DIM2*2.0)*(DIM2*6.0-DIM3*6.0))/6.0+DIM2*(DIM2+DIM4)*(DIM2*(-1.0/2.0)+DIM3/2.0+DIM4/2.0)*2.0)+(DIM2*dDIM1*((DIM2*DIM2)*4.0+(DIM3*DIM3)*3.0-DIM2*DIM3*6.0))/6.0+(DIM2*dDIM4*pow(DIM3+DIM4*2.0,2.0))/2.0;
}

void MAST::Solid1DHatSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip){
    Ip = (DIM2*(DIM1*(DIM2*DIM2*DIM2)*-4.8E+1-(DIM1*DIM1*DIM1)*DIM2*1.6E+1+DIM1*(DIM3*DIM3*DIM3)*8.0+(DIM1*DIM1*DIM1)*DIM3*8.0+DIM1*(DIM4*DIM4*DIM4)*1.6E+1-DIM2*(DIM3*DIM3*DIM3)*8.0+(DIM1*DIM1*DIM1)*DIM4*1.6E+1-(DIM2*DIM2*DIM2)*DIM3*3.6E+1-DIM2*(DIM4*DIM4*DIM4)*1.6E+1-(DIM2*DIM2*DIM2)*DIM4*7.2E+1+DIM3*(DIM4*DIM4*DIM4)*3.2E+1+(DIM3*DIM3*DIM3)*DIM4*8.0+(DIM1*DIM1*DIM1*DIM1)*4.0+(DIM2*DIM2*DIM2*DIM2)*2.0E+1+DIM3*DIM3*DIM3*DIM3+(DIM4*DIM4*DIM4*DIM4)*1.6E+1+(DIM1*DIM1)*(DIM2*DIM2)*4.0E+1+(DIM1*DIM1)*(DIM3*DIM3)*1.2E+1+(DIM2*DIM2)*(DIM3*DIM3)*2.5E+1+(DIM2*DIM2)*(DIM4*DIM4)*4.0+(DIM3*DIM3)*(DIM4*DIM4)*2.4E+1-DIM1*DIM2*(DIM3*DIM3)*3.6E+1+DIM1*(DIM2*DIM2)*DIM3*6.4E+1-(DIM1*DIM1)*DIM2*DIM3*3.6E+1+DIM1*(DIM2*DIM2)*DIM4*1.28E+2-(DIM1*DIM1)*DIM2*DIM4*7.2E+1+DIM1*DIM3*(DIM4*DIM4)*2.4E+1+DIM1*(DIM3*DIM3)*DIM4*2.4E+1+(DIM1*DIM1)*DIM3*DIM4*2.4E+1-DIM2*DIM3*(DIM4*DIM4)*2.4E+1-DIM2*(DIM3*DIM3)*DIM4*2.4E+1+(DIM2*DIM2)*DIM3*DIM4*5.2E+1-DIM1*DIM2*DIM3*DIM4*7.2E+1))/(DIM1*2.4E+1-DIM2*2.4E+1+DIM3*1.2E+1+DIM4*2.4E+1);
}

void MAST::Solid1DHatSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp){
    dIp = (dDIM2*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2*DIM2)*4.88E+2-(DIM1*DIM1*DIM1*DIM1)*DIM2*6.4E+1+DIM1*(DIM3*DIM3*DIM3*DIM3)*1.0E+1+(DIM1*DIM1*DIM1*DIM1)*DIM3*2.0E+1+DIM1*(DIM4*DIM4*DIM4*DIM4)*6.4E+1-DIM2*(DIM3*DIM3*DIM3*DIM3)*1.6E+1+(DIM1*DIM1*DIM1*DIM1)*DIM4*4.0E+1+(DIM2*DIM2*DIM2*DIM2)*DIM3*3.16E+2-DIM2*(DIM4*DIM4*DIM4*DIM4)*6.4E+1+(DIM2*DIM2*DIM2*DIM2)*DIM4*6.32E+2+DIM3*(DIM4*DIM4*DIM4*DIM4)*8.0E+1+(DIM3*DIM3*DIM3*DIM3)*DIM4*1.0E+1+(DIM1*DIM1*DIM1*DIM1*DIM1)*8.0-(DIM2*DIM2*DIM2*DIM2*DIM2)*1.6E+2+DIM3*DIM3*DIM3*DIM3*DIM3+(DIM4*DIM4*DIM4*DIM4*DIM4)*3.2E+1-(DIM1*DIM1)*(DIM2*DIM2*DIM2)*5.44E+2+(DIM1*DIM1*DIM1)*(DIM2*DIM2)*2.72E+2+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*2.8E+1+(DIM1*DIM1*DIM1)*(DIM3*DIM3)*3.2E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4)*3.2E+1+(DIM1*DIM1*DIM1)*(DIM4*DIM4)*3.2E+1+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*9.1E+1-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*2.44E+2+(DIM2*DIM2)*(DIM4*DIM4*DIM4)*5.6E+1-(DIM2*DIM2*DIM2)*(DIM4*DIM4)*5.92E+2+(DIM3*DIM3)*(DIM4*DIM4*DIM4)*8.0E+1+(DIM3*DIM3*DIM3)*(DIM4*DIM4)*4.0E+1-DIM1*DIM2*(DIM3*DIM3*DIM3)*1.04E+2-DIM1*(DIM2*DIM2*DIM2)*DIM3*7.36E+2-(DIM1*DIM1*DIM1)*DIM2*DIM3*1.76E+2-DIM1*DIM2*(DIM4*DIM4*DIM4)*6.4E+1-DIM1*(DIM2*DIM2*DIM2)*DIM4*1.472E+3-(DIM1*DIM1*DIM1)*DIM2*DIM4*3.52E+2+DIM1*DIM3*(DIM4*DIM4*DIM4)*1.28E+2+DIM1*(DIM3*DIM3*DIM3)*DIM4*5.6E+1+(DIM1*DIM1*DIM1)*DIM3*DIM4*8.0E+1-DIM2*DIM3*(DIM4*DIM4*DIM4)*1.28E+2-DIM2*(DIM3*DIM3*DIM3)*DIM4*8.0E+1-(DIM2*DIM2*DIM2)*DIM3*DIM4*7.84E+2+DIM1*(DIM2*DIM2)*(DIM3*DIM3)*4.14E+2-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*2.16E+2+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*5.76E+2+DIM1*(DIM2*DIM2)*(DIM4*DIM4)*7.92E+2-(DIM1*DIM1)*DIM2*(DIM4*DIM4)*2.88E+2+(DIM1*DIM1)*(DIM2*DIM2)*DIM4*1.152E+3+DIM1*(DIM3*DIM3)*(DIM4*DIM4)*1.2E+2+(DIM1*DIM1)*DIM3*(DIM4*DIM4)*9.6E+1+(DIM1*DIM1)*(DIM3*DIM3)*DIM4*9.6E+1-DIM2*(DIM3*DIM3)*(DIM4*DIM4)*1.44E+2+(DIM2*DIM2)*DIM3*(DIM4*DIM4)*3.72E+2+(DIM2*DIM2)*(DIM3*DIM3)*DIM4*3.54E+2-DIM1*DIM2*DIM3*(DIM4*DIM4)*3.84E+2-DIM1*DIM2*(DIM3*DIM3)*DIM4*3.84E+2+DIM1*(DIM2*DIM2)*DIM3*DIM4*1.224E+3-(DIM1*DIM1)*DIM2*DIM3*DIM4*5.76E+2))/1.2E+1+(DIM2*dDIM3*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-1.52E+2-(DIM1*DIM1*DIM1)*DIM2*7.2E+1+DIM1*(DIM3*DIM3*DIM3)*2.4E+1+(DIM1*DIM1*DIM1)*DIM3*4.8E+1+DIM1*(DIM4*DIM4*DIM4)*9.6E+1-DIM2*(DIM3*DIM3*DIM3)*2.4E+1+(DIM1*DIM1*DIM1)*DIM4*4.8E+1-(DIM2*DIM2*DIM2)*DIM3*1.0E+2-DIM2*(DIM4*DIM4*DIM4)*9.6E+1-(DIM2*DIM2*DIM2)*DIM4*1.04E+2+DIM3*(DIM4*DIM4*DIM4)*9.6E+1+(DIM3*DIM3*DIM3)*DIM4*2.4E+1+(DIM1*DIM1*DIM1*DIM1)*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*5.2E+1+(DIM3*DIM3*DIM3*DIM3)*3.0+(DIM4*DIM4*DIM4*DIM4)*4.8E+1+(DIM1*DIM1)*(DIM2*DIM2)*1.6E+2+(DIM1*DIM1)*(DIM3*DIM3)*6.0E+1+(DIM1*DIM1)*(DIM4*DIM4)*9.6E+1+(DIM2*DIM2)*(DIM3*DIM3)*7.3E+1+(DIM2*DIM2)*(DIM4*DIM4)*1.48E+2+(DIM3*DIM3)*(DIM4*DIM4)*7.2E+1-DIM1*DIM2*(DIM3*DIM3)*1.32E+2+DIM1*(DIM2*DIM2)*DIM3*2.44E+2-(DIM1*DIM1)*DIM2*DIM3*1.92E+2-DIM1*DIM2*(DIM4*DIM4)*2.4E+2+DIM1*(DIM2*DIM2)*DIM4*2.48E+2-(DIM1*DIM1)*DIM2*DIM4*1.92E+2+DIM1*DIM3*(DIM4*DIM4)*1.92E+2+DIM1*(DIM3*DIM3)*DIM4*1.2E+2+(DIM1*DIM1)*DIM3*DIM4*1.44E+2-DIM2*DIM3*(DIM4*DIM4)*1.92E+2-DIM2*(DIM3*DIM3)*DIM4*1.2E+2+(DIM2*DIM2)*DIM3*DIM4*1.96E+2-DIM1*DIM2*DIM3*DIM4*3.36E+2))/1.2E+1+(DIM2*dDIM1*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-8.0E+1-(DIM1*DIM1*DIM1)*DIM2*4.8E+1+DIM1*(DIM3*DIM3*DIM3)*1.2E+1+(DIM1*DIM1*DIM1)*DIM3*2.4E+1-DIM2*(DIM3*DIM3*DIM3)*1.8E+1+(DIM1*DIM1*DIM1)*DIM4*4.8E+1-(DIM2*DIM2*DIM2)*DIM3*5.2E+1-(DIM2*DIM2*DIM2)*DIM4*1.04E+2+(DIM3*DIM3*DIM3)*DIM4*1.2E+1+(DIM1*DIM1*DIM1*DIM1)*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*2.8E+1+(DIM3*DIM3*DIM3*DIM3)*3.0+(DIM1*DIM1)*(DIM2*DIM2)*8.8E+1+(DIM1*DIM1)*(DIM3*DIM3)*2.4E+1+(DIM1*DIM1)*(DIM4*DIM4)*4.8E+1+(DIM2*DIM2)*(DIM3*DIM3)*4.3E+1+(DIM2*DIM2)*(DIM4*DIM4)*1.24E+2+(DIM3*DIM3)*(DIM4*DIM4)*1.2E+1-DIM1*DIM2*(DIM3*DIM3)*6.0E+1+DIM1*(DIM2*DIM2)*DIM3*1.12E+2-(DIM1*DIM1)*DIM2*DIM3*8.4E+1-DIM1*DIM2*(DIM4*DIM4)*1.44E+2+DIM1*(DIM2*DIM2)*DIM4*2.24E+2-(DIM1*DIM1)*DIM2*DIM4*1.68E+2+DIM1*DIM3*(DIM4*DIM4)*4.8E+1+DIM1*(DIM3*DIM3)*DIM4*4.8E+1+(DIM1*DIM1)*DIM3*DIM4*7.2E+1-DIM2*DIM3*(DIM4*DIM4)*7.2E+1-DIM2*(DIM3*DIM3)*DIM4*7.2E+1+(DIM2*DIM2)*DIM3*DIM4*1.48E+2-DIM1*DIM2*DIM3*DIM4*1.92E+2))/6.0+(DIM2*dDIM4*1.0/pow(DIM1*2.0-DIM2*2.0+DIM3+DIM4*2.0,2.0)*(DIM1*(DIM2*DIM2*DIM2)*-1.52E+2-(DIM1*DIM1*DIM1)*DIM2*7.2E+1+DIM1*(DIM3*DIM3*DIM3)*1.2E+1+(DIM1*DIM1*DIM1)*DIM3*2.4E+1+DIM1*(DIM4*DIM4*DIM4)*9.6E+1-DIM2*(DIM3*DIM3*DIM3)*1.2E+1-(DIM2*DIM2*DIM2)*DIM3*5.2E+1-DIM2*(DIM4*DIM4*DIM4)*9.6E+1-(DIM2*DIM2*DIM2)*DIM4*8.0+DIM3*(DIM4*DIM4*DIM4)*9.6E+1+(DIM3*DIM3*DIM3)*DIM4*2.4E+1+(DIM1*DIM1*DIM1*DIM1)*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*5.2E+1+(DIM3*DIM3*DIM3*DIM3)*3.0+(DIM4*DIM4*DIM4*DIM4)*4.8E+1+(DIM1*DIM1)*(DIM2*DIM2)*1.6E+2+(DIM1*DIM1)*(DIM3*DIM3)*2.4E+1+(DIM1*DIM1)*(DIM4*DIM4)*4.8E+1+(DIM2*DIM2)*(DIM3*DIM3)*2.5E+1+(DIM2*DIM2)*(DIM4*DIM4)*5.2E+1+(DIM3*DIM3)*(DIM4*DIM4)*7.2E+1-DIM1*DIM2*(DIM3*DIM3)*4.8E+1+DIM1*(DIM2*DIM2)*DIM3*1.24E+2-(DIM1*DIM1)*DIM2*DIM3*9.6E+1-DIM1*DIM2*(DIM4*DIM4)*9.6E+1+DIM1*(DIM2*DIM2)*DIM4*8.0+DIM1*DIM3*(DIM4*DIM4)*1.44E+2+DIM1*(DIM3*DIM3)*DIM4*7.2E+1+(DIM1*DIM1)*DIM3*DIM4*4.8E+1-DIM2*DIM3*(DIM4*DIM4)*1.44E+2-DIM2*(DIM3*DIM3)*DIM4*7.2E+1+(DIM2*DIM2)*DIM3*DIM4*5.2E+1-DIM1*DIM2*DIM3*DIM4*9.6E+1))/6.0;
}

void MAST::Solid1DHatSectionProperty::calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_w){
    J1_w = (DIM2*DIM2*DIM2)*(DIM1-DIM2*2.0);
}

void MAST::Solid1DHatSectionProperty::calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_w){
    dJ1_w = (DIM2*DIM2*DIM2)*dDIM1+(DIM2*DIM2)*dDIM2*(DIM1*3.0-DIM2*8.0);
}

void MAST::Solid1DHatSectionProperty::calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_w){
    J2_w = DIM2*pow(DIM1-DIM2*2.0,3.0);
}

void MAST::Solid1DHatSectionProperty::calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_w){
    dJ2_w = DIM2*dDIM1*pow(DIM1-DIM2*2.0,2.0)*3.0+dDIM2*pow(DIM1-DIM2*2.0,2.0)*(DIM1-DIM2*8.0);
}

void MAST::Solid1DHatSectionProperty::calcJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_ft){
    J1_ft = (DIM2*DIM2*DIM2)*DIM3;
}

void MAST::Solid1DHatSectionProperty::calcdJ1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_ft){
    dJ1_ft = (DIM2*DIM2*DIM2)*dDIM3+(DIM2*DIM2)*DIM3*dDIM2*3.0;
}

void MAST::Solid1DHatSectionProperty::calcJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_ft){
    J2_ft = DIM2*(DIM3*DIM3*DIM3);
}

void MAST::Solid1DHatSectionProperty::calcdJ2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_ft){
    dJ2_ft = (DIM3*DIM3*DIM3)*dDIM2+DIM2*(DIM3*DIM3)*dDIM3*3.0;
}

void MAST::Solid1DHatSectionProperty::calcJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_fb){
    J1_fb = (DIM2*DIM2*DIM2)*(DIM2+DIM4);
}

void MAST::Solid1DHatSectionProperty::calcdJ1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_fb){
    dJ1_fb = (DIM2*DIM2*DIM2)*dDIM4+(DIM2*DIM2)*dDIM2*(DIM2*4.0+DIM4*3.0);
}

void MAST::Solid1DHatSectionProperty::calcJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_fb){
    J2_fb = DIM2*pow(DIM2+DIM4,3.0);
}

void MAST::Solid1DHatSectionProperty::calcdJ2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_fb){
    dJ2_fb = DIM2*dDIM4*pow(DIM2+DIM4,2.0)*3.0+dDIM2*(DIM2*4.0+DIM4)*pow(DIM2+DIM4,2.0);
}

void MAST::Solid1DHatSectionProperty::calck1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_fb){
    k1_fb = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk1_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_fb){
    dk1_fb = 0.0;
}

void MAST::Solid1DHatSectionProperty::calck2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_fb){
    k2_fb = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk2_fb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_fb){
    dk2_fb = 0.0;
}

void MAST::Solid1DHatSectionProperty::calck1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_ft){
    k1_ft = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk1_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_ft){
    dk1_ft = 0.0;
}

void MAST::Solid1DHatSectionProperty::calck2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_ft){
    k2_ft = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk2_ft(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_ft){
    dk2_ft = 0.0;
}

void MAST::Solid1DHatSectionProperty::calck1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_w){
    k1_w = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_w){
    dk1_w = 0.0;
}

void MAST::Solid1DHatSectionProperty::calck2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_w){
    k2_w = 1.0/3.0;
}

void MAST::Solid1DHatSectionProperty::calcdk2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_w){
    dk2_w = 0.0;
}

void MAST::Solid1DHatSectionProperty::calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Jc){
    Jc = DIM2*((DIM2*DIM2*DIM2)*3.73826991E+8-sqrt(2.0)*pow(DIM2*DIM2,3.0/2.0)*2.64043776E+8)*8.935219657483246E-7;
}

void MAST::Solid1DHatSectionProperty::calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJc){
    dJc = dDIM2*((DIM2*DIM2*DIM2)*1.24608997E+8-sqrt(2.0)*pow(DIM2*DIM2,3.0/2.0)*8.8014592E+7)*1.07222635889799E-5;
}



// TODO: Improve accuracy of calculation of torsion constant and corresponding derivative
void MAST::Solid1DHatSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J)
{
    Real t_fb = DIM2;
    Real w_fb = DIM4+DIM2;
    
    Real t_ft = DIM2;
    Real w_ft = DIM3;
    
    Real h_w = DIM1-2*DIM2;
    Real t_w = DIM2;
    
    Real k_ft, k_fb, k_w, c, J_w, J_ft, J_fb;
    Real wf = t_w/t_ft;
    
    // Sum of rectangles with correction of El Darwish and Johnson
    if ( (wf>0.5) and (wf<1.0) )
    {
        k_ft = 0.3333333333333333;
        k_fb = 0.3333333333333333; 
        k_w = 0.3333333333333333;
        calcJc(DIM1, DIM2, DIM3, DIM4, c);
    }
    
    // Sum of rectangles method, with rectangles corrected for aspect ratios
    else
    {
        if (w_ft>t_ft){
            calck1_ft(DIM1, DIM2, DIM3, DIM4, k_ft);
        }
        else{
            calck2_ft(DIM1, DIM2, DIM3, DIM4, k_ft);
        }
        
        if (w_fb>t_fb){
            calck1_fb(DIM1, DIM2, DIM3, DIM4, k_fb);
        }
        else{
            calck2_fb(DIM1, DIM2, DIM3, DIM4, k_ft);
        }
        
        if (h_w>t_w){
            calck1_w(DIM1, DIM2, DIM3, DIM4, k_w);
        }
        else{
            calck2_w(DIM1, DIM2, DIM3, DIM4, k_w);
        }
        c = 0.0;
    }
    
    
    if (w_ft>t_ft){
        calcJ1_ft(DIM1, DIM2, DIM3, DIM4, J_ft);
    }
    else{
        calcJ2_ft(DIM1, DIM2, DIM3, DIM4, J_ft);
    }
    
    if (w_fb>t_fb){
        calcJ1_fb(DIM1, DIM2, DIM3, DIM4, J_fb);
    }
    else{
        calcJ2_fb(DIM1, DIM2, DIM3, DIM4, J_fb);
    }
    
    if (h_w>t_w){
        calcJ1_w(DIM1, DIM2, DIM3, DIM4, J_w);
    }
    else{
        calcJ2_w(DIM1, DIM2, DIM3, DIM4, J_w);
    }
    
    J = 2.0*k_fb*J_fb + 2.0*k_w*J_w + k_ft*J_ft + c;
}



void MAST::Solid1DHatSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ)
{
    Real t_fb = DIM2;
    Real w_fb = DIM4+DIM2;
    
    Real t_ft = DIM2;
    Real w_ft = DIM3;
    
    Real h_w = DIM1-2*DIM2;
    Real t_w = DIM2;
    
    Real k_ft, k_fb, k_w, c, J_w, J_ft, J_fb;
    Real wf = t_w/t_ft;
    
    Real dk_ft, dk_fb, dk_w, dc, dJ_w, dJ_ft, dJ_fb;
    
    
    // Sum of rectangles with correction of El Darwish and Johnson
    if ( (wf>0.5) and (wf<1.0) )
    {
        k_ft = 0.3333333333333333;
        k_fb = 0.3333333333333333;
        k_w = 0.3333333333333333;
        dk_ft = 0.0;
        dk_fb = 0.0;
        dk_w = 0.0;
        calcdJc(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dc);
    }
    
    // Sum of rectangles method, with rectangles corrected for aspect ratios
    else
    {
        if (w_ft>t_ft){
            calck1_ft(DIM1, DIM2, DIM3, DIM4, k_ft);
            calcdk1_ft(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_ft);
        }
        else{
            calck2_ft(DIM1, DIM2, DIM3, DIM4, k_ft);
            calcdk2_ft(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_ft);
        }
        
        if (w_fb>t_fb){
            calck1_fb(DIM1, DIM2, DIM3, DIM4, k_fb);
            calcdk1_fb(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_fb);
        }
        else{
            calck2_fb(DIM1, DIM2, DIM3, DIM4, k_fb);
            calcdk2_fb(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_fb);
        }
        
        if (h_w>t_w){
            calck1_w(DIM1, DIM2, DIM3, DIM4, k_w);
            calcdk1_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_w);
        }
        else{
            calck2_w(DIM1, DIM2, DIM3, DIM4, k_w);
            calcdk2_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_w);
        }
        dc = 0.0;
    }
    
    
    if (w_ft>t_ft){
        calcJ1_ft(DIM1, DIM2, DIM3, DIM4, J_ft);
        calcdJ1_ft(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_ft);
    }
    else{
        calcJ2_ft(DIM1, DIM2, DIM3, DIM4, J_ft);
        calcdJ2_ft(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_ft);
    }
    
    if (w_fb>t_fb){
        calcJ1_fb(DIM1, DIM2, DIM3, DIM4, J_fb);
        calcdJ1_fb(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_fb);
    }
    else{
        calcJ2_fb(DIM1, DIM2, DIM3, DIM4, J_fb);
        calcdJ2_fb(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_fb);
    }
    
    if (h_w>t_w){
        calcJ1_w(DIM1, DIM2, DIM3, DIM4, J_w);
        calcdJ1_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_w);
    }
    else{
        calcJ2_w(DIM1, DIM2, DIM3, DIM4, J_w);
        calcdJ2_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_w);
    }
    
    dJ = 2.0*(dk_fb*J_fb + k_fb*dJ_fb) + 2.0*(dk_w*J_w + k_w*dJ_w) + (dk_ft*J_ft + k_ft*dJ_ft) + dc;
}


void MAST::Solid1DHatSectionElementPropertyCard::init() {
    
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
    else if (DIM2v>=(0.5*DIM1v)){
        libmesh_error_msg("DIM2>=0.5*DIM1");
    }

    _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DHatSectionProperty::calcA,
                                                              MAST::Solid1DHatSectionProperty::calcdA,
                                                              DIM1, DIM2, DIM3, 
                                                              DIM4));
    
    _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DHatSectionProperty::calcA,
                                                                MAST::Solid1DHatSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DHatSectionProperty::calcA,
                                                                MAST::Solid1DHatSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DHatSectionProperty::calcJ,
                                                                MAST::Solid1DHatSectionProperty::calcdJ,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4));
    
    _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DHatSectionProperty::calcIp,
                                                                MAST::Solid1DHatSectionProperty::calcdIp,
                                                                MAST::Solid1DHatSectionProperty::calcA,
                                                                MAST::Solid1DHatSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DHatSectionProperty::calcIz,
                                                                MAST::Solid1DHatSectionProperty::calcdIz,
                                                                MAST::Solid1DHatSectionProperty::calcIy,
                                                                MAST::Solid1DHatSectionProperty::calcdIy,
                                                                MAST::Solid1DHatSectionProperty::calcA,
                                                                MAST::Solid1DHatSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}
