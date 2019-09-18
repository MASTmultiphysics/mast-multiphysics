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
#include "property_cards/solid_1d_gbox_section_element_property_card.h"


// GBOX Section (GBOX in Astros 21.2, looks like Roman Numeral II)
void MAST::Solid1DGBOXSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& A){
    A = DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4;
}

void MAST::Solid1DGBOXSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dA){
    dA = DIM5*dDIM2*2.0+dDIM3*(DIM1-DIM5*2.0)+dDIM4*(DIM1-DIM5*2.0)-dDIM5*(DIM2*-2.0+DIM3*2.0+DIM4*2.0)+dDIM1*(DIM3+DIM4);
}

void MAST::Solid1DGBOXSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iz){
    Iz = (DIM1*(DIM3*DIM3*DIM3))/1.2E+1+(DIM1*(DIM4*DIM4*DIM4))/1.2E+1-(DIM5*pow(-DIM2+DIM3+DIM4,3.0))/6.0+DIM1*DIM3*pow(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4),2.0)+DIM1*DIM4*pow(DIM2/2.0-DIM3/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4),2.0)-DIM5*pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*(-DIM2+DIM3+DIM4)*2.0;
}

void MAST::Solid1DGBOXSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIz){
    dIz = (dDIM1*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1)*(DIM4*DIM4*DIM4*DIM4*DIM4)+(DIM3*DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*4.0+(DIM4*DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*1.0E+1+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*1.0E+1+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*(DIM5*DIM5)*2.4E+1+(DIM2*DIM2)*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*(DIM4*DIM4)*(DIM5*DIM5)*2.4E+1+(DIM3*DIM3)*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0E+1+(DIM3*DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*4.0E+1-DIM1*(DIM3*DIM3*DIM3*DIM3*DIM3)*DIM5*4.0-DIM1*(DIM4*DIM4*DIM4*DIM4*DIM4)*DIM5*4.0+(DIM1*DIM1)*DIM3*(DIM4*DIM4*DIM4*DIM4)*5.0+(DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)*DIM4*5.0-DIM2*(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*8.0+(DIM2*DIM2*DIM2*DIM2)*DIM3*(DIM5*DIM5)*1.2E+1-DIM2*(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*8.0+(DIM2*DIM2*DIM2*DIM2)*DIM4*(DIM5*DIM5)*1.2E+1+DIM3*(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*2.0E+1+(DIM3*DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*2.0E+1-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4*DIM4)*1.2E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*DIM4*1.2E+1-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*DIM5*4.0E+1-DIM1*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*DIM5*4.0E+1-DIM2*DIM3*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*8.0E+1-DIM2*(DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*8.0E+1-(DIM2*DIM2*DIM2)*DIM3*DIM4*(DIM5*DIM5)*9.6E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*(DIM4*DIM4)*2.4E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.2E+1+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*1.2E+1-DIM2*(DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*1.44E+2+(DIM2*DIM2)*DIM3*(DIM4*DIM4)*(DIM5*DIM5)*1.44E+2+(DIM2*DIM2)*(DIM3*DIM3)*DIM4*(DIM5*DIM5)*1.44E+2+DIM1*DIM2*(DIM3*DIM3*DIM3*DIM3)*DIM5*4.0+DIM1*DIM2*(DIM4*DIM4*DIM4*DIM4)*DIM5*4.0-DIM1*DIM3*(DIM4*DIM4*DIM4*DIM4)*DIM5*2.0E+1-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*DIM5*2.0E+1+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*DIM5*6.4E+1+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*DIM5*6.4E+1+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*DIM5*4.8E+1+DIM1*DIM2*(DIM3*DIM3)*(DIM4*DIM4)*DIM5*1.2E+2-DIM1*(DIM2*DIM2)*DIM3*(DIM4*DIM4)*DIM5*9.6E+1-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*DIM5*9.6E+1))/1.2E+1+(dDIM2*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*((DIM2*DIM2*DIM2*DIM2)*(DIM5*DIM5*DIM5)*4.0+(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5)*4.0+(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5)*4.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*4.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM5*DIM5*DIM5)*2.4E+1+(DIM2*DIM2)*(DIM4*DIM4)*(DIM5*DIM5*DIM5)*2.4E+1+(DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5*DIM5)*2.4E+1-(DIM1*DIM1*DIM1)*DIM3*(DIM4*DIM4*DIM4)*2.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3*DIM3)*DIM4*2.0-DIM1*(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*4.0+(DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)*DIM5-DIM1*(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0-DIM2*(DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5)*1.6E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4*DIM4)*DIM5-(DIM2*DIM2*DIM2)*DIM3*(DIM5*DIM5*DIM5)*1.6E+1-DIM2*(DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*DIM4*(DIM5*DIM5*DIM5)*1.6E+1+DIM3*(DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5)*1.6E+1+(DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5*DIM5)*1.6E+1+(DIM1*DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*4.0+(DIM1*DIM1*DIM1)*DIM2*(DIM3*DIM3)*DIM4*4.0+DIM1*DIM2*(DIM3*DIM3*DIM3)*(DIM5*DIM5)*1.6E+1+DIM1*(DIM2*DIM2*DIM2)*DIM3*(DIM5*DIM5)*8.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*DIM5*4.0+DIM1*DIM2*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.6E+1+DIM1*(DIM2*DIM2*DIM2)*DIM4*(DIM5*DIM5)*8.0-(DIM1*DIM1)*DIM2*(DIM4*DIM4*DIM4)*DIM5*4.0-DIM1*DIM3*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*2.4E+1-DIM1*(DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*2.4E+1+(DIM1*DIM1)*DIM3*(DIM4*DIM4*DIM4)*DIM5*1.2E+1+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*DIM4*DIM5*1.2E+1-DIM2*DIM3*(DIM4*DIM4)*(DIM5*DIM5*DIM5)*4.8E+1-DIM2*(DIM3*DIM3)*DIM4*(DIM5*DIM5*DIM5)*4.8E+1+(DIM2*DIM2)*DIM3*DIM4*(DIM5*DIM5*DIM5)*4.8E+1-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*(DIM5*DIM5)*2.0E+1+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*DIM5*4.0-DIM1*(DIM2*DIM2)*(DIM4*DIM4)*(DIM5*DIM5)*2.0E+1+(DIM1*DIM1)*(DIM2*DIM2)*(DIM4*DIM4)*DIM5*4.0-DIM1*(DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*4.0E+1+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*DIM5*2.2E+1+DIM1*DIM2*DIM3*(DIM4*DIM4)*(DIM5*DIM5)*6.4E+1+DIM1*DIM2*(DIM3*DIM3)*DIM4*(DIM5*DIM5)*6.4E+1-DIM1*(DIM2*DIM2)*DIM3*DIM4*(DIM5*DIM5)*4.8E+1-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*DIM5*2.8E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*DIM4*DIM5*2.8E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*DIM4*DIM5*1.2E+1))/2.0-(dDIM5*(-DIM2+DIM3+DIM4)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1)*(DIM4*DIM4*DIM4*DIM4)+(DIM2*DIM2*DIM2*DIM2)*(DIM5*DIM5)*4.0+(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*4.0+(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*4.0+(DIM1*DIM1)*(DIM2*DIM2)*(DIM4*DIM4)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4)*6.0+(DIM2*DIM2)*(DIM3*DIM3)*(DIM5*DIM5)*2.4E+1+(DIM2*DIM2)*(DIM4*DIM4)*(DIM5*DIM5)*2.4E+1+(DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*2.4E+1-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM5*4.0-DIM1*(DIM4*DIM4*DIM4*DIM4)*DIM5*4.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*2.0-(DIM1*DIM1)*DIM2*(DIM4*DIM4*DIM4)*2.0+(DIM1*DIM1)*DIM3*(DIM4*DIM4*DIM4)*4.0+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*DIM4*4.0-DIM2*(DIM3*DIM3*DIM3)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*DIM3*(DIM5*DIM5)*1.6E+1-DIM2*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*DIM4*(DIM5*DIM5)*1.6E+1+DIM3*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.6E+1+(DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*1.6E+1-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4)*6.0-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*DIM4*6.0-(DIM1*DIM1)*(DIM2*DIM2)*DIM3*DIM4*4.0-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM5*1.2E+1-DIM1*(DIM2*DIM2)*(DIM4*DIM4)*DIM5*1.2E+1-DIM1*(DIM3*DIM3)*(DIM4*DIM4)*DIM5*2.4E+1-DIM2*DIM3*(DIM4*DIM4)*(DIM5*DIM5)*4.8E+1-DIM2*(DIM3*DIM3)*DIM4*(DIM5*DIM5)*4.8E+1+(DIM2*DIM2)*DIM3*DIM4*(DIM5*DIM5)*4.8E+1+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM5*1.2E+1+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM5*4.0+DIM1*DIM2*(DIM4*DIM4*DIM4)*DIM5*1.2E+1+DIM1*(DIM2*DIM2*DIM2)*DIM4*DIM5*4.0-DIM1*DIM3*(DIM4*DIM4*DIM4)*DIM5*1.6E+1-DIM1*(DIM3*DIM3*DIM3)*DIM4*DIM5*1.6E+1+DIM1*DIM2*DIM3*(DIM4*DIM4)*DIM5*3.6E+1+DIM1*DIM2*(DIM3*DIM3)*DIM4*DIM5*3.6E+1-DIM1*(DIM2*DIM2)*DIM3*DIM4*DIM5*2.4E+1))/6.0+(dDIM3*(DIM1-DIM5*2.0)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*pow(DIM1*(DIM3*DIM3)+DIM1*(DIM4*DIM4)-(DIM2*DIM2)*DIM5*2.0-(DIM3*DIM3)*DIM5*2.0-(DIM4*DIM4)*DIM5*2.0-DIM1*DIM2*DIM4*2.0+DIM1*DIM3*DIM4*2.0+DIM2*DIM3*DIM5*4.0+DIM2*DIM4*DIM5*4.0-DIM3*DIM4*DIM5*4.0,2.0))/4.0+(dDIM4*(DIM1-DIM5*2.0)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*pow(DIM1*(DIM3*DIM3)+DIM1*(DIM4*DIM4)-(DIM2*DIM2)*DIM5*2.0-(DIM3*DIM3)*DIM5*2.0-(DIM4*DIM4)*DIM5*2.0-DIM1*DIM2*DIM3*2.0+DIM1*DIM3*DIM4*2.0+DIM2*DIM3*DIM5*4.0+DIM2*DIM4*DIM5*4.0-DIM3*DIM4*DIM5*4.0,2.0))/4.0;
}

void MAST::Solid1DGBOXSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM3)/1.2E+1+((DIM1*DIM1*DIM1)*DIM4)/1.2E+1-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0;
}

void MAST::Solid1DGBOXSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIy){
    dIy = dDIM3*((DIM1*DIM1*DIM1)/1.2E+1-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0)+dDIM4*((DIM1*DIM1*DIM1)/1.2E+1-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0)+((DIM1*DIM1)*dDIM1*(DIM3+DIM4))/4.0+(DIM5*dDIM2*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0-(dDIM5*pow(DIM5*2.0+DIM6,2.0)*(-DIM2+DIM3+DIM4))/2.0-(DIM5*dDIM6*(DIM5*6.0+DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0;
}

void MAST::Solid1DGBOXSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Ip){
    Ip = (DIM1*(DIM3*DIM3*DIM3))/1.2E+1+((DIM1*DIM1*DIM1)*DIM3)/1.2E+1+(DIM1*(DIM4*DIM4*DIM4))/1.2E+1+((DIM1*DIM1*DIM1)*DIM4)/1.2E+1-(DIM5*pow(-DIM2+DIM3+DIM4,3.0))/6.0-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0+DIM1*DIM3*pow(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4),2.0)+DIM1*DIM4*pow(DIM2/2.0-DIM3/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4),2.0)-DIM5*pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*(-DIM2+DIM3+DIM4)*2.0;
}

void MAST::Solid1DGBOXSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dIp){
    dIp = dDIM2*((DIM5*pow(-DIM2+DIM3+DIM4,2.0))/2.0+(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0+DIM5*pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*2.0+(DIM5*DIM5)*pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,3.0)*(-DIM2+DIM3+DIM4)*8.0+DIM1*DIM4*(DIM2/2.0-DIM3/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4))*(((DIM1*DIM3)/2.0-(DIM1*DIM4)/2.0)/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4)-DIM5*(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*2.0+1.0/2.0)*2.0-DIM1*DIM3*(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4))*(-((DIM1*DIM3)/2.0-(DIM1*DIM4)/2.0)/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4)+DIM5*(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*2.0+1.0/2.0)*2.0-DIM5*((DIM1*DIM3)/2.0-(DIM1*DIM4)/2.0)*(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*(-DIM2+DIM3+DIM4)*4.0)+dDIM4*((DIM1*(DIM4*DIM4))/4.0-(DIM5*pow(-DIM2+DIM3+DIM4,2.0))/2.0+(DIM1*DIM1*DIM1)/1.2E+1+(DIM1*pow(DIM2-DIM3+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0),2.0))/4.0-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0-((DIM1*DIM1)*(DIM2*DIM2)*DIM5*pow(DIM3-DIM4,2.0)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0))/2.0+(DIM1*DIM1)*(DIM2*DIM2)*DIM5*(DIM3-DIM4)*(-DIM2+DIM3+DIM4)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)-(DIM1*DIM1)*DIM2*DIM4*(DIM1*DIM3+DIM2*DIM5-DIM3*DIM5*2.0)*(DIM2/2.0-DIM3/2.0+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3*2.0+DIM1*DIM4*2.0+DIM2*DIM5*4.0-DIM3*DIM5*4.0-DIM4*DIM5*4.0))*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*2.0+DIM1*DIM3*(DIM1-DIM5*2.0)*(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3*2.0+DIM1*DIM4*2.0+DIM2*DIM5*4.0-DIM3*DIM5*4.0-DIM4*DIM5*4.0))*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*(DIM1*(DIM3*DIM3)+DIM1*(DIM4*DIM4)-(DIM2*DIM2)*DIM5*2.0-(DIM3*DIM3)*DIM5*2.0-(DIM4*DIM4)*DIM5*2.0-DIM1*DIM2*DIM3*2.0+DIM1*DIM3*DIM4*2.0+DIM2*DIM3*DIM5*4.0+DIM2*DIM4*DIM5*4.0-DIM3*DIM4*DIM5*4.0)+(DIM1*DIM1)*(DIM2*DIM2)*DIM5*(DIM1-DIM5*2.0)*pow(DIM3-DIM4,2.0)*(-DIM2+DIM3+DIM4)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,3.0))+dDIM3*((DIM1*(DIM3*DIM3))/4.0-(DIM5*pow(-DIM2+DIM3+DIM4,2.0))/2.0+(DIM1*DIM1*DIM1)/1.2E+1+(DIM1*pow(-DIM2+DIM4+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0),2.0))/4.0-(DIM5*((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0))/6.0-((DIM1*DIM1)*(DIM2*DIM2)*DIM5*pow(DIM3-DIM4,2.0)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0))/2.0-(DIM1*DIM1)*(DIM2*DIM2)*DIM5*(DIM3-DIM4)*(-DIM2+DIM3+DIM4)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)+(DIM1*DIM1)*DIM2*DIM3*(DIM1*DIM4+DIM2*DIM5-DIM4*DIM5*2.0)*(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3*2.0+DIM1*DIM4*2.0+DIM2*DIM5*4.0-DIM3*DIM5*4.0-DIM4*DIM5*4.0))*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*2.0-DIM1*DIM4*(DIM1-DIM5*2.0)*(DIM2/2.0-DIM3/2.0+(DIM1*DIM2*(DIM3-DIM4))/(DIM1*DIM3*2.0+DIM1*DIM4*2.0+DIM2*DIM5*4.0-DIM3*DIM5*4.0-DIM4*DIM5*4.0))*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*(DIM1*(DIM3*DIM3)+DIM1*(DIM4*DIM4)-(DIM2*DIM2)*DIM5*2.0-(DIM3*DIM3)*DIM5*2.0-(DIM4*DIM4)*DIM5*2.0-DIM1*DIM2*DIM4*2.0+DIM1*DIM3*DIM4*2.0+DIM2*DIM3*DIM5*4.0+DIM2*DIM4*DIM5*4.0-DIM3*DIM4*DIM5*4.0)+(DIM1*DIM1)*(DIM2*DIM2)*DIM5*(DIM1-DIM5*2.0)*pow(DIM3-DIM4,2.0)*(-DIM2+DIM3+DIM4)*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,3.0))-dDIM5*((((DIM5*DIM5)*4.0+(DIM6*DIM6)*3.0+DIM5*DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0+pow(-DIM2+DIM3+DIM4,3.0)/6.0+pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*(-DIM2+DIM3+DIM4)*2.0+(DIM5*(DIM5*8.0+DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0+DIM5*pow(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0),2.0)*(DIM2*-2.0+DIM3*2.0+DIM4*2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,3.0)*(-DIM2+DIM3+DIM4)*4.0-DIM1*DIM3*(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))*(DIM2*(-1.0/2.0)+DIM4/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4))*(DIM2*-2.0+DIM3*2.0+DIM4*2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*2.0-DIM1*DIM4*(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))*(DIM2/2.0-DIM3/2.0+(DIM1*DIM3*(DIM2/2.0-DIM4/2.0)-DIM1*DIM4*(DIM2/2.0-DIM3/2.0))/(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4))*(DIM2*-2.0+DIM3*2.0+DIM4*2.0)*1.0/pow(DIM5*(-DIM2+DIM3+DIM4)*-2.0+DIM1*DIM3+DIM1*DIM4,2.0)*2.0)+(dDIM1*1.0/pow(DIM1*DIM3+DIM1*DIM4+DIM2*DIM5*2.0-DIM3*DIM5*2.0-DIM4*DIM5*2.0,2.0)*((DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3*DIM3)+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3*DIM3)*3.0+(DIM1*DIM1)*(DIM4*DIM4*DIM4*DIM4*DIM4)+(DIM1*DIM1*DIM1*DIM1)*(DIM4*DIM4*DIM4)*3.0+(DIM3*DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*4.0+(DIM4*DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0+(DIM1*DIM1)*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*1.0E+1+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*1.0E+1+(DIM1*DIM1)*(DIM3*DIM3*DIM3)*(DIM5*DIM5)*1.2E+1+(DIM1*DIM1)*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.2E+1+(DIM2*DIM2)*(DIM3*DIM3*DIM3)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*(DIM3*DIM3)*(DIM5*DIM5)*2.4E+1+(DIM2*DIM2)*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*1.6E+1-(DIM2*DIM2*DIM2)*(DIM4*DIM4)*(DIM5*DIM5)*2.4E+1+(DIM3*DIM3)*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*4.0E+1+(DIM3*DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*4.0E+1-DIM1*(DIM3*DIM3*DIM3*DIM3*DIM3)*DIM5*4.0-DIM1*(DIM4*DIM4*DIM4*DIM4*DIM4)*DIM5*4.0+(DIM1*DIM1)*DIM3*(DIM4*DIM4*DIM4*DIM4)*5.0+(DIM1*DIM1)*(DIM3*DIM3*DIM3*DIM3)*DIM4*5.0+(DIM1*DIM1*DIM1*DIM1)*DIM3*(DIM4*DIM4)*9.0+(DIM1*DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*9.0-(DIM1*DIM1*DIM1)*(DIM3*DIM3*DIM3)*DIM5*1.2E+1-DIM2*(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*8.0-(DIM1*DIM1*DIM1)*(DIM4*DIM4*DIM4)*DIM5*1.2E+1+(DIM2*DIM2*DIM2*DIM2)*DIM3*(DIM5*DIM5)*1.2E+1-DIM2*(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*8.0+(DIM2*DIM2*DIM2*DIM2)*DIM4*(DIM5*DIM5)*1.2E+1+DIM3*(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*2.0E+1+(DIM3*DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*2.0E+1-(DIM1*DIM1)*DIM2*DIM3*(DIM4*DIM4*DIM4)*1.2E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3*DIM3)*DIM4*1.2E+1+(DIM1*DIM1*DIM1)*DIM2*(DIM3*DIM3)*DIM5*1.2E+1+(DIM1*DIM1*DIM1)*DIM2*(DIM4*DIM4)*DIM5*1.2E+1-DIM1*(DIM3*DIM3)*(DIM4*DIM4*DIM4)*DIM5*4.0E+1-DIM1*(DIM3*DIM3*DIM3)*(DIM4*DIM4)*DIM5*4.0E+1-(DIM1*DIM1*DIM1)*DIM3*(DIM4*DIM4)*DIM5*3.6E+1-(DIM1*DIM1*DIM1)*(DIM3*DIM3)*DIM4*DIM5*3.6E+1-DIM2*DIM3*(DIM4*DIM4*DIM4)*(DIM5*DIM5)*8.0E+1-DIM2*(DIM3*DIM3*DIM3)*DIM4*(DIM5*DIM5)*8.0E+1-(DIM2*DIM2*DIM2)*DIM3*DIM4*(DIM5*DIM5)*9.6E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*(DIM4*DIM4)*2.4E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*(DIM4*DIM4)*1.2E+1+(DIM1*DIM1)*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*1.2E+1-(DIM1*DIM1)*DIM2*(DIM3*DIM3)*(DIM5*DIM5)*2.4E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM3*(DIM5*DIM5)*1.2E+1-(DIM1*DIM1)*DIM2*(DIM4*DIM4)*(DIM5*DIM5)*2.4E+1+(DIM1*DIM1)*(DIM2*DIM2)*DIM4*(DIM5*DIM5)*1.2E+1+(DIM1*DIM1)*DIM3*(DIM4*DIM4)*(DIM5*DIM5)*3.6E+1+(DIM1*DIM1)*(DIM3*DIM3)*DIM4*(DIM5*DIM5)*3.6E+1-DIM2*(DIM3*DIM3)*(DIM4*DIM4)*(DIM5*DIM5)*1.44E+2+(DIM2*DIM2)*DIM3*(DIM4*DIM4)*(DIM5*DIM5)*1.44E+2+(DIM2*DIM2)*(DIM3*DIM3)*DIM4*(DIM5*DIM5)*1.44E+2+DIM1*DIM2*(DIM3*DIM3*DIM3*DIM3)*DIM5*4.0+DIM1*DIM2*(DIM4*DIM4*DIM4*DIM4)*DIM5*4.0-DIM1*DIM3*(DIM4*DIM4*DIM4*DIM4)*DIM5*2.0E+1-DIM1*(DIM3*DIM3*DIM3*DIM3)*DIM4*DIM5*2.0E+1+DIM1*DIM2*DIM3*(DIM4*DIM4*DIM4)*DIM5*6.4E+1+DIM1*DIM2*(DIM3*DIM3*DIM3)*DIM4*DIM5*6.4E+1+DIM1*(DIM2*DIM2*DIM2)*DIM3*DIM4*DIM5*4.8E+1+(DIM1*DIM1*DIM1)*DIM2*DIM3*DIM4*DIM5*2.4E+1+DIM1*DIM2*(DIM3*DIM3)*(DIM4*DIM4)*DIM5*1.2E+2-DIM1*(DIM2*DIM2)*DIM3*(DIM4*DIM4)*DIM5*9.6E+1-DIM1*(DIM2*DIM2)*(DIM3*DIM3)*DIM4*DIM5*9.6E+1-(DIM1*DIM1)*DIM2*DIM3*DIM4*(DIM5*DIM5)*4.8E+1))/1.2E+1-(DIM5*dDIM6*(DIM5*6.0+DIM6*6.0)*(-DIM2+DIM3+DIM4))/6.0;
}

void MAST::Solid1DGBOXSectionProperty::calcJ_box(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J_box){
    J_box = (DIM5*(DIM3/2.0+DIM4/2.0)*pow(DIM5+DIM6,2.0)*pow(-DIM2+DIM3/2.0+DIM4/2.0,2.0)*-2.0)/(pow(DIM3/2.0+DIM4/2.0,2.0)-DIM2*(DIM3/2.0+DIM4/2.0)+DIM5*DIM5-DIM5*(DIM5*2.0+DIM6));
}

void MAST::Solid1DGBOXSectionProperty::calcdJ_box(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ_box){
    dJ_box = DIM5*dDIM3*pow(DIM5+DIM6,2.0)*(DIM2*-2.0+DIM3+DIM4)*1.0/pow(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*4.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*4.0,2.0)*(DIM2*(DIM3*DIM3)*2.0+DIM2*(DIM4*DIM4)*2.0-DIM2*(DIM5*DIM5)*8.0-DIM3*(DIM4*DIM4)*3.0-(DIM3*DIM3)*DIM4*3.0+DIM3*(DIM5*DIM5)*1.2E+1+DIM4*(DIM5*DIM5)*1.2E+1-DIM3*DIM3*DIM3-DIM4*DIM4*DIM4+DIM2*DIM3*DIM4*4.0-DIM2*DIM5*DIM6*8.0+DIM3*DIM5*DIM6*1.2E+1+DIM4*DIM5*DIM6*1.2E+1)+DIM5*dDIM4*pow(DIM5+DIM6,2.0)*(DIM2*-2.0+DIM3+DIM4)*1.0/pow(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*4.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*4.0,2.0)*(DIM2*(DIM3*DIM3)*2.0+DIM2*(DIM4*DIM4)*2.0-DIM2*(DIM5*DIM5)*8.0-DIM3*(DIM4*DIM4)*3.0-(DIM3*DIM3)*DIM4*3.0+DIM3*(DIM5*DIM5)*1.2E+1+DIM4*(DIM5*DIM5)*1.2E+1-DIM3*DIM3*DIM3-DIM4*DIM4*DIM4+DIM2*DIM3*DIM4*4.0-DIM2*DIM5*DIM6*8.0+DIM3*DIM5*DIM6*1.2E+1+DIM4*DIM5*DIM6*1.2E+1)-dDIM5*(DIM3+DIM4)*(DIM5+DIM6)*pow(DIM2*-2.0+DIM3+DIM4,2.0)*1.0/pow(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*4.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*4.0,2.0)*((DIM3*DIM3)*DIM5*3.0+(DIM3*DIM3)*DIM6+(DIM4*DIM4)*DIM5*3.0+(DIM4*DIM4)*DIM6-(DIM5*DIM5)*DIM6*4.0-(DIM5*DIM5*DIM5)*4.0-DIM2*DIM3*DIM5*6.0-DIM2*DIM3*DIM6*2.0-DIM2*DIM4*DIM5*6.0-DIM2*DIM4*DIM6*2.0+DIM3*DIM4*DIM5*6.0+DIM3*DIM4*DIM6*2.0)+DIM5*dDIM6*(DIM3+DIM4)*(DIM5+DIM6)*pow(DIM2*-2.0+DIM3+DIM4,2.0)*(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*2.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*2.0)*1.0/pow(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*4.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*4.0,2.0)*2.0-DIM5*dDIM2*(DIM3+DIM4)*pow(DIM5+DIM6,2.0)*(DIM2*-2.0+DIM3+DIM4)*1.0/pow(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*4.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*4.0,2.0)*(-DIM3*DIM3-DIM4*DIM4+(DIM5*DIM5)*8.0+DIM2*DIM3*2.0+DIM2*DIM4*2.0-DIM3*DIM4*2.0+DIM5*DIM6*8.0)*2.0;
}

void MAST::Solid1DGBOXSectionProperty::calcJ1_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_fendst){
    J1_fendst = -(DIM3*DIM3*DIM3)*(DIM1*(-1.0/2.0)+DIM5+DIM6/2.0);
}

void MAST::Solid1DGBOXSectionProperty::calcdJ1_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_fendst){
    dJ1_fendst = ((DIM3*DIM3*DIM3)*dDIM1)/2.0-(DIM3*DIM3*DIM3)*dDIM5-((DIM3*DIM3*DIM3)*dDIM6)/2.0-(DIM3*DIM3)*dDIM3*(DIM1*(-1.0/2.0)+DIM5+DIM6/2.0)*3.0;
}

void MAST::Solid1DGBOXSectionProperty::calcJ2_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_fendst){
    J2_fendst = DIM3*pow(-DIM1+DIM5*2.0+DIM6,3.0)*(-1.0/8.0);
}

void MAST::Solid1DGBOXSectionProperty::calcdJ2_fendst(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_fendst){
    dJ2_fendst = pow(-DIM1+DIM5*2.0+DIM6,2.0)*(-DIM1*dDIM3-DIM3*dDIM1*3.0+DIM3*dDIM5*6.0+DIM5*dDIM3*2.0+DIM3*dDIM6*3.0+DIM6*dDIM3)*(-1.0/8.0);
}

void MAST::Solid1DGBOXSectionProperty::calcJ1_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J1_fendsb){
    J1_fendsb = -(DIM4*DIM4*DIM4)*(DIM1*(-1.0/2.0)+DIM5+DIM6/2.0);
}

void MAST::Solid1DGBOXSectionProperty::calcdJ1_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ1_fendsb){
    dJ1_fendsb = ((DIM4*DIM4*DIM4)*dDIM1)/2.0-(DIM4*DIM4*DIM4)*dDIM5-((DIM4*DIM4*DIM4)*dDIM6)/2.0-(DIM4*DIM4)*dDIM4*(DIM1*(-1.0/2.0)+DIM5+DIM6/2.0)*3.0;
}

void MAST::Solid1DGBOXSectionProperty::calcJ2_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J2_fendsb){
    J2_fendsb = DIM4*pow(-DIM1+DIM5*2.0+DIM6,3.0)*(-1.0/8.0);
}

void MAST::Solid1DGBOXSectionProperty::calcdJ2_fendsb(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ2_fendsb){
    dJ2_fendsb = pow(-DIM1+DIM5*2.0+DIM6,2.0)*(-DIM1*dDIM4-DIM4*dDIM1*3.0+DIM4*dDIM5*6.0+DIM5*dDIM4*2.0+DIM4*dDIM6*6.0+DIM6*dDIM4)*(-1.0/8.0);
}

void MAST::Solid1DGBOXSectionProperty::calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& Jc){
    Jc = (DIM3*DIM3*DIM3*DIM3)*(-2.1E+1/1.0E+2)-(DIM4*DIM4*DIM4*DIM4)*(2.1E+1/1.0E+2)-(1.0/(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*pow(DIM3*DIM3+(DIM5*DIM5)/4.0,4.0)*((DIM3*DIM3)*4.2E+2+(DIM5*DIM5)*7.25E+2-DIM3*DIM5*2.204E+3))/5.0E+3-(1.0/(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*pow(DIM4*DIM4+(DIM5*DIM5)/4.0,4.0)*((DIM4*DIM4)*4.2E+2+(DIM5*DIM5)*7.25E+2-DIM4*DIM5*2.204E+3))/5.0E+3;
}

void MAST::Solid1DGBOXSectionProperty::calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJc){
    dJc = dDIM5*((1.0/(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*pow(DIM3*DIM3+(DIM5*DIM5)/4.0,4.0)*(DIM3*2.204E+3-DIM5*1.45E+3))/5.0E+3+(1.0/(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*pow(DIM4*DIM4+(DIM5*DIM5)/4.0,4.0)*(DIM4*2.204E+3-DIM5*1.45E+3))/5.0E+3-(1.0/(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*DIM5*pow(DIM3*DIM3+(DIM5*DIM5)/4.0,3.0)*((DIM3*DIM3)*4.2E+2+(DIM5*DIM5)*7.25E+2-DIM3*DIM5*2.204E+3))/2.5E+3-(1.0/(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*DIM5*pow(DIM4*DIM4+(DIM5*DIM5)/4.0,3.0)*((DIM4*DIM4)*4.2E+2+(DIM5*DIM5)*7.25E+2-DIM4*DIM5*2.204E+3))/2.5E+3)-(1.0/(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*dDIM3*(DIM3*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*5.51E+3-(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*DIM5*8.46336E+5+pow(DIM3,1.0E+1)*7.5264E+5-pow(DIM5,1.0E+1)*2.175E+3-(DIM3*DIM3)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*2.404E+4+(DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*5.2896E+4-(DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*7.632E+4+(DIM3*DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5*DIM5*DIM5)*1.05792E+5-(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5*DIM5)*2.82112E+5+(DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3*DIM3)*(DIM5*DIM5)*2.9312E+5))/6.4E+5-(1.0/(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*dDIM4*(DIM4*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*5.51E+3-(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*DIM5*8.46336E+5+pow(DIM4,1.0E+1)*7.5264E+5-pow(DIM5,1.0E+1)*2.175E+3-(DIM4*DIM4)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*2.404E+4+(DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*5.2896E+4-(DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5*DIM5*DIM5*DIM5)*7.632E+4+(DIM4*DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5*DIM5*DIM5)*1.05792E+5-(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5*DIM5)*2.82112E+5+(DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4*DIM4)*(DIM5*DIM5)*2.9312E+5))/6.4E+5;
}



// TODO: Improve accuracy of calculation of torsion constant and corresponding derivative
void MAST::Solid1DGBOXSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& J)
{
    // NOTE: Not currently using El Darwish & Johnson correction.
    
    Real t_ft = DIM3;
    Real t_fb = DIM4;
    Real w_ft = DIM1;
    Real w_fb = DIM1;
    Real h_w = DIM2 - DIM3 - DIM4;
    Real t_w = DIM5;
    
    Real k_ft, k_fb, k_w, ct, cb, J_box, J_fb, J_ft;
    Real wf_t = t_w/t_ft;
    Real wf_b = t_w/t_fb;
    
    // torsion constant of the box section
    calcJ_box(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_box);
    
    if (w_ft>t_ft){
        calcJ1_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_ft);
    }
    else{
        calcJ2_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_ft);
    }
    
    if (w_fb>t_fb){
        calcJ1_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_fb);
    }
    else{
        calcJ2_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_fb);
    }
    
    J = J_box + 2.0*(k_ft*J_ft) + 2.0*(k_fb*J_fb);
}



void MAST::Solid1DGBOXSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& DIM5, Real& DIM6, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dDIM5, Real& dDIM6, Real& dJ)
{
    // NOTE: Not currently using El Darwish & Johnson correction.
    
    Real t_ft = DIM3;
    Real t_fb = DIM4;
    Real w_ft = DIM1;
    Real w_fb = DIM1;
    Real h_w = DIM2 - DIM3 - DIM4;
    Real t_w = DIM5;
    
    Real k_ft, k_fb, k_w, ct, cb, J_box, J_fb, J_ft;
    Real dk_ft, dk_fb, dk_w, dct, dcb, dJ_box, dJ_fb, dJ_ft;
    Real wf_t = t_w/t_ft;
    Real wf_b = t_w/t_fb;
    
    // torsion constant of the box section
    calcdJ_box(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6,
               dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dJ_box);
    
    if (w_ft>t_ft){
        calcJ1_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_ft);
        calcdJ1_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6,
                       dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dJ_ft);
    }
    else{
        calcJ2_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_ft);
        calcdJ2_fendst(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6,
                       dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dJ_ft);
    }
    
    if (w_fb>t_fb){
        calcJ1_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_fb);
        calcdJ1_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6,
                       dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dJ_fb);
    }
    else{
        calcJ2_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, J_fb);
        calcdJ2_fendsb(DIM1, DIM2, DIM3, DIM4, DIM5, DIM6,
                       dDIM1, dDIM2, dDIM3, dDIM4, dDIM5, dDIM6, dJ_fb);
    }
    
    
    //J = J_box + 2.0*(k_ft*J_ft) + 2.0*(k_fb*J_fb);
    dJ = dJ_box + 2.0*(dk_ft*J_ft + k_ft*dJ_ft) + 2.0*(dk_fb*J_fb + k_fb*dJ_fb);
}





void MAST::Solid1DGboxSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &DIM3     =  this->get<MAST::FieldFunction<Real> >("DIM3"),
    &DIM4     =  this->get<MAST::FieldFunction<Real> >("DIM4"),
    &DIM5     =  this->get<MAST::FieldFunction<Real> >("DIM5"),
    &DIM6     =  this->get<MAST::FieldFunction<Real> >("DIM6"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v, DIM2v, DIM3v, DIM4v, DIM5v, DIM6v;
    DIM1(DIM1v); DIM2(DIM2v); DIM3(DIM3v); DIM4(DIM4v); DIM5(DIM5v); DIM6(DIM6v);
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
    else if (DIM5v<=0){
        libmesh_error_msg("DIM5<=0");
    }
    else if (DIM6v<=0){
        libmesh_error_msg("DIM6<=0");
    }
    else if ((DIM6v+2.0*DIM5v)>=DIM1v){
        libmesh_error_msg("(DIM6+2*DIM5)>=DIM1");
    }
    else if ((DIM3v+DIM4v)>=DIM2v){
        libmesh_error_msg("(DIM3+DIM4)>=DIM2");
    }
    
    
    _A.reset(new MAST::Solid1D6ParameterSectionProperty::Area(MAST::Solid1DGBOXSectionProperty::calcA,
                                                              MAST::Solid1DGBOXSectionProperty::calcdA,
                                                              DIM1, DIM2, DIM3, 
                                                              DIM4, DIM5, DIM6));
    
    _Ay.reset(new MAST::Solid1D6ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DGBOXSectionProperty::calcA,
                                                                MAST::Solid1DGBOXSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4, 
                                                                DIM5, DIM6,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D6ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DGBOXSectionProperty::calcA,
                                                                MAST::Solid1DGBOXSectionProperty::calcdA,
                                                                DIM1, DIM2, 
                                                                DIM3, DIM4,
                                                                DIM5, DIM6,
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D6ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DGBOXSectionProperty::calcJ,
                                                                MAST::Solid1DGBOXSectionProperty::calcdJ,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                DIM5, DIM6));
    
    _Ip.reset(new MAST::Solid1D6ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DGBOXSectionProperty::calcIp,
                                                                MAST::Solid1DGBOXSectionProperty::calcdIp,
                                                                MAST::Solid1DGBOXSectionProperty::calcA,
                                                                MAST::Solid1DGBOXSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                DIM5, DIM6,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D6ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DGBOXSectionProperty::calcIz,
                                                                MAST::Solid1DGBOXSectionProperty::calcdIz,
                                                                MAST::Solid1DGBOXSectionProperty::calcIy,
                                                                MAST::Solid1DGBOXSectionProperty::calcdIy,
                                                                MAST::Solid1DGBOXSectionProperty::calcA,
                                                                MAST::Solid1DGBOXSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                DIM3, DIM4,
                                                                DIM5, DIM6,
                                                                hy_off,
                                                                hz_off));
    
    _initialized = true;
}

