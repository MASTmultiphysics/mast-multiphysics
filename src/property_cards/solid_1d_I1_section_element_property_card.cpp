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
#include "property_cards/solid_1d_I1_section_element_property_card.h"


// Bi-Symmetric I_Beam (I1 in Nastran)
void MAST::Solid1DI1SectionProperty::calcA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& A){
    A = DIM4*(DIM1+DIM2)-DIM1*DIM3;
}

void MAST::Solid1DI1SectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dA){
    dA = -DIM1*dDIM3+DIM4*dDIM2-dDIM1*(DIM3-DIM4)+dDIM4*(DIM1+DIM2);
}

void MAST::Solid1DI1SectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iz){
    Iz = (DIM2*(DIM3*DIM3*DIM3))/1.2E+1-((DIM3*DIM3*DIM3-DIM4*DIM4*DIM4)*(DIM1+DIM2))/1.2E+1;
}

void MAST::Solid1DI1SectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIz){
    dIz = ((DIM4*DIM4*DIM4)*dDIM2)/1.2E+1-(dDIM1*(DIM3*DIM3*DIM3-DIM4*DIM4*DIM4))/1.2E+1+((DIM4*DIM4)*dDIM4*(DIM1+DIM2))/4.0-(DIM1*(DIM3*DIM3)*dDIM3)/4.0;
}

void MAST::Solid1DI1SectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Iy){
    Iy = ((DIM2*DIM2*DIM2)*DIM3)/1.2E+1-((DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,3.0))/6.0;
}

void MAST::Solid1DI1SectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIy){
    dIy = (dDIM4*pow(DIM1+DIM2,3.0))/1.2E+1+dDIM2*(((DIM2*DIM2)*DIM3)/4.0-((DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0))/2.0)-dDIM3*(pow(DIM1+DIM2,3.0)/1.2E+1-(DIM2*DIM2*DIM2)/1.2E+1)-(dDIM1*(DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0))/2.0;
}

void MAST::Solid1DI1SectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Ip){
    Ip = (DIM2*(DIM3*DIM3*DIM3))/1.2E+1+((DIM2*DIM2*DIM2)*DIM3)/1.2E+1-((DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,3.0))/6.0-((DIM3*DIM3*DIM3-DIM4*DIM4*DIM4)*(DIM1+DIM2))/1.2E+1;
}

void MAST::Solid1DI1SectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dIp){
    dIp = -dDIM1*(((DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0))/2.0+(DIM3*DIM3*DIM3)/1.2E+1-(DIM4*DIM4*DIM4)/1.2E+1)+dDIM4*(pow(DIM1+DIM2,3.0)/1.2E+1+((DIM4*DIM4)*(DIM1+DIM2))/4.0)+dDIM2*(((DIM2*DIM2)*DIM3)/4.0-((DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0))/2.0+(DIM4*DIM4*DIM4)/1.2E+1)-(DIM1*dDIM3*(DIM1*DIM1+(DIM2*DIM2)*3.0+(DIM3*DIM3)*3.0+DIM1*DIM2*3.0))/1.2E+1;
}

void MAST::Solid1DI1SectionProperty::calcJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_w){
    J1_w = (DIM2*DIM2*DIM2)*DIM3;
}

void MAST::Solid1DI1SectionProperty::calcdJ1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_w){
    dJ1_w = (DIM2*DIM2*DIM2)*dDIM3+(DIM2*DIM2)*DIM3*dDIM2*3.0;
}

void MAST::Solid1DI1SectionProperty::calcJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_w){
    J2_w = DIM2*(DIM3*DIM3*DIM3);
}

void MAST::Solid1DI1SectionProperty::calcdJ2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_w){
    dJ2_w = (DIM3*DIM3*DIM3)*dDIM2+DIM2*(DIM3*DIM3)*dDIM3*3.0;
}

void MAST::Solid1DI1SectionProperty::calcJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J1_f){
    J1_f = pow(DIM3-DIM4,3.0)*(DIM1+DIM2)*(-1.0/8.0);
}

void MAST::Solid1DI1SectionProperty::calcdJ1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ1_f){
    dJ1_f = dDIM1*pow(DIM3-DIM4,3.0)*(-1.0/8.0)-(dDIM2*pow(DIM3-DIM4,3.0))/8.0-dDIM3*pow(DIM3-DIM4,2.0)*(DIM1+DIM2)*(3.0/8.0)+dDIM4*pow(DIM3-DIM4,2.0)*(DIM1+DIM2)*(3.0/8.0);
}

void MAST::Solid1DI1SectionProperty::calcJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J2_f){
    J2_f = -(DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,3.0);
}

void MAST::Solid1DI1SectionProperty::calcdJ2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ2_f){
    dJ2_f = dDIM3*pow(DIM1+DIM2,3.0)*(-1.0/2.0)+(dDIM4*pow(DIM1+DIM2,3.0))/2.0-dDIM1*(DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0)*3.0-dDIM2*(DIM3/2.0-DIM4/2.0)*pow(DIM1+DIM2,2.0)*3.0;
}

void MAST::Solid1DI1SectionProperty::calck1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_f){
    k1_f = 1.0/3.0;
}

void MAST::Solid1DI1SectionProperty::calcdk1_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_f){
    dk1_f = 0.0;
}

void MAST::Solid1DI1SectionProperty::calck2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_f){
    k2_f = 1.0/3.0;
}

void MAST::Solid1DI1SectionProperty::calcdk2_f(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_f){
    dk2_f = 0.0;
}

void MAST::Solid1DI1SectionProperty::calck1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k1_w){
    k1_w = 1.0/3.0;
}

void MAST::Solid1DI1SectionProperty::calcdk1_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk1_w){
    dk1_w = 0.0;
}

void MAST::Solid1DI1SectionProperty::calck2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& k2_w){
    k2_w = 1.0/3.0;
}

void MAST::Solid1DI1SectionProperty::calcdk2_w(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dk2_w){
    dk2_w = 0.0;
}

void MAST::Solid1DI1SectionProperty::calcJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& Jc){
    Jc = pow(DIM3-DIM4,4.0)*(-2.1E+1/8.0E+2)-(1.0/pow(DIM3-DIM4,6.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,4.0)*((DIM2*DIM2)*7.25E+2+(DIM3*DIM3)*1.05E+2+(DIM4*DIM4)*1.05E+2+DIM2*DIM3*1.102E+3-DIM2*DIM4*1.102E+3-DIM3*DIM4*2.1E+2))/2.0E+4;
}

void MAST::Solid1DI1SectionProperty::calcdJc(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJc){
    dJc = -dDIM3*(pow(DIM3-DIM4,3.0)*(2.1E+1/2.0E+2)-1.0/pow(DIM3-DIM4,7.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,4.0)*((DIM2*DIM2)*7.25E+2+(DIM3*DIM3)*1.05E+2+(DIM4*DIM4)*1.05E+2+DIM2*DIM3*1.102E+3-DIM2*DIM4*1.102E+3-DIM3*DIM4*2.1E+2)*3.0E-4+(1.0/pow(DIM3-DIM4,6.0)*(DIM2*1.102E+3+DIM3*2.1E+2-DIM4*2.1E+2)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,4.0))/2.0E+4+(1.0/pow(DIM3-DIM4,6.0)*(DIM3*2.0-DIM4*2.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,3.0)*((DIM2*DIM2)*7.25E+2+(DIM3*DIM3)*1.05E+2+(DIM4*DIM4)*1.05E+2+DIM2*DIM3*1.102E+3-DIM2*DIM4*1.102E+3-DIM3*DIM4*2.1E+2))/5.0E+3)+dDIM4*(pow(DIM3-DIM4,3.0)*(2.1E+1/2.0E+2)-1.0/pow(DIM3-DIM4,7.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,4.0)*((DIM2*DIM2)*7.25E+2+(DIM3*DIM3)*1.05E+2+(DIM4*DIM4)*1.05E+2+DIM2*DIM3*1.102E+3-DIM2*DIM4*1.102E+3-DIM3*DIM4*2.1E+2)*3.0E-4+(1.0/pow(DIM3-DIM4,6.0)*(DIM2*1.102E+3+DIM3*2.1E+2-DIM4*2.1E+2)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,4.0))/2.0E+4+(1.0/pow(DIM3-DIM4,6.0)*(DIM3*2.0-DIM4*2.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,3.0)*((DIM2*DIM2)*7.25E+2+(DIM3*DIM3)*1.05E+2+(DIM4*DIM4)*1.05E+2+DIM2*DIM3*1.102E+3-DIM2*DIM4*1.102E+3-DIM3*DIM4*2.1E+2))/5.0E+3)-(dDIM2*1.0/pow(DIM3-DIM4,6.0)*pow(DIM2*DIM2+DIM3*DIM3+DIM4*DIM4-DIM3*DIM4*2.0,3.0)*(DIM2*(DIM3*DIM3)*1.145E+3+(DIM2*DIM2)*DIM3*4.959E+3+DIM2*(DIM4*DIM4)*1.145E+3-(DIM2*DIM2)*DIM4*4.959E+3+DIM3*(DIM4*DIM4)*1.653E+3-(DIM3*DIM3)*DIM4*1.653E+3+(DIM2*DIM2*DIM2)*3.625E+3+(DIM3*DIM3*DIM3)*5.51E+2-(DIM4*DIM4*DIM4)*5.51E+2-DIM2*DIM3*DIM4*2.29E+3))/1.0E+4;
}


// TODO: Improve accuracy of calculation of torsion constant and corresponding derivative
void MAST::Solid1DI1SectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& J)
{
    Real t_f = 0.5*(DIM4-DIM3);
    Real w_f = DIM2+DIM1;
    Real h_w = DIM3;
    Real t_w = DIM2;
    
    Real k_f, k_w, c, J_w, J_f;
    Real wf = t_w/t_f;
    
    // Sum of rectangles with correction of El Darwish and Johnson
    if ( (wf>0.5) and (wf<1.0) )
    {
        k_f = 0.3333333333333333;
        k_w = 0.3333333333333333;
        calcJc(DIM1, DIM2, DIM3, DIM4, c);
    }
    
    // Sum of rectangles method, with rectangles corrected for aspect ratios
    else
    {
        if (w_f>t_f){
            calck1_f(DIM1, DIM2, DIM3, DIM4, k_f);
        }
        else{
            calck2_f(DIM1, DIM2, DIM3, DIM4, k_f);
        }
        
        if (h_w>t_w){
            calck1_w(DIM1, DIM2, DIM3, DIM4, k_w);
        }
        else{
            calck2_w(DIM1, DIM2, DIM3, DIM4, k_w);
        }
        c = 0.0;
    }
    
    
    if (w_f>t_f){
        calcJ1_f(DIM1, DIM2, DIM3, DIM4, J_f);
    }
    else{
        calcJ2_f(DIM1, DIM2, DIM3, DIM4, J_f);
    }
    
    if (h_w>t_w){
        calcJ1_w(DIM1, DIM2, DIM3, DIM4, J_w);
    }
    else{
        calcJ2_w(DIM1, DIM2, DIM3, DIM4, J_w);
    }
    
    J = 2.0*(k_f*J_f) + (k_w*J_w) + c;
}



void MAST::Solid1DI1SectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& DIM3, Real& DIM4, Real& dDIM1, Real& dDIM2, Real& dDIM3, Real& dDIM4, Real& dJ)
{
    Real t_f = 0.5*(DIM4-DIM3);
    Real w_f = DIM2+DIM1;
    Real h_w = DIM3;
    Real t_w = DIM2;
    
    Real k_f, k_w, c, dk_f, dk_w, dc, J_w, J_f, dJ_w, dJ_f;
    Real wf = t_w/t_f;
    
    // Sum of rectangles with correction of El Darwish and Johnson
    if ( (wf>0.5) and (wf<1.0) )
    {
        k_f = 0.3333333333333333;
        k_w = 0.3333333333333333;
        dk_f = 0.0;
        dk_w = 0.0;
        calcdJc(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dc);
    }
    
    // Sum of rectangles method, with rectangles corrected for aspect ratios
    else
    {
        if (w_f>t_f){
            calck1_f(DIM1, DIM2, DIM3, DIM4, k_f);
            calcdk1_f(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_f);
        }
        else{
            calck2_f(DIM1, DIM2, DIM3, DIM4, k_f);
            calcdk2_f(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dk_f);
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
    
    
    if (w_f>t_f){
        calcJ1_f(DIM1, DIM2, DIM3, DIM4, J_f);
        calcdJ1_f(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_f);
    }
    else{
        calcJ2_f(DIM1, DIM2, DIM3, DIM4, J_f);
        calcdJ2_f(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_f);
    }
    
    if (h_w>t_w){
        calcJ1_w(DIM1, DIM2, DIM3, DIM4, J_w);
        calcdJ1_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_w);
    }
    else{
        calcJ2_w(DIM1, DIM2, DIM3, DIM4, J_w);
        calcdJ2_w(DIM1, DIM2, DIM3, DIM4, dDIM1, dDIM2, dDIM3, dDIM4, dJ_w);
    }
    
    dJ = 2.0*(dk_f*J_f + k_f*dJ_f) + (dk_w*J_w + k_w*dJ_w) + dc;
}


const std::vector<libMesh::Point> 
MAST::Solid1DI1SectionElementPropertyCard::get_geom_points(const libMesh::Point& p, const Real t, const uint n) const 
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
    
    libMesh::Point shift(0.5*(DIM1v+DIM2v), 0.5*DIM4v);
    libMesh::Point offset(offset_z, offset_y);
    
    Real th = 0.5*(DIM4v - DIM3v);
    std::vector<libMesh::Point> points = {
        libMesh::Point(0., 0.) - shift + offset,
        libMesh::Point(DIM1v+DIM2v, 0.) - shift + offset,
        libMesh::Point(DIM1v+DIM2v, th) - shift + offset,
        libMesh::Point(0.5*DIM1v+DIM2v, th) - shift + offset,
        libMesh::Point(0.5*DIM1v+DIM2v, th+DIM3v) - shift + offset,
        libMesh::Point(DIM1v+DIM2v, th+DIM3v) - shift + offset,
        libMesh::Point(DIM1v+DIM2v, DIM4v) - shift + offset,
        libMesh::Point(0., DIM4v) - shift + offset,
        libMesh::Point(0., th+DIM3v) - shift + offset,
        libMesh::Point(0.5*DIM1v, th+DIM3v) - shift + offset,
        libMesh::Point(0.5*DIM1v, th) - shift + offset,
        libMesh::Point(0., th) - shift + offset
    };
    
    return points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DI1SectionElementPropertyCard::get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const
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
    
    libMesh::Point dshift(0.5*(dDIM1v+dDIM2v), 0.5*dDIM4v);
    libMesh::Point doffset(doffset_z, doffset_y);
    
    Real dt = 0.5*(dDIM4v - dDIM3v);
    std::vector<libMesh::Point> points = {
        libMesh::Point(0., 0.) - dshift + doffset,
        libMesh::Point(dDIM1v+dDIM2v, 0.) - dshift + doffset,
        libMesh::Point(dDIM1v+dDIM2v, dt) - dshift + doffset,
        libMesh::Point(0.5*dDIM1v+dDIM2v, dt) - dshift + doffset,
        libMesh::Point(0.5*dDIM1v+dDIM2v, dt+dDIM3v) - dshift + doffset,
        libMesh::Point(dDIM1v+dDIM2v, dt+dDIM3v) - dshift + doffset,
        libMesh::Point(dDIM1v+dDIM2v, dDIM4v) - dshift + doffset,
        libMesh::Point(0., dDIM4v) - dshift + doffset,
        libMesh::Point(0., dt+dDIM3v) - dshift + doffset,
        libMesh::Point(0.5*dDIM1v, dt+dDIM3v) - dshift + doffset,
        libMesh::Point(0.5*dDIM1v, dt) - dshift + doffset,
        libMesh::Point(0., dt) - dshift + doffset
    };
    
    return points;
}


const libMesh::Point 
MAST::Solid1DI1SectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real offset_y, offset_z;
    hy_off(p, t, offset_y);     
    hz_off(p, t, offset_z);
    
    /*!
     * I1 section is bisymmetric and thus the shear center is located at the
     * centroid.  In this case, the origin lies on the centroid.
     */
    libMesh::Point offset(offset_z, offset_y);
    return libMesh::Point(0., 0., 0.) + offset;
}


const libMesh::Point 
MAST::Solid1DI1SectionElementPropertyCard::get_shear_center_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real doffset_y, doffset_z;
    hy_off.derivative(f, p, t, doffset_y);
    hz_off.derivative(f, p, t, doffset_z);
    
    /*!
     * I1 section is bisymmetric and thus the shear center is located at the
     * centroid. In this case, the origin lies on the centroid.
     */
    libMesh::Point doffset(doffset_z, doffset_y);
    return libMesh::Point(0., 0., 0.) + doffset;
}


const libMesh::Point 
MAST::Solid1DI1SectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real offset_y, offset_z;
    hy_off(p, t, offset_y);
    hz_off(p, t, offset_z);
    
    /*!
     * In this case, the origin lies on the centroid.
     */
    libMesh::Point offset(offset_z, offset_y);
    return libMesh::Point(0., 0., 0.) + offset;
}



const libMesh::Point 
MAST::Solid1DI1SectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real doffset_y, doffset_z;
    hy_off.derivative(f, p, t, doffset_y);
    hz_off.derivative(f, p, t, doffset_z);
    
    /*!
     * In this case, the origin lies on the centroid.
     */
    libMesh::Point doffset(doffset_z, doffset_y);
    return libMesh::Point(0., 0., 0.) + doffset;
}

const std::vector<libMesh::Point> 
MAST::Solid1DI1SectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
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
    
    // Points
    libMesh::Point shift(0., 0., 0.);
    std::vector<libMesh::Point> points = {
        libMesh::Point(0.5*(DIM1v+DIM2v), 0.5*DIM4v) - shift - ps + offset,
        libMesh::Point(0.5*(DIM1v+DIM2v), -0.5*DIM4v) - shift - ps + offset,
        libMesh::Point(-0.5*(DIM1v+DIM2v), -0.5*DIM4v) - shift - ps + offset,
        libMesh::Point(-0.5*(DIM1v+DIM2v), 0.5*DIM4v) - shift - ps + offset
    };
    
    return points;
};

const std::vector<libMesh::Point> 
MAST::Solid1DI1SectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
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
    
    // Points
    std::vector<libMesh::Point> points = {
        libMesh::Point(0.5*(dDIM1v+dDIM2v), 0.5*dDIM4v) - dps + doffset,
        libMesh::Point(0.5*(dDIM1v+dDIM2v), -0.5*dDIM4v) - dps + doffset,
        libMesh::Point(-0.5*(dDIM1v+dDIM2v), -0.5*dDIM4v) - dps + doffset,
        libMesh::Point(-0.5*(dDIM1v+dDIM2v), 0.5*dDIM4v) - dps + doffset
    };
}

void MAST::Solid1DI1SectionElementPropertyCard::init(const libMesh::LibMeshInit& init) {
    
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
    else if (DIM4v<=DIM3v){
        libmesh_error_msg("DIM4<=DIM3");
    }
    
    // Create a cross section model of this section
    cross_section.reset(new MAST::CrossSection(init, 3500, *this, libMesh::TRI6));
    
//     _A.reset(new MAST::Solid1D4ParameterSectionProperty::Area(MAST::Solid1DI1SectionProperty::calcA,
//                                                               MAST::Solid1DI1SectionProperty::calcdA,
//                                                               DIM1, DIM2, DIM3, 
//                                                               DIM4));
//     
//     _Ay.reset(new MAST::Solid1D4ParameterSectionProperty::AreaYMoment(
//                                                                 MAST::Solid1DI1SectionProperty::calcA,
//                                                                 MAST::Solid1DI1SectionProperty::calcdA,
//                                                                 DIM1, DIM2, 
//                                                                 DIM3, DIM4,
//                                                                 hz_off));
//     
//     _Az.reset(new MAST::Solid1D4ParameterSectionProperty::AreaZMoment(
//                                                                 MAST::Solid1DI1SectionProperty::calcA,
//                                                                 MAST::Solid1DI1SectionProperty::calcdA,
//                                                                 DIM1, DIM2, 
//                                                                 DIM3, DIM4,
//                                                                 hy_off));
//     
//     _J.reset(new MAST::Solid1D4ParameterSectionProperty::TorsionalConstant(
//                                                                 MAST::Solid1DI1SectionProperty::calcJ,
//                                                                 MAST::Solid1DI1SectionProperty::calcdJ,
//                                                                 DIM1, DIM2,
//                                                                 DIM3, DIM4));
//     
//     _Ip.reset(new MAST::Solid1D4ParameterSectionProperty::PolarInertia(
//                                                                 MAST::Solid1DI1SectionProperty::calcIp,
//                                                                 MAST::Solid1DI1SectionProperty::calcdIp,
//                                                                 MAST::Solid1DI1SectionProperty::calcA,
//                                                                 MAST::Solid1DI1SectionProperty::calcdA,
//                                                                 DIM1, DIM2,
//                                                                 DIM3, DIM4,
//                                                                 hy_off,
//                                                                 hz_off));
//     
//     _AI.reset(new MAST::Solid1D4ParameterSectionProperty::AreaInertiaMatrix(
//                                                                 MAST::Solid1DI1SectionProperty::calcIz,
//                                                                 MAST::Solid1DI1SectionProperty::calcdIz,
//                                                                 MAST::Solid1DI1SectionProperty::calcIy,
//                                                                 MAST::Solid1DI1SectionProperty::calcdIy,
//                                                                 MAST::Solid1DI1SectionProperty::calcA,
//                                                                 MAST::Solid1DI1SectionProperty::calcdA,
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


void MAST::Solid1DI1SectionElementPropertyCard::create_cross_section(
    const libMesh::LibMeshInit& init, const uint n_target_elems, 
    const libMesh::ElemType element_type)
{
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, 
                                               *this, element_type));
}


void MAST::Solid1DI1SectionElementPropertyCard::calculate_properties_pilkey()
{
    cross_section->calculate_geometric_properties();
    cross_section->calculate_warping_properties();
}
