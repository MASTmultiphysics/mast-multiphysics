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
#include "property_cards/solid_1d_nparameter_section_element_property_card.h"
#include "property_cards/solid_1d_bar_section_element_property_card.h"


// Rectangle (BAR in Siemens NX Nastran and Astros 21.2)
void MAST::Solid1DBarSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = DIM1*DIM2;
}

void MAST::Solid1DBarSectionProperty::calcdA(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dA){
    dA = DIM1*dDIM2+DIM2*dDIM1;
}

void MAST::Solid1DBarSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (DIM1*(DIM2*DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIz(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIz){
    dIz = ((DIM2*DIM2*DIM2)*dDIM1)/1.2E+1+(DIM1*(DIM2*DIM2)*dDIM2)/4.0;
}

void MAST::Solid1DBarSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = ((DIM1*DIM1*DIM1)*DIM2)/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIy(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIy){
    dIy = ((DIM1*DIM1*DIM1)*dDIM2)/1.2E+1+((DIM1*DIM1)*DIM2*dDIM1)/4.0;
}

void MAST::Solid1DBarSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (DIM1*DIM2*(DIM1*DIM1+DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcdIp(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dIp){
    dIp = (DIM1*dDIM2*(DIM1*DIM1+(DIM2*DIM2)*3.0))/1.2E+1+(DIM2*dDIM1*((DIM1*DIM1)*3.0+DIM2*DIM2))/1.2E+1;
}

void MAST::Solid1DBarSectionProperty::calcJ1(Real& DIM1, Real& DIM2, Real& J1){
    J1 = (1.0/(DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2*DIM2)*((DIM1*DIM1*DIM1*DIM1)*DIM2*-2.52E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*4.0E+2+(DIM2*DIM2*DIM2*DIM2*DIM2)*2.1E+1))/1.2E+3;
}

void MAST::Solid1DBarSectionProperty::calcdJ1(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ1){
    dJ1 = (1.0/(DIM1*DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2*DIM2)*dDIM1*((DIM1*DIM1*DIM1*DIM1*DIM1)*1.0E+2-(DIM2*DIM2*DIM2*DIM2*DIM2)*2.1E+1))/3.0E+2+(1.0/(DIM1*DIM1*DIM1*DIM1)*(DIM2*DIM2)*dDIM2*((DIM1*DIM1*DIM1*DIM1)*DIM2*-4.2E+1+(DIM1*DIM1*DIM1*DIM1*DIM1)*5.0E+1+(DIM2*DIM2*DIM2*DIM2*DIM2)*7.0))/5.0E+1;
}

void MAST::Solid1DBarSectionProperty::calcJ2(Real& DIM1, Real& DIM2, Real& J2){
    J2 = ((DIM1*DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2)*(DIM1*(DIM2*DIM2*DIM2*DIM2)*-2.52E+2+(DIM1*DIM1*DIM1*DIM1*DIM1)*2.1E+1+(DIM2*DIM2*DIM2*DIM2*DIM2)*4.0E+2))/1.2E+3;
}

void MAST::Solid1DBarSectionProperty::calcdJ2(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ2){
    dJ2 = (DIM1*DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2*DIM2)*dDIM2*((DIM1*DIM1*DIM1*DIM1*DIM1)*2.1E+1-(DIM2*DIM2*DIM2*DIM2*DIM2)*1.0E+2)*(-1.0/3.0E+2)+((DIM1*DIM1)*1.0/(DIM2*DIM2*DIM2*DIM2)*dDIM1*(DIM1*(DIM2*DIM2*DIM2*DIM2)*-4.2E+1+(DIM1*DIM1*DIM1*DIM1*DIM1)*7.0+(DIM2*DIM2*DIM2*DIM2*DIM2)*5.0E+1))/5.0E+1;
}

void MAST::Solid1DBarSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    if (DIM1>DIM2)
    {
        MAST::Solid1DBarSectionProperty::calcJ1(DIM1, DIM2, J);
    }
    else
    {
        MAST::Solid1DBarSectionProperty::calcJ2(DIM1, DIM2, J);
    }
}


void MAST::Solid1DBarSectionProperty::calcdJ(Real& DIM1, Real& DIM2, Real& dDIM1, Real& dDIM2, Real& dJ){
    if (DIM1>DIM2)
    {
        MAST::Solid1DBarSectionProperty::calcdJ1(DIM1, DIM2, dDIM1, dDIM2, dJ);
    }
    else
    {
        MAST::Solid1DBarSectionProperty::calcdJ2(DIM1, DIM2, dDIM1, dDIM2, dJ);
    }
}


const std::vector<libMesh::Point> 
MAST::Solid1DBarSectionElementPropertyCard::get_geom_points(const libMesh::Point& p, const Real t, const uint n) const 
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v); 
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
        
    std::vector<libMesh::Point> points = {
        libMesh::Point(-0.5*DIM1v, -0.5*DIM2v) + offset,
        libMesh::Point( 0.5*DIM1v, -0.5*DIM2v) + offset,
        libMesh::Point( 0.5*DIM1v,  0.5*DIM2v) + offset,
        libMesh::Point(-0.5*DIM1v,  0.5*DIM2v) + offset
    };
    
    return points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DBarSectionElementPropertyCard::get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const 
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1.derivative(f, p, t, dDIM1v);
    DIM2.derivative(f, p, t, dDIM2v);
    hy_off.derivative(f, p, t, doffset_y);
    hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    
    std::vector<libMesh::Point> points = {
        libMesh::Point(-0.5*DIM1v, -0.5*DIM2v) + offset,
        libMesh::Point(0.5*DIM1v, -0.5*DIM2v) + offset,
        libMesh::Point(0.5*DIM1v, 0.5*DIM2v) + offset,
        libMesh::Point(-0.5*DIM1v, 0.5*DIM2v) + offset
    };
    
    std::vector<libMesh::Point> dpoints = {
        libMesh::Point(-0.5*dDIM1v, -0.5*dDIM2v) + doffset,
        libMesh::Point(0.5*dDIM1v, -0.5*dDIM2v) + doffset,
        libMesh::Point(0.5*dDIM1v, 0.5*dDIM2v) + doffset,
        libMesh::Point(-0.5*dDIM1v, 0.5*dDIM2v) + doffset
    };
    
    return dpoints;
}


const libMesh::Point
MAST::Solid1DBarSectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_centroid(p, t);
}


const libMesh::Point
MAST::Solid1DBarSectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_shear_center(p, t);
}


const libMesh::Point
MAST::Solid1DBarSectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const 
{
    return cross_section->get_centroid_derivative(f, p, t);
}


// const libMesh::Point 
// MAST::Solid1DBarSectionElementPropertyCard::get_shear_center_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const 
// {
//     return cross_section->get_shear_center_derivative(f, p, t);
// }


const libMesh::Point 
MAST::Solid1DBarSectionElementPropertyCard::get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t) 
{
    return cross_section->get_shear_center_derivative(f, p, t);
}


const std::vector<libMesh::Point> 
MAST::Solid1DBarSectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    // Ordered C, D, E, F as defined in MSC Nastran manual (See PBARL or PBEAML)
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    // Points Relative to Point ps
    std::vector<libMesh::Point> points = {
        libMesh::Point(0.5*DIM1v, 0.5*DIM2v) + offset - ps,
        libMesh::Point(0.5*DIM1v, -0.5*DIM2v) + offset - ps,
        libMesh::Point(-0.5*DIM1v, -0.5*DIM2v) + offset - ps,
        libMesh::Point(-0.5*DIM1v, 0.5*DIM2v) + offset - ps
    };
    
    return points;
};


const std::vector<libMesh::Point> 
MAST::Solid1DBarSectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
{
    // Ordered C, D, E, F as defined in MSC Nastran manual (See PBARL or PBEAML)
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v); DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1.derivative(f, p, t, dDIM1v);
    DIM2.derivative(f, p, t, dDIM2v);
    hy_off.derivative(f, p, t, doffset_y);
    hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    
    // Derivative of Stress Points Relative to ps
    // FIXME: Need to account for changes in ps, i.e. when ps is centroid or shear center
    std::vector<libMesh::Point> points = {
        libMesh::Point( 0.5*dDIM1v,  0.5*dDIM2v) + doffset - dps,
        libMesh::Point( 0.5*dDIM1v, -0.5*dDIM2v) + doffset - dps,
        libMesh::Point(-0.5*dDIM1v, -0.5*dDIM2v) + doffset - dps,
        libMesh::Point(-0.5*dDIM1v,  0.5*dDIM2v) + doffset - dps
    };
    
    return points;
};


void MAST::Solid1DBarSectionElementPropertyCard::init(const libMesh::LibMeshInit& init) {
    
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
    
    // Create a cross section model of this section
    cross_section.reset(new MAST::CrossSection(init, 3500, *this, libMesh::TRI6));
    
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

void MAST::Solid1DBarSectionElementPropertyCard::create_cross_section(
    const libMesh::LibMeshInit& init, const uint n_target_elems, 
    const libMesh::ElemType element_type)
{
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, 
                                               *this, element_type));
}


void MAST::Solid1DBarSectionElementPropertyCard::calculate_properties_pilkey()
{
    cross_section->calculate_geometric_properties();
    cross_section->calculate_warping_properties();
}
