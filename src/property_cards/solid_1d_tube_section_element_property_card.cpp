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
#include "property_cards/solid_1d_tube_section_element_property_card.h"


#define PI 3.1415926535897932


// Annulus (TUBE in Siemens NX Nastran and Astros 21.2)
void MAST::Solid1DTubeSectionProperty::calcA(Real& DIM1, Real& DIM2, Real& A){
    A = PI*(DIM1*DIM1-DIM2*DIM2);
}

void MAST::Solid1DTubeSectionProperty::calcdA(Real& DIM1, Real& DIM2, 
                                              Real& dDIM1, Real& dDIM2,
                                              Real& dA){
    dA = PI*(DIM1*dDIM1-DIM2*dDIM2)*2.0;
}

void MAST::Solid1DTubeSectionProperty::calcIz(Real& DIM1, Real& DIM2, Real& Iz){
    Iz = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/4.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIz(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2, 
                                               Real& dIz){
    dIz = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2);
}

void MAST::Solid1DTubeSectionProperty::calcIy(Real& DIM1, Real& DIM2, Real& Iy){
    Iy = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/4.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIy(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2, 
                                               Real& dIy){
    dIy = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2);
}

void MAST::Solid1DTubeSectionProperty::calcIp(Real& DIM1, Real& DIM2, Real& Ip){
    Ip = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/2.0;
}

void MAST::Solid1DTubeSectionProperty::calcdIp(Real& DIM1, Real& DIM2, 
                                               Real& dDIM1, Real& dDIM2,
                                               Real& dIp){
    dIp = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2)*2.0;
}

void MAST::Solid1DTubeSectionProperty::calcJ(Real& DIM1, Real& DIM2, Real& J){
    J = (PI*(DIM1*DIM1*DIM1*DIM1-DIM2*DIM2*DIM2*DIM2))/2.0;
}

void MAST::Solid1DTubeSectionProperty::calcdJ(Real& DIM1, Real& DIM2, 
                                              Real& dDIM1, Real& dDIM2, 
                                              Real& dJ){
    dJ = PI*((DIM1*DIM1*DIM1)*dDIM1-(DIM2*DIM2*DIM2)*dDIM2)*2.0;
}


const std::vector<libMesh::Point>
MAST::Solid1DTubeSectionElementPropertyCard::get_geom_points(const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v);      DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    std::vector<libMesh::Point> geom_points(n);
    for (uint i=0; i<n; i++)
    {
        Real theta = (2.*PI)/Real(n+1) * i;
        geom_points[i] = libMesh::Point(DIM1v*cos(theta), DIM1v*sin(theta), 0.) + offset;
    }
    
    return geom_points;
}


const std::vector<libMesh::Point>
MAST::Solid1DTubeSectionElementPropertyCard::get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    DIM2(p, t, DIM2v);      DIM2.derivative(f, p, t, dDIM2v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    std::vector<libMesh::Point> dgeom_points(n);
    for (uint i=0; i<n; i++)
    {
        Real theta = (2.*PI)/Real(n+1) * i;
        dgeom_points[i] = libMesh::Point(dDIM1v*cos(theta), dDIM1v*sin(theta), 0.) + doffset;
    }
    
    return dgeom_points;
}


const std::vector<std::vector<libMesh::Point>> 
MAST::Solid1DTubeSectionElementPropertyCard::get_holes_points(const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v);      DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    std::vector<libMesh::Point> hole_points(n);
    for (uint i=0; i<n; i++)
    {
        Real theta = (2.*PI)/Real(n+1) * i;
        hole_points[i] = libMesh::Point(DIM2v*cos(theta), DIM2v*sin(theta), 0.) + offset;
    }
    
    return {hole_points};
}


const std::vector<std::vector<libMesh::Point>>
MAST::Solid1DTubeSectionElementPropertyCard::get_holes_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    DIM2(p, t, DIM2v);      DIM2.derivative(f, p, t, dDIM2v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    std::vector<libMesh::Point> dholes_points(n);
    for (uint i=0; i<n; i++)
    {
        Real theta = (2.*PI)/Real(n+1) * i;
        dholes_points[i] = libMesh::Point(dDIM2v*cos(theta), dDIM2v*sin(theta), 0.) + doffset;
    }
    
    return {dholes_points};
}


const libMesh::Point 
MAST::Solid1DTubeSectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v);      DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    return libMesh::Point(offset_z, offset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DTubeSectionElementPropertyCard::get_shear_center_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    DIM2(p, t, DIM2v);      DIM2.derivative(f, p, t, dDIM2v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    return libMesh::Point(doffset_z, doffset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DTubeSectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v);      DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    return libMesh::Point(offset_z, offset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DTubeSectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    DIM2(p, t, DIM2v);      DIM2.derivative(f, p, t, dDIM2v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    return libMesh::Point(doffset_z, doffset_y, 0.);
}


const std::vector<libMesh::Point> 
MAST::Solid1DTubeSectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    DIM1(p, t, DIM1v);      DIM2(p, t, DIM2v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    return {libMesh::Point(0., DIM1v, 0.) + offset - ps,
            libMesh::Point(DIM1v, 0., 0.) + offset - ps,
            libMesh::Point(0. -DIM1v, 0.) + offset - ps,
            libMesh::Point(-DIM1v, 0., 0.) + offset - ps
    };
}


const std::vector<libMesh::Point> 
MAST::Solid1DTubeSectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &DIM2     =  this->get<MAST::FieldFunction<Real> >("DIM2"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, DIM2v, offset_y, offset_z;
    Real dDIM1v, dDIM2v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    DIM2(p, t, DIM2v);      DIM2.derivative(f, p, t, dDIM2v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    libMesh::Point doffset(doffset_z, doffset_y);
    
    return {libMesh::Point(0., dDIM1v, 0.) + doffset - dps,
            libMesh::Point(dDIM1v, 0., 0.) + doffset - dps,
            libMesh::Point(0. -dDIM1v, 0.) + doffset - dps,
            libMesh::Point(-dDIM1v, 0., 0.) + doffset - dps
    };
}


void MAST::Solid1DTubeSectionElementPropertyCard::init(const libMesh::LibMeshInit& init) 
{
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
    else if (DIM2v>=DIM1v){
        libmesh_error_msg("DIM2>=DIM1");
    }
    
    // Create a cross section model of this section
    cross_section.reset(new MAST::CrossSection(init, 3500, *this, libMesh::TRI6));
    
//     _A.reset(new MAST::Solid1DnParameterSectionProperty::Area(*cross_section));
//     
//     _Ay.reset(new MAST::Solid1DnParameterSectionProperty::AreaYMoment(*cross_section));
//     
//     _Az.reset(new MAST::Solid1DnParameterSectionProperty::AreaZMoment(*cross_section));
//     
//     _AI.reset(new MAST::Solid1DnParameterSectionProperty::AreaInertiaMatrix(*cross_section));
//     
//     _Ip.reset(new MAST::Solid1DnParameterSectionProperty::PolarInertia(*cross_section));
//     
//     _J.reset(new MAST::Solid1DnParameterSectionProperty::TorsionalConstant(*cross_section));
    
    _A.reset(new MAST::Solid1D2ParameterSectionProperty::Area(MAST::Solid1DTubeSectionProperty::calcA,
                                                              MAST::Solid1DTubeSectionProperty::calcdA,
                                                              DIM1, DIM2));
    
    _Ay.reset(new MAST::Solid1D2ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D2ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off));
    
    _AI.reset(new MAST::Solid1D2ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DTubeSectionProperty::calcIz,
                                                                MAST::Solid1DTubeSectionProperty::calcdIz,
                                                                MAST::Solid1DTubeSectionProperty::calcIy,
                                                                MAST::Solid1DTubeSectionProperty::calcdIy,
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _Ip.reset(new MAST::Solid1D2ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DTubeSectionProperty::calcIp,
                                                                MAST::Solid1DTubeSectionProperty::calcdIp,
                                                                MAST::Solid1DTubeSectionProperty::calcA,
                                                                MAST::Solid1DTubeSectionProperty::calcdA,
                                                                DIM1, DIM2,
                                                                hy_off,
                                                                hz_off));
    
    _J.reset(new MAST::Solid1D2ParameterSectionProperty::TorsionalConstant(
                            MAST::Solid1DTubeSectionProperty::calcJ,
                            MAST::Solid1DTubeSectionProperty::calcdJ,
                            DIM1, DIM2));
    
    _Kappa.reset(new MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix(*cross_section));
    
    _Gamma.reset(new MAST::Solid1DnParameterSectionProperty::WarpingConstant(*cross_section));
    
    _initialized = true;
}

void MAST::Solid1DTubeSectionElementPropertyCard::create_cross_section(
    const libMesh::LibMeshInit& init, const uint n_target_elems, 
    const libMesh::ElemType element_type)
{
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, 
                                               *this, element_type));
}


void MAST::Solid1DTubeSectionElementPropertyCard::calculate_properties_pilkey()
{
    cross_section->calculate_geometric_properties();
    cross_section->calculate_warping_properties();
}
