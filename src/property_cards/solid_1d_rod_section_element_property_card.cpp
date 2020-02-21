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

// C++ includes
#include <math.h>

// MAST includes
#include "property_cards/solid_1d_rod_section_element_property_card.h"


#define PI 3.1415926535897932


void MAST::Solid1DRodSectionProperty::calcA(Real& DIM1, Real& A){
    A = PI*DIM1*DIM1;
}

void MAST::Solid1DRodSectionProperty::calcdA(Real& DIM1, Real& dDIM1, Real& dA){
    dA = 2.*PI*DIM1*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIz(Real& DIM1, Real& Iz){
    Iz = PI*pow(DIM1,4.)/4.0;
}

void MAST::Solid1DRodSectionProperty::calcdIz(Real& DIM1, Real& dDIM1, Real& dIz){
    dIz = PI*pow(DIM1,3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIy(Real& DIM1, Real& Iy){
    Iy = PI*pow(DIM1,4.)/4.0; 
}

void MAST::Solid1DRodSectionProperty::calcdIy(Real& DIM1, Real& dDIM1, Real& dIy){
    dIy = PI*pow(DIM1,3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcIp(Real& DIM1, Real& Ip){
    Ip = PI*pow(DIM1, 4.)/2.0;
}

void MAST::Solid1DRodSectionProperty::calcdIp(Real& DIM1, Real& dDIM1, Real& dIp){
    dIp = 2.*PI*pow(DIM1, 3.)*dDIM1;
}

void MAST::Solid1DRodSectionProperty::calcJ(Real& DIM1, Real& J){
    J = PI*pow(DIM1, 4.)/2.0;
}

void MAST::Solid1DRodSectionProperty::calcdJ(Real& DIM1, Real& dDIM1, Real& dJ){
    dJ = 2.*PI*pow(DIM1, 3.)*dDIM1;
}


MAST::Solid1DRodSectionProperty::WarpingConstant::
WarpingConstant():
MAST::FieldFunction<Real>("WarpingConstant")
{
}

void MAST::Solid1DRodSectionProperty::WarpingConstant::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = 0.0;
}


void MAST::Solid1DRodSectionProperty::WarpingConstant::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m) const
{
    m = 0.0;
}


void MAST::Solid1DRodSectionProperty::WarpingConstant::
derivative (MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m)
{
    m = 0.0;
}

const std::vector<libMesh::Point>
MAST::Solid1DRodSectionElementPropertyCard::get_geom_points(const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    DIM1(p, t, DIM1v);
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
MAST::Solid1DRodSectionElementPropertyCard::get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    Real dDIM1v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
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


const libMesh::Point 
MAST::Solid1DRodSectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    DIM1(p, t, DIM1v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    return libMesh::Point(offset_z, offset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DRodSectionElementPropertyCard::get_shear_center_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    Real dDIM1v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    return libMesh::Point(doffset_z, doffset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DRodSectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    DIM1(p, t, DIM1v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    return libMesh::Point(offset_z, offset_y, 0.);
}


const libMesh::Point 
MAST::Solid1DRodSectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    Real dDIM1v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
    hy_off(p, t, offset_y); hy_off.derivative(f, p, t, doffset_y);
    hz_off(p, t, offset_z); hz_off.derivative(f, p, t, doffset_z);
    
    return libMesh::Point(doffset_z, doffset_y, 0.);
}


const std::vector<libMesh::Point> 
MAST::Solid1DRodSectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    DIM1(p, t, DIM1v);
    hy_off(p, t, offset_y); hz_off(p, t, offset_z);
    
    libMesh::Point offset(offset_z, offset_y);
    
    return {libMesh::Point(0., DIM1v, 0.) + offset - ps,
            libMesh::Point(DIM1v, 0., 0.) + offset - ps,
            libMesh::Point(0. -DIM1v, 0.) + offset - ps,
            libMesh::Point(-DIM1v, 0., 0.) + offset - ps
    };
}


const std::vector<libMesh::Point> 
MAST::Solid1DRodSectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
{
    const MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real DIM1v, offset_y, offset_z;
    Real dDIM1v, doffset_y, doffset_z;
    DIM1(p, t, DIM1v);      DIM1.derivative(f, p, t, dDIM1v);
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


void MAST::Solid1DRodSectionElementPropertyCard::init(const libMesh::LibMeshInit& init) {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &DIM1     =  this->get<MAST::FieldFunction<Real> >("DIM1"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    // Check that dimensions are physically correct
    Real DIM1v; DIM1(DIM1v);
    if (DIM1v<=0){
        libmesh_error_msg("DIM1<=0");
    }
    
    // Create a cross section model of this section
    cross_section.reset(new MAST::CrossSection(init, 3500, *this, libMesh::TRI6));
    
    _A.reset(new MAST::Solid1D1ParameterSectionProperty::Area(MAST::Solid1DRodSectionProperty::calcA,
                                                              MAST::Solid1DRodSectionProperty::calcdA,
                                                              DIM1));
    
    _Ay.reset(new MAST::Solid1D1ParameterSectionProperty::AreaYMoment(
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1, 
                                                                hz_off));
    
    _Az.reset(new MAST::Solid1D1ParameterSectionProperty::AreaZMoment(
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1, 
                                                                hy_off));
    
    _J.reset(new MAST::Solid1D1ParameterSectionProperty::TorsionalConstant(
                                                                MAST::Solid1DRodSectionProperty::calcJ,
                                                                MAST::Solid1DRodSectionProperty::calcdJ,
                                                                DIM1));
    
    _Ip.reset(new MAST::Solid1D1ParameterSectionProperty::PolarInertia(
                                                                MAST::Solid1DRodSectionProperty::calcIp,
                                                                MAST::Solid1DRodSectionProperty::calcdIp,
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _AI.reset(new MAST::Solid1D1ParameterSectionProperty::AreaInertiaMatrix(
                                                                MAST::Solid1DRodSectionProperty::calcIz,
                                                                MAST::Solid1DRodSectionProperty::calcdIz,
                                                                MAST::Solid1DRodSectionProperty::calcIy,
                                                                MAST::Solid1DRodSectionProperty::calcdIy,
                                                                MAST::Solid1DRodSectionProperty::calcA,
                                                                MAST::Solid1DRodSectionProperty::calcdA,
                                                                DIM1,
                                                                hy_off,
                                                                hz_off));
    
    _Gamma.reset(new MAST::Solid1DRodSectionProperty::WarpingConstant());
    
    _Kappa.reset(new MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix(*cross_section));
    
    _initialized = true;
}


void MAST::Solid1DRodSectionElementPropertyCard::create_cross_section(
    const libMesh::LibMeshInit& init, const uint n_target_elems, 
    const libMesh::ElemType element_type)
{
    cross_section.reset(new MAST::CrossSection(init, n_target_elems, 
                                               *this, element_type));
}


void MAST::Solid1DRodSectionElementPropertyCard::calculate_properties_pilkey()
{
    cross_section->calculate_geometric_properties();
    cross_section->calculate_warping_properties();
}
