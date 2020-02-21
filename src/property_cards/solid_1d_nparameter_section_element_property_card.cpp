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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

// MAST includes
#include "property_cards/solid_1d_nparameter_section_element_property_card.h"

MAST::Solid1DnParameterSectionProperty::Area::
Area(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("Area"),
_cross_section(cross_section)
{
}
    
void MAST::Solid1DnParameterSectionProperty::Area::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = _cross_section.get_area(p, t);
}
    
void MAST::Solid1DnParameterSectionProperty::Area::
derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
           Real& m) const 
{
    m = _cross_section.get_area_derivative(f, p, t);
}


MAST::Solid1DnParameterSectionProperty::AreaYMoment::
AreaYMoment(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("AreaYMoment"),
_cross_section(cross_section)
{
}

void MAST::Solid1DnParameterSectionProperty::AreaYMoment::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = _cross_section.get_first_area_moment_y(p, t);
}

void MAST::Solid1DnParameterSectionProperty::AreaYMoment::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            Real& m) const 
{
    m = _cross_section.get_first_area_moment_y_derivative(f, p, t);
}
    


MAST::Solid1DnParameterSectionProperty::AreaZMoment::
AreaZMoment(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("AreaZMoment"),
_cross_section(cross_section)
{
}
    
void MAST::Solid1DnParameterSectionProperty::AreaZMoment::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = _cross_section.get_first_area_moment_z(p, t);
}
    
void MAST::Solid1DnParameterSectionProperty::AreaZMoment::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            Real& m) const
{
    m = _cross_section.get_first_area_moment_z_derivative(f, p, t);
}



MAST::Solid1DnParameterSectionProperty::AreaInertiaMatrix::
AreaInertiaMatrix(MAST::CrossSection& cross_section):
MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
_cross_section(cross_section)
{
}

void MAST::Solid1DnParameterSectionProperty::AreaInertiaMatrix::
operator() (const libMesh::Point& p, const Real t, RealMatrixX& m) const 
{
    m = RealMatrixX::Zero(2,2);
    
    RealVectorX I = _cross_section.get_second_area_moments(p, t);
    m(0,0) = I(0);
    m(1,1) = I(1);
    m(0,1) = m(1,0) = I(2);
}

void MAST::Solid1DnParameterSectionProperty::AreaInertiaMatrix::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            RealMatrixX& m) const 
{
    m = RealMatrixX::Zero(2,2);
    
    RealVectorX dI = _cross_section.get_second_area_moments_derivative(f, p, t);
    
    m(0,0) = dI(0);
    m(1,1) = dI(1);
    m(0,1) = m(1,0) = dI(2);
}


MAST::Solid1DnParameterSectionProperty::PolarInertia::
PolarInertia(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("PolarInertia"),
_cross_section(cross_section)
{
}

void MAST::Solid1DnParameterSectionProperty::PolarInertia::
operator() (const libMesh::Point& p,const Real t, Real& m) const 
{
    m = _cross_section.get_second_polar_area_moment(p, t);
}
    

void MAST::Solid1DnParameterSectionProperty::PolarInertia::
derivative (const MAST::FunctionBase& f, const libMesh::Point& p,
            const Real t, Real& m) const 
{
    m = _cross_section.get_second_polar_area_moment_derivative(f, p, t);
}


MAST::Solid1DnParameterSectionProperty::TorsionalConstant::
TorsionalConstant(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("TorsionalConstant"),
_cross_section(cross_section)
{
}

void MAST::Solid1DnParameterSectionProperty::TorsionalConstant::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = _cross_section.get_torsion_constant(p, t);
}

void MAST::Solid1DnParameterSectionProperty::TorsionalConstant::
derivative (MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m)
{
    m = _cross_section.get_torsion_constant_derivative(f, p, t);
}


MAST::Solid1DnParameterSectionProperty::WarpingConstant::
WarpingConstant(MAST::CrossSection& cross_section):
MAST::FieldFunction<Real>("WarpingConstant"),
_cross_section(cross_section)
{
}
    

void MAST::Solid1DnParameterSectionProperty::WarpingConstant::
operator() (const libMesh::Point& p, const Real t, Real& m) const 
{
    m = _cross_section.get_warping_constant(p, t);
}

    
void MAST::Solid1DnParameterSectionProperty::WarpingConstant::
derivative (MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
Real& m)
{
    m = _cross_section.get_warping_constant_derivative(f, p, t);
}


MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix::
ShearCoefficientMatrix(MAST::CrossSection& cross_section):
MAST::FieldFunction<RealMatrixX>("ShearCoefficientMatrix"),
_cross_section(cross_section)
{
}

void MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix::
operator() (const libMesh::Point& p, const Real t, RealMatrixX& m) const 
{
    m = RealMatrixX::Zero(2,2);
    
    RealVectorX kappa = _cross_section.get_shear_coefficients(p, t);
    m(0,0) = kappa(0);
    m(1,1) = kappa(1);
    m(0,1) = m(1,0) = kappa(2);
}

void MAST::Solid1DnParameterSectionProperty::ShearCoefficientMatrix::
derivative (MAST::FunctionBase& f, const libMesh::Point& p, const Real t,
            RealMatrixX& m)
{
    m = RealMatrixX::Zero(2,2);
    
    RealVectorX dkappa = _cross_section.get_shear_coefficients_derivative(f, p, t);
    
    m(0,0) = dkappa(0);
    m(1,1) = dkappa(1);
    m(0,1) = m(1,0) = dkappa(2);
}
