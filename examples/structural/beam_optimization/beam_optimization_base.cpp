/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "examples/base/multilinear_interpolation.h"
#include "examples/structural/beam_optimization/beam_optimization_base.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"



MAST::BeamWeight::BeamWeight(MAST::PhysicsDisciplineBase& discipline):
MAST::FieldFunction<Real>("BeamWeight"),
_discipline(discipline)
{ }



MAST::BeamWeight::~BeamWeight() { }



void
MAST::BeamWeight::operator() (const libMesh::Point& p,
                              Real t,
                              Real& v) const {
    
    // get a reference to the mesh
    const libMesh::MeshBase&
    mesh    = _discipline.get_equation_systems().get_mesh();
    libMesh::MeshBase::const_element_iterator
    eit     = mesh.active_local_elements_begin(),
    eend    = mesh.active_local_elements_end();
    
    Real h, rho, x0, x1, dx;
    v = 0.;
    
    libMesh::Point elem_p;
    const unsigned int n_sec = 3; // number of quadrature divs
    
    for ( ; eit != eend; eit++ ) {
        
        // get a pointer to the element and then as the discipline
        // for the element property
        const libMesh::Elem* e = *eit;
        
        const MAST::ElementPropertyCardBase& prop =
        _discipline.get_property_card(*e);
        
        // assuming that the element is one-dimensional, we need
        // its section area value
        
        // before that, convert the property to a 1D section property
        // card
        const MAST::Solid1DSectionElementPropertyCard& prop1d =
        dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
        
        // get a reference to the section area
        const MAST::FieldFunction<Real>& area = prop1d.A();
        
        // get a reference to the density variable
        const MAST::MaterialPropertyCardBase& mat = prop.get_material();
        const MAST::FieldFunction<Real> &rhof =
        mat.get<MAST::FieldFunction<Real> >("rho");
        
        
        // for each element iterate over the length and calculate the
        // BeamWeight from the section area and section density
        // use three point trapezoidal rule to calculate the integral
        x0 = e->point(0)(0);
        x1 = e->point(1)(0);
        dx = (x1-x0)/n_sec;
        for (unsigned int i=0; i<n_sec; i++) {
            elem_p(0) = x0 + dx*(i+0.5);
            area(elem_p, 0., h);
            rhof(elem_p, 0., rho);
            v += h * rho * dx;
        }
        
    }
}





void
MAST::BeamWeight::derivative(const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             Real t,
                             Real& v) const {
    
    // get a reference to the mesh
    const libMesh::MeshBase&
    mesh    = _discipline.get_equation_systems().get_mesh();
    libMesh::MeshBase::const_element_iterator
    eit     = mesh.active_local_elements_begin(),
    eend    = mesh.active_local_elements_end();
    
    Real h, rho, dh, drho, x0, x1, dx;
    v = 0.;
    
    libMesh::Point elem_p;
    const unsigned int n_sec = 3; // number of quadrature divs
    
    for ( ; eit != eend; eit++ ) {
        
        // get a pointer to the element and then as the discipline
        // for the element property
        const libMesh::Elem* e = *eit;
        
        const MAST::ElementPropertyCardBase& prop =
        _discipline.get_property_card(*e);
        
        // assuming that the element is one-dimensional, we need
        // its section area value
        
        // before that, convert the property to a 1D section property
        // card
        const MAST::Solid1DSectionElementPropertyCard& prop1d =
        dynamic_cast<const MAST::Solid1DSectionElementPropertyCard&>(prop);
        
        // get a reference to the section area
        const MAST::FieldFunction<Real>&
        area = prop1d.A();
        
        // get a reference to the density variable
        const MAST::MaterialPropertyCardBase& mat =
        prop.get_material();
        const MAST::FieldFunction<Real> &rhof =
        mat.get<MAST::FieldFunction<Real> >("rho");
        
        
        // for each element iterate over the length and calculate the
        // BeamWeight from the section area and section density
        // use three point trapezoidal rule to calculate the integral
        x0 = e->point(0)(0);
        x1 = e->point(1)(0);
        dx = (x1-x0)/n_sec;
        for (unsigned int i=0; i<n_sec; i++) {
            elem_p(0) = x0 + dx*(i+0.5);
            area(elem_p, 0., h);
            area.derivative(d, f, elem_p, 0., dh);
            rhof(elem_p, 0., rho);
            rhof.derivative(d, f, elem_p, 0., drho);
            v += (dh * rho + h * drho) * dx;
        }
        
    }
}


