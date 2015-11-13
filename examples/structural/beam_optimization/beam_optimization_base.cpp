/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "examples/structural/beam_optimization/beam_optimization_base.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"


MAST::MultilinearInterpolation::MultilinearInterpolation
(const std::string& nm,
 std::map<Real, MAST::FieldFunction<Real>*>& values):
MAST::FieldFunction<Real>(nm),
_values(values) {
    
    // make sure that the size of the provided values is finite
    libmesh_assert(values.size() > 0);
    
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it = values.begin(), end = values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        _functions.insert(it->second->master());
}

MAST::MultilinearInterpolation::MultilinearInterpolation
(const MAST::MultilinearInterpolation& o):
MAST::FieldFunction<Real>(o),
_values(o._values) {
    std::map<Real, MAST::FieldFunction<Real>*>::iterator
    it = _values.begin(), end = _values.end();
    
    // tell the function that it is dependent on the provided functions
    for ( ; it != end; it++)
        _functions.insert(it->second->master());
}


std::auto_ptr<MAST::FieldFunction<Real> >
MAST::MultilinearInterpolation::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::MultilinearInterpolation(*this));
}

MAST::MultilinearInterpolation::~MultilinearInterpolation() {
    
}


void
MAST::MultilinearInterpolation::operator() (const libMesh::Point& p,
                                            Real t,
                                            Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        (*it1->second)(p, t, f0);
        (*it2->second)(p, t, f1);
        // now interpolate
        v =  (f0 +
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}

void
MAST::MultilinearInterpolation::derivative(const MAST::DerivativeType d,
                                           const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           Real t,
                                           Real& v) const {
    
    //
    // the following is used for calculation of the return value
    //   f(x) is defined for x for each x0 < x < x1
    //   if   x <= x0,      f(x) = f(x0)
    //   if   x0 < x < x1,  f(x) is interpolated
    //   if   x >= x1,      f(x) = f(x1)
    //
    
    std::map<Real, MAST::FieldFunction<Real>*>::const_iterator
    it1, it2;
    std::map<Real, MAST::FieldFunction<Real>*>::const_reverse_iterator
    rit = _values.rbegin();
    it1  = _values.begin();
    
    // check the lower bound
    if (p(0) <=  it1->first) {
        (*it1->second)(p, t, v);
    }
    // check the upper bound
    else if (p(0) >=  rit->first) {
        (*rit->second)(p, t, v);
    }
    else {
        // if it gets here, the ordinate is in between the provided range
        it2 = _values.lower_bound(p(0));
        // this cannot be the first element of the map
        libmesh_assert(it2 != _values.begin());
        // it2 provides the upper bound. The lower bound is provided by the
        // preceding iterator
        it1 = it2--;
        Real f0 = 0., f1 = 0.;
        it1->second->derivative(d, f, p, t, f0);
        it2->second->derivative(d, f, p, t, f1);
        // now interpolate
        v =  (f0 +
              (p(0) - it1->first)/(it2->first - it1->first) *
              (f1-f0));
    }
}




MAST::BeamOffset::BeamOffset(const std::string& nm,
                             MAST::FieldFunction<Real> *thickness):
MAST::FieldFunction<Real>(nm),
_dim(thickness) {
    
    _functions.insert(thickness->master());
}

MAST::BeamOffset::BeamOffset(const MAST::BeamOffset& o):
MAST::FieldFunction<Real>(o),
_dim(o._dim->clone().release()) {
    
    _functions.insert(_dim->master());
}


std::auto_ptr<MAST::FieldFunction<Real> >
MAST::BeamOffset::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::BeamOffset(*this));
}


MAST::BeamOffset::~BeamOffset() {
    delete _dim;
}


void
MAST::BeamOffset::operator() (const libMesh::Point& p,
                              Real t,
                              Real& v) const {
    
    (*_dim)(p, t, v);
    v *= 0.5;
}


void
MAST::BeamOffset::derivative(const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             Real t,
                             Real& v) const {
    _dim->derivative(d, f, p, t, v);
    v *= 0.5;
}




MAST::Weight::Weight(MAST::PhysicsDisciplineBase& discipline):
MAST::FieldFunction<Real>("Weight"),
_discipline(discipline)
{ }



MAST::Weight::Weight(const MAST::Weight& w):
MAST::FieldFunction<Real>(w),
_discipline(w._discipline)
{ }



std::auto_ptr<MAST::FieldFunction<Real> >
MAST::Weight::clone() const {
    
    return std::auto_ptr<MAST::FieldFunction<Real> >
    (new MAST::Weight(*this));
}



MAST::Weight::~Weight() { }



void
MAST::Weight::operator() (const libMesh::Point& p,
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
        // weight from the section area and section density
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
MAST::Weight::derivative(const MAST::DerivativeType d,
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
        // weight from the section area and section density
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


