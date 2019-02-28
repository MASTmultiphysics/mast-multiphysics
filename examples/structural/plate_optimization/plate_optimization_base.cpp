///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2019  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//// MAST includes
//#include "examples/base/multilinear_interpolation.h"
//#include "examples/structural/plate_optimization/plate_optimization_base.h"
//#include "property_cards/solid_2d_section_element_property_card.h"
//#include "property_cards/material_property_card_base.h"
//
//
//
//
//MAST::PlateWeight::PlateWeight(MAST::PhysicsDisciplineBase& discipline):
//MAST::FieldFunction<Real>("PlateWeight"),
//_discipline(discipline)
//{ }
//
//
//
//MAST::PlateWeight::~PlateWeight() { }
//
//
//
//void
//MAST::PlateWeight::operator() (const libMesh::Point& p,
//                               Real t,
//                               Real& v) const {
//    
//    // get a reference to the mesh
//    const libMesh::MeshBase&
//    mesh    = _discipline.get_equation_systems().get_mesh();
//    libMesh::MeshBase::const_element_iterator
//    eit     = mesh.active_local_elements_begin(),
//    eend    = mesh.active_local_elements_end();
//    
//    Real h, rho;
//    v = 0.;
//    libMesh::Point elem_p;
//    
//    for ( ; eit != eend; eit++ ) {
//        
//        // get a pointer to the element and then as the discipline
//        // for the element property
//        const libMesh::Elem* e = *eit;
//        
//        const MAST::ElementPropertyCardBase& prop =
//        _discipline.get_property_card(*e);
//        
//        // assuming that the element is one-dimensional, we need
//        // its section area value
//        
//        // before that, convert the property to a 1D section property
//        // card
//        const MAST::Solid2DSectionElementPropertyCard& prop2d =
//        dynamic_cast<const MAST::Solid2DSectionElementPropertyCard&>(prop);
//        
//        // get a reference to the section area
//        const MAST::FieldFunction<Real>& th =
//        prop2d.get<MAST::FieldFunction<Real> >("h");
//        
//        // get a reference to the density variable
//        const MAST::MaterialPropertyCardBase& mat = prop.get_material();
//        const MAST::FieldFunction<Real> &rhof =
//        mat.get<MAST::FieldFunction<Real> >("rho");
//        
//        elem_p  = e->centroid();
//        th  ( elem_p, 0.,   h);
//        rhof( elem_p, 0., rho);
//        v += e->volume() * h * rho;
//    }
//}
//
//
//
//
//
//void
//MAST::PlateWeight::derivative( const MAST::FunctionBase& f,
//                              const libMesh::Point& p,
//                              Real t,
//                              Real& v) const {
//    
//    // get a reference to the mesh
//    const libMesh::MeshBase&
//    mesh    = _discipline.get_equation_systems().get_mesh();
//    libMesh::MeshBase::const_element_iterator
//    eit     = mesh.active_local_elements_begin(),
//    eend    = mesh.active_local_elements_end();
//    
//    Real h, rho, dh, drho;
//    v = 0.;
//    libMesh::Point elem_p;
//    
//    for ( ; eit != eend; eit++ ) {
//        
//        // get a pointer to the element and then as the discipline
//        // for the element property
//        const libMesh::Elem* e = *eit;
//        
//        const MAST::ElementPropertyCardBase& prop =
//        _discipline.get_property_card(*e);
//        
//        // assuming that the element is one-dimensional, we need
//        // its section area value
//        
//        // before that, convert the property to a 1D section property
//        // card
//        const MAST::Solid2DSectionElementPropertyCard& prop2d =
//        dynamic_cast<const MAST::Solid2DSectionElementPropertyCard&>(prop);
//        
//        // get a reference to the section area
//        const MAST::FieldFunction<Real>& th =
//        prop2d.get<MAST::FieldFunction<Real> >("h");
//        
//        // get a reference to the density variable
//        const MAST::MaterialPropertyCardBase& mat =
//        prop.get_material();
//        const MAST::FieldFunction<Real> &rhof =
//        mat.get<MAST::FieldFunction<Real> >("rho");
//        
//        elem_p = e->centroid();
//        // for each element iterate over the length and calculate the
//        // PlateWeight from the section area and section density
//        // use three point trapezoidal rule to calculate the integral
//        th  ( elem_p, 0., h);
//        rhof( elem_p, 0., rho);
//        th.derivative  ( f, elem_p, 0.,   dh);
//        rhof.derivative( f, elem_p, 0., drho);
//        v += e->volume() * (dh * rho + h * drho);
//        
//    }
//}
//
//
