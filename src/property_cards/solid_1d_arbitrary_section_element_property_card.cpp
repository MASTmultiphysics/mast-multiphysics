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

//TODO: Arbitrary cross sections currently only support sensitivites w.r.t. 
// offsets. Other sensitivies will be returned as zero.

// MAST includes
#include "property_cards/solid_1d_arbitrary_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"
#include "base/elem_base.h"

#define PI 3.1415926535897932

namespace MAST {
    namespace Solid1DArbitrarySectionProperty {
        
        class Area: public MAST::FieldFunction<Real> {
        public:
            Area(const Real& A):
            MAST::FieldFunction<Real>("Area"), _A(A) {}
            
            virtual ~Area() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                m = _A;
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                m = 0.0;
            }
        protected:
            const Real& _A;
            
        };
        
        
        
        class TorsionalConstant: public MAST::FieldFunction<Real> {
        public:
            TorsionalConstant(const Real& T):
            MAST::FieldFunction<Real>("TorsionalConstant"), _T(T) {}
            
            virtual ~TorsionalConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {             
                m = _T;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                m = 0.0;
            }
        protected:
            const Real& _T;
            
        };
        
        
        
        class PolarInertia: public MAST::FieldFunction<Real> {
        public:
            PolarInertia(const Real& J, const Real& A,
                         const MAST::FieldFunction<Real>&  hy_offset,
                         const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("PolarInertia"),
            _J(J),  _A(A),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~PolarInertia() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real offy, offz;
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m = _J + _A*pow(offy,2) + _A*pow(offz,2);
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real dA, offy, offz, doffy, doffz;
                _hy_offset (p, t, offy);  _hy_offset.derivative( f, p, t, doffy);
                _hz_offset (p, t, offz);  _hz_offset.derivative( f, p, t, doffz);
                dA = 0.0;
                
                m = 0.0 + pow(offy,2)*dA + 2.*_A*offy*doffy + pow(offz,2)*dA + 2.*_A*offz*doffz;
            }
            
        protected:
            const Real& _J, _A;
            const MAST::FieldFunction<Real> &_hy_offset, &_hz_offset;
        };
        
        
        /*!
         *   calculates the area moment about the Y-axis due to an offset 
         *   along the Z-axis
         */
        class AreaYMoment: public MAST::FieldFunction<Real> {
        public:
            AreaYMoment(const Real& A,
                        const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _A(A),
            _hz_offset(hz_offset) {
                _functions.insert(&hz_offset);
            }
            
            
            virtual ~AreaYMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real off;
                _hz_offset(p, t, off);
                
                m = _A * off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real off, doff;
                _hz_offset (p, t, off); _hz_offset.derivative( f, p, t, doff);
                
                m = 0.0 + _A*doff;
            }
            
        protected:
            const Real& _A;
            const MAST::FieldFunction<Real> &_hz_offset;
        };
        
        
        
        /*!
         *   calculates the area moment about the Z-axis due to an offset
         *   along the Y-axis
         */
        class AreaZMoment: public MAST::FieldFunction<Real> {
        public:
            AreaZMoment(const Real& A,
                        const MAST::FieldFunction<Real>&  hy_offset):
            MAST::FieldFunction<Real>("AreaZMoment"),
            _A(A),
            _hy_offset(hy_offset) {
                _functions.insert(&hy_offset);
            }
            
            virtual ~AreaZMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real off;
                _hy_offset(p, t, off);
                
                m = _A*off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real off, doff;
                _hy_offset(p, t, off); _hy_offset.derivative( f, p, t, doff);
                
                m = 0.0 + _A*doff;
            }
            
        protected:
            const Real& _A;
            const MAST::FieldFunction<Real> &_hy_offset;
        };
        
        
        
        /*!
         *   calculates the 2x2 matrix of area inertia for the section with 
         *   individual entries as 
         *
         *   0 x 0 = int_omega  (y+yoff)^2 dy dz
         *   0 x 1 = int_omega  (y+yoff) (z+zoff) dy dz
         *   1 x 0 = int_omega  (y+yoff) (z+zoff) dy dz
         *   1 x 1 = int_omega  (z+zoff)^2 dy dz
         */
        class AreaInertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            AreaInertiaMatrix(const Real& Izz, const Real& Iyy, const Real& A,
                              const MAST::FieldFunction<Real>&  hy_offset,
                              const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
            _Izz(Izz), _Iyy(Iyy), _A(A),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~AreaInertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real offy, offz;
                m = RealMatrixX::Zero(2,2);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m(0,0) = _Izz + _A*pow(offy,2); // Izz for v-bending
                m(0,1) = _A * offy * offz;
                m(1,0) = m(0,1);
                m(1,1) = _Iyy + _A*pow(offz,2); // Iyy for w-bending
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real dA, offy, offz, doffy, doffz;
                m = RealMatrixX::Zero(2,2);
                _hy_offset(p, t, offy); _hy_offset.derivative( f, p, t, doffy);
                _hz_offset(p, t, offz); _hz_offset.derivative( f, p, t, doffz);
                dA = 0.0;
                
                m(0,0) = 0.0 + dA*pow(offy,2) + _A*2.0*offy;
                m(0,1) = dA*offy*offz + _A*doffy*offz + _A*offy*doffz;
                m(1,0) = m(0,1);
                m(1,1) = 0.0 + dA*pow(offz,2) + _A*2.0*offz;
            }
            
        protected:
            const Real& _Izz, _Iyy, _A;
            const MAST::FieldFunction<Real> &_hy_offset, &_hz_offset;
        };
    }
}
 

void MAST::Solid1DArbitrarySectionElementPropertyCard::init(libMesh::MeshBase& mesh){
    
    libmesh_assert(!_initialized);
    
    this->calculateGeometricProperties(mesh);
    
    MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(_A_val));
    _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
    _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
    _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
    _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
                                                                      hy_off,
                                                                      hz_off));
    _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
                                                                           _Iyy_val,
                                                                    _A_val, hy_off,
                                                                  hz_off));
    
    _initialized = true;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::init(RealMatrixX& vertices){
    
    libmesh_assert(!_initialized);
    
    this->calculateGeometricProperties(vertices);
    
    MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(_A_val));
    _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
    _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
    _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
    _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
                                                                      hy_off,
                                                                      hz_off));
    _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
                                                                           _Iyy_val,
                                                                    _A_val, hy_off,
                                                                  hz_off));
    
    _initialized = true;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::init(Real A, Real Izz, 
                                                            Real Iyy, Real Ip, 
                                                            Real J){
    libmesh_assert(!_initialized);
    
    _A_val = A;
    _Izz_val = Izz;
    _Iyy_val = Iyy;
    _Ip_val = Ip;
    _J_val = J;
    
    MAST::FieldFunction<Real>
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(_A_val));
    _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
    _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
    _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
    _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
                                                                      hy_off,
                                                                      hz_off));
    _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
                                                                           _Iyy_val,
                                                                    _A_val, hy_off,
                                                                  hz_off));
    
    _initialized = true;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::clear() {
    
    libmesh_assert(!_initialized);
    
    _A.reset();
    _Ay.reset();
    _Az.reset();
    _J.reset();
    _Ip.reset();
    _AI.reset();
    
    _A_val=0;
    _Izz_val=0;
    _Iyy_val=0;
    _Ip_val=0;
    _J_val=0;
    
    _torsionConstantSet = false;
    
    _initialized = false;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::setTorsionalConstant(Real J){
    _J_val = J;
    _torsionConstantSet = true;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::calculateGeometricProperties(libMesh::MeshBase& mesh){
    Real C_x, C_y, A, Ix, Iy, Ixy;
    C_x = 0.0; C_y = 0.0; A = 0.0; Ix = 0.0; Iy = 0.0; Ixy = 0.0;
    libMesh::MeshBase::const_element_iterator el_it = mesh.elements_begin();
    libMesh::MeshBase::const_element_iterator end_el_it = mesh.elements_end();
    for ( ; el_it != end_el_it; el_it++)
    {
        for (uint i=0; i<(*el_it)->n_nodes(); i++)
        {            
            Real xi = (*el_it)->point(i)(0);
            Real yi = (*el_it)->point(i)(1);
            Real xi1, yi1;
            
            if (i==(*el_it)->n_nodes()-1)
            {
                xi1 = (*el_it)->point(0)(0);
                yi1 = (*el_it)->point(0)(1);
            }
            else
            {
                xi1 = (*el_it)->point(i+1)(0);
                yi1 = (*el_it)->point(i+1)(1);
            }
            Real Ibase = (xi*yi1 - xi1*yi);
            Ix += Ibase*(yi*yi + yi*yi1 + yi1*yi1);
            Iy += Ibase*(xi*xi + xi*xi1 + xi1*xi1);
            Ixy += Ibase*(xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi);
        }
        A += (*el_it)->volume(); // Returns the area for 2D elements.
        C_x += (*el_it)->centroid()(0)*(*el_it)->volume();
        C_y += (*el_it)->centroid()(1)*(*el_it)->volume();
    }
    C_x /= A;
    C_y /= A;
    Ix /= 12.0;
    Iy /= 12.0;
    Ixy /= 24.0;
    _A_val = A;
    _Izz_val = Ix-A*C_y*C_y;
    _Iyy_val= Iy - A*C_x*C_x;
    _Ip_val = _Izz_val + _Iyy_val;  // Perpendicular Axis Theory
    
    if (not _torsionConstantSet)
    {
        // A very very rough approximation to the torsion constant
        Real a = sqrt(A); // Length of side of square with equivalent area as mesh
        _J_val = 0.189*pow(a,4);
        // FIXME: For arbitrary sections, it may be necessary to perform FEA 
    }
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::calculateGeometricProperties(RealMatrixX& vertices){
    Real C_x, C_y, A, Ix, Iy, Ixy;
    C_x = 0.0; C_y = 0.0; A = 0.0; Ix = 0.0; Iy = 0.0; Ixy = 0.0;
    for (uint i=0; i<vertices.rows(); i++)
    {        
        Real xi = vertices(i,0);
        Real yi = vertices(i,1);
        Real xi1, yi1, yin1;
        
        if (i==0)
        {
            yin1 = vertices(vertices.rows()-1,1);
        }
        else
        {
            yin1 = vertices(i-1,1);
        }
        
        if (i==vertices.rows()-1)
        {
            xi1 = vertices(0,0);
            yi1 = vertices(0,1);
        }
        else
        {
            xi1 = vertices(i+1,0);
            yi1 = vertices(i+1,1);
        }
        C_x += xi;
        C_y += yi;
        A += xi*(yi1 - yin1); // Shoelace Formula
        Real Ibase = (xi*yi1 - xi1*yi);
        Ix += Ibase*(yi*yi + yi*yi1 + yi1*yi1);
        Iy += Ibase*(xi*xi + xi*xi1 + xi1*xi1);
        Ixy += Ibase*(xi*yi1 + 2*xi*yi + 2*xi1*yi1 + xi1*yi);
    }
    C_x /= double(vertices.rows());
    C_y /= double(vertices.rows());
    A = abs(A)/2.0; // Shoelace Formula
    Ix /= 12.0;
    Iy /= 12.0;
    Ixy /= 24.0;
    _A_val = A;
    _Izz_val = Ix-A*C_y*C_y;
    _Iyy_val= Iy - A*C_x*C_x;
    _Ip_val = _Izz_val + _Iyy_val;   // Perpendicular Axis Theory
    
    if (not _torsionConstantSet)
    {
        // A very very rough approximation to the torsion constant
        Real a = sqrt(A); // Length of side of square with equivalent area as mesh
        _J_val = 0.189*pow(a,4);
        // FIXME: For arbitrary sections, it may be necessary to perform FEA 
    }
}
