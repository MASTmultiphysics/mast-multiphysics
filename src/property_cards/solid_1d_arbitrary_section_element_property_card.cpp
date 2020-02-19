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
            Area(const MAST::FieldFunction<Real>& A):
            MAST::FieldFunction<Real>("Area"), _A(A){
                _functions.insert(&A);
            }
            
            virtual ~Area() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A;
                _A(p, t, A);
                m = A;
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, dA;
                _A(p, t, A); _A.derivative(f, p, t, dA);
                m = dA;
            }
        protected:
            const MAST::FieldFunction<Real>& _A;
            
        };
        
        
        class TorsionalConstant: public MAST::FieldFunction<Real> {
        public:
            TorsionalConstant(const MAST::FieldFunction<Real>& T):
            MAST::FieldFunction<Real>("TorsionalConstant"), _T(T){
                _functions.insert(&T);
            }
            
            virtual ~TorsionalConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {          
                Real T;
                _T(p, t, T);
                m = T;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real T, dT;
                _T(p, t, T); _T.derivative(f, p, t, dT);
                m = dT;
            }
        protected:
            const MAST::FieldFunction<Real>& _T;
            
        };
        
        
        class WarpingConstant: public MAST::FieldFunction<Real> {
        public:
            WarpingConstant(const MAST::FieldFunction<Real>& W):
            MAST::FieldFunction<Real>("WarpingConstant"), _W(W)
            {
                _functions.insert(&W);
            }
            
            virtual ~WarpingConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real W;
                _W(p, t, W);
                m = W;
            }
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real W, dW;
                _W(p, t, W); _W.derivative(f, p, t, dW);
                m = dW;
            }
        protected:
            const MAST::FieldFunction<Real>& _W;
            
        };
        
        
        class ShearCoefficientMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ShearCoefficientMatrix(const MAST::FieldFunction<Real>& Kzz, 
                                   const MAST::FieldFunction<Real>& Kyy):
            MAST::FieldFunction<RealMatrixX>("ShearCoefficientMatrix"), 
            _Kzz(Kzz),
            _Kyy(Kyy)
            {
                _functions.insert(&Kzz);
                _functions.insert(&Kyy);
            }
            
            virtual ~ShearCoefficientMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real Kzz;
                Real Kyy;
                _Kzz(p, t, Kzz);
                _Kyy(p, t, Kyy);
                m = RealMatrixX::Zero(2,2);
                m(0,0) = Kzz;
                m(1,1) = Kyy;
            }
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real dKzz, dKyy;
                _Kyy.derivative(f, p, t, dKyy);
                _Kzz.derivative(f, p, t, dKzz);
                m = RealMatrixX::Zero(2,2);
                m(0,0) = dKzz;
                m(1,1) = dKyy;
            }
        protected:
            const MAST::FieldFunction<Real>& _Kyy;
            const MAST::FieldFunction<Real>& _Kzz;
        };
        
        
        class PolarInertia: public MAST::FieldFunction<Real> {
        public:
            PolarInertia(const MAST::FieldFunction<Real>& J, 
                         const MAST::FieldFunction<Real>& A,
                         const MAST::FieldFunction<Real>&  hy_offset,
                         const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("PolarInertia"),
            _J(J),  _A(A),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&J);
                _functions.insert(&A);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~PolarInertia() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, J, offy, offz;
                _A(p, t, A);
                _J(p, t, J);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m = J + A*pow(offy,2) + A*pow(offz,2);
            }
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, dA, J, dJ, offy, offz, doffy, doffz;
                _A(p, t, A);              _A.derivative(f, p, t, dA);
                _J(p, t, J);              _J.derivative(f, p, t, dJ);
                _hy_offset (p, t, offy);  _hy_offset.derivative( f, p, t, doffy);
                _hz_offset (p, t, offz);  _hz_offset.derivative( f, p, t, doffz);
                
                m = dJ + pow(offy,2)*dA + 2.*A*offy*doffy + pow(offz,2)*dA + 2.*A*offz*doffz;
            }
            
        protected:
            const MAST::FieldFunction<Real> &_J, &_A, &_hy_offset, &_hz_offset;
        };
        
        
        /*!
         *   calculates the area moment about the Y-axis due to an offset 
         *   along the Z-axis
         */
        class AreaYMoment: public MAST::FieldFunction<Real> {
        public:
            AreaYMoment(const MAST::FieldFunction<Real>& A,
                        const MAST::FieldFunction<Real>& hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _A(A),
            _hz_offset(hz_offset) {
                _functions.insert(&A);
                _functions.insert(&hz_offset);
            }
            
            
            virtual ~AreaYMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, off;
                _A(p, t, A);
                _hz_offset(p, t, off);
                
                m = A * off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, dA, off, doff;
                _A(p, t, A);            _A.derivative(f, p, t, dA);
                _hz_offset (p, t, off); _hz_offset.derivative( f, p, t, doff);
                
                m = dA*off + A*doff;
            }
            
        protected:
            const MAST::FieldFunction<Real>& _A;
            const MAST::FieldFunction<Real> &_hz_offset;
        };
        
        
        
        /*!
         *   calculates the area moment about the Z-axis due to an offset
         *   along the Y-axis
         */
        class AreaZMoment: public MAST::FieldFunction<Real> {
        public:
            AreaZMoment(const MAST::FieldFunction<Real>& A,
                        const MAST::FieldFunction<Real>& hy_offset):
            MAST::FieldFunction<Real>("AreaZMoment"),
            _A(A),
            _hy_offset(hy_offset) {
                _functions.insert(&A);
                _functions.insert(&hy_offset);
            }
            
            virtual ~AreaZMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, off;
                _A(p, t, A);
                _hy_offset(p, t, off);
                
                m = A*off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real A, dA, off, doff;
                _A(p, t, A);           _A.derivative(f, p, t, dA);
                _hy_offset(p, t, off); _hy_offset.derivative( f, p, t, doff);
                
                m = dA*off + A*doff;
            }
            
        protected:
            const MAST::FieldFunction<Real>& _A;
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
            AreaInertiaMatrix(const MAST::FieldFunction<Real>& Izz, 
                              const MAST::FieldFunction<Real>& Iyy, 
                              const MAST::FieldFunction<Real>& Izy,
                              const MAST::FieldFunction<Real>& A,
                              const MAST::FieldFunction<Real>&  hy_offset,
                              const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
            _Izz(Izz), _Iyy(Iyy), _Izy(Izy), _A(A),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&Izz);
                _functions.insert(&Iyy);
                _functions.insert(&Izy);
                _functions.insert(&A);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~AreaInertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real Izz, Iyy, Izy, A, offy, offz;
                m = RealMatrixX::Zero(2,2);
                _A(p, t, A);
                _Izz(p, t, Izz);
                _Iyy(p, t, Iyy);
                _Izy(p, t, Izy);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m(0,0) = Izz + A*pow(offy,2); // Izz for v-bending
                m(0,1) = Izy + A * offy * offz; // Izy
                m(1,0) = m(0,1);
                m(1,1) = Iyy + A*pow(offz,2); // Iyy for w-bending
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real A, dA, Izz, dIzz, Iyy, dIyy, Izy, dIzy;
                Real offy, offz, doffy, doffz;
                m = RealMatrixX::Zero(2,2);
                _A(p, t, A);            _A.derivative(f, p, t, dA);
                _Izz(p, t, Izz);        _Izz.derivative(f, p, t, dIzz);
                _Iyy(p, t, Iyy);        _Iyy.derivative(f, p, t, dIyy);
                _Izy(p, t, Izy);        _Izy.derivative(f, p, t, dIzy);
                _hy_offset(p, t, offy); _hy_offset.derivative(f, p, t, doffy);
                _hz_offset(p, t, offz); _hz_offset.derivative(f, p, t, doffz);
                
                m(0,0) = dIzz + dA*pow(offy,2) + A*2.0*offy*doffy;
                m(0,1) = dIzy + dA*offy*offz + A*doffy*offz + A*offy*doffz;
                m(1,0) = m(0,1);
                m(1,1) = dIyy + dA*pow(offz,2) + A*2.0*offz*doffz;
            }
            
        protected:
            const MAST::FieldFunction<Real> &_Izz, &_Iyy, &_Izy,  &_A;
            const MAST::FieldFunction<Real> &_hy_offset, &_hz_offset;
        };
    }
}
 

//  // FIXME: Need this to autocreate FieldFunction<Real> from calculated properties
// void MAST::Solid1DArbitrarySectionElementPropertyCard::init(libMesh::MeshBase& mesh){
//     
//     libmesh_assert(!_initialized);
//     
//     this->calculateGeometricProperties(mesh);
//     
//     MAST::FieldFunction<Real>
//     &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
//     &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
//     
// //     _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(_A_val));
//     _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
//     _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
//     _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
//     _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
//                                                                       hy_off,
//                                                                       hz_off));
//     _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
//                                                                            _Iyy_val,
//                                                                     _A_val, hy_off,
//                                                                   hz_off));
//     
//     _initialized = true;
// }


// // FIXME: Need this to autocreate FieldFunction<Real> from calculated properties
// void MAST::Solid1DArbitrarySectionElementPropertyCard::init(RealMatrixX& vertices){
//     
//     libmesh_assert(!_initialized);
//     
//     this->calculateGeometricProperties(vertices);
//     
//     MAST::FieldFunction<Real>
//     &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
//     &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
//     
// //     _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(_A_val));
//     _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
//     _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
//     _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
//     _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
//                                                                       hy_off,
//                                                                       hz_off));
//     _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
//                                                                            _Iyy_val,
//                                                                     _A_val, hy_off,
//                                                                   hz_off));
//     
//     _initialized = true;
// }


// // FIXME: Need this to autocreate FieldFunction<Real> from specified properties
// void MAST::Solid1DArbitrarySectionElementPropertyCard::init(Real A, Real Izz, 
//                                                             Real Iyy, Real Ip, 
//                                                             Real J){
//     libmesh_assert(!_initialized);
//     
//     _A_val = A;
//     _Izz_val = Izz;
//     _Iyy_val = Iyy;
//     _Ip_val = Ip;
//     _J_val = J;
//     
//     MAST::FieldFunction<Real>
//     &Af      =  this->get<MAST::FieldFunction<Real> >("A"),
//     &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
//     &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
//     
//     _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(Af));
//     _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(_A_val, hz_off));
//     _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(_A_val, hy_off));
//     _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(_J_val));
//     _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(_Ip_val, _A_val,
//                                                                       hy_off,
//                                                                       hz_off));
//     _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(_Izz_val,
//                                                                            _Iyy_val,
//                                                                     _A_val, hy_off,
//                                                                   hz_off));
//     
//     _initialized = true;
// }


void MAST::Solid1DArbitrarySectionElementPropertyCard::init()
{
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &A      =  this->get<MAST::FieldFunction<Real> >("A"),
    &Izz    =  this->get<MAST::FieldFunction<Real> >("Izz"),
    &Iyy    =  this->get<MAST::FieldFunction<Real> >("Iyy"),
    &Izy    =  this->get<MAST::FieldFunction<Real> >("Izy"),
    &Ip     =  this->get<MAST::FieldFunction<Real> >("Ip"),
    &J      =  this->get<MAST::FieldFunction<Real> >("J"),
    &W      =  this->get<MAST::FieldFunction<Real> >("W"),
    &Kappazz    =  this->get<MAST::FieldFunction<Real> >("Kappazz"),
    &Kappayy    =  this->get<MAST::FieldFunction<Real> >("Kappayy"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DArbitrarySectionProperty::Area(A));
    _Ay.reset(new MAST::Solid1DArbitrarySectionProperty::AreaYMoment(A, hz_off));
    _Az.reset(new MAST::Solid1DArbitrarySectionProperty::AreaZMoment(A, hy_off));
    _J.reset(new MAST::Solid1DArbitrarySectionProperty::TorsionalConstant(J));
    _Ip.reset(new MAST::Solid1DArbitrarySectionProperty::PolarInertia(Ip, A,
                                                                      hy_off,
                                                                      hz_off));
    _AI.reset(new MAST::Solid1DArbitrarySectionProperty::AreaInertiaMatrix(Izz,
                                                  Iyy, Izy, A, hy_off, hz_off));
    _Gamma.reset(new MAST::Solid1DArbitrarySectionProperty::WarpingConstant(W));
    _Kappa.reset(new MAST::Solid1DArbitrarySectionProperty::ShearCoefficientMatrix(Kappazz, Kappayy));
    
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
    _Kappa.reset();
    _Gamma.reset();
    
    _A_val=0;
    _Izz_val=0;
    _Iyy_val=0;
    _Ip_val=0;
    _J_val=0;
    _W_val=0;
    _Kxx_val=0;
    _Kyy_val=0;
    
    _torsionConstantSet = false;
    
    _initialized = false;
}


void MAST::Solid1DArbitrarySectionElementPropertyCard::setTorsionalConstant(Real J){
    _J_val = J;
    _torsionConstantSet = true;
}

const std::vector<libMesh::Point> MAST::Solid1DArbitrarySectionElementPropertyCard::get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    return _stress_points;
}

void MAST::Solid1DArbitrarySectionElementPropertyCard::add_stress_point(const libMesh::Point stress_point)
{
    _stress_points.push_back(stress_point);
}

const libMesh::Point MAST::Solid1DArbitrarySectionElementPropertyCard::get_centroid(const libMesh::Point& p, const Real t) const
{
    return libMesh::Point(0., 0., 0.);
}

const libMesh::Point MAST::Solid1DArbitrarySectionElementPropertyCard::get_shear_center(const libMesh::Point& p, const Real t) const
{
    return libMesh::Point(0., 0., 0.);
}

const libMesh::Point MAST::Solid1DArbitrarySectionElementPropertyCard::get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t) const
{
    return libMesh::Point(0., 0., 0.);
}

const libMesh::Point MAST::Solid1DArbitrarySectionElementPropertyCard::get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, const Real t)
{
    return libMesh::Point(0., 0., 0.);
}

const std::vector<libMesh::Point> MAST::Solid1DArbitrarySectionElementPropertyCard::get_stress_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, const Real t, const libMesh::Point dps) const
{
    std::vector<libMesh::Point> dstress_points(_stress_points.size());
    for (uint i=0; i<_stress_points.size(); i++)
    {
        dstress_points[i] = libMesh::Point(0., 0., 0.);
    }
    return dstress_points;
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
