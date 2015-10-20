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
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"
#include "base/elem_base.h"
#include "mesh/local_elem_base.h"



namespace MAST {
    namespace Solid1DSectionProperty {
        
        class Area: public MAST::FieldFunction<Real> {
        public:
            Area(MAST::FieldFunction<Real> *hy,
                 MAST::FieldFunction<Real>* hz):
            MAST::FieldFunction<Real>("Area"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
            }
            
            Area(const MAST::Solid1DSectionProperty::Area &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionProperty::Area(*this));
            }
            
            virtual ~Area() {
                delete _hy;
                delete _hz;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                
                m = hy*hz;
            }
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz;
                (*_hy)(p, t, hy); _hy->derivative(d, f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->derivative(d, f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz;
        };
        
        
        
        class TorsionalConstant: public MAST::FieldFunction<Real> {
        public:
            TorsionalConstant(MAST::FieldFunction<Real> *hy,
                              MAST::FieldFunction<Real>* hz):
            MAST::FieldFunction<Real>("TorsionalConstant"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
            }
            
            TorsionalConstant(const MAST::Solid1DSectionProperty::TorsionalConstant &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionProperty::TorsionalConstant(*this));
            }
            
            virtual ~TorsionalConstant() {
                delete _hy;
                delete _hz;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, a, b;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                
                
                // shorter side is b, and longer side is a
                if (hy > hz) {
                    a = hy;
                    b = hz;
                }
                else {
                    a = hz;
                    b = hy;
                }
                
                m = a*pow(b,3)*(1./3.-.21*b/a*(1.-pow(b/a,4)/12.));
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz, a, b, da, db;
                (*_hy)(p, t, hy); _hy->derivative(d, f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->derivative(d, f, p, t, dhz);
                
                // shorter side is b, and longer side is a
                if (hy > hz) {
                    a = hy; da = dhy;
                    b = hz; db = dhz;
                }
                else {
                    a = hz; da = dhz;
                    b = hy; db = dhy;
                }
                
                m =
                da*pow(b,3)*(1./3.-.21*b/a*(1.-pow(b/a,4)/12.)) +
                a*3.*pow(b,2)*db*(1./3.-.21*b/a*(1.-pow(b/a,4)/12.)) +
                a*pow(b,3)*(-.21*db/a*(1.-pow(b/a,4)/12.) +
                            (.21*b/pow(a,2)*da*(1.-pow(b/a,4)/12.)) +
                            (-.21*b/a*(-4.*pow(b,3)*db/pow(a,4)/12.+
                                       4.*pow(b,4)/pow(a,5)*da/12.)));
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz;
        };
        
        
        
        class PolarInertia: public MAST::FieldFunction<Real> {
        public:
            PolarInertia(MAST::FieldFunction<Real> *hy,
                         MAST::FieldFunction<Real>* hz,
                         MAST::FieldFunction<Real>* hy_offset,
                         MAST::FieldFunction<Real>* hz_offset):
            MAST::FieldFunction<Real>("PolarInertia"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
                _functions.insert(hy_offset->master());
                _functions.insert(hz_offset->master());
            }
            
            PolarInertia(const MAST::Solid1DSectionProperty::PolarInertia &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hy_offset(f._hy_offset->clone().release()),
            _hz_offset(f._hz_offset->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
                _functions.insert(_hy_offset->master());
                _functions.insert(_hz_offset->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionProperty::PolarInertia(*this));
            }
            
            virtual ~PolarInertia() {
                delete _hy;
                delete _hz;
                delete _hy_offset;
                delete _hz_offset;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, offy, offz;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hy_offset)(p, t, offy);
                (*_hz_offset)(p, t, offz);
                
                m = hy*hz*((pow(hy,2) + pow(hz,2))/12. +
                           pow(offy,2) + pow(offz,2));
            }
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz, offy, offz, doffy, doffz;
                (*_hy)        (p, t, hy);           _hy->derivative(d, f, p, t, dhy);
                (*_hz)        (p, t, hz);           _hz->derivative(d, f, p, t, dhz);
                (*_hy_offset) (p, t, offy);  _hy_offset->derivative(d, f, p, t, doffy);
                (*_hz_offset) (p, t, offz);  _hz_offset->derivative(d, f, p, t, doffz);
                
                
                m =
                (dhy*hz + hy*dhz) * ((pow(hy,2) + pow(hz,2))/12. +
                                     pow(offy,2) + pow(offz,2)) +
                2.*hy*hz*((hy*dhy + hz*dhz)/12. +
                          offy*doffy + offz*doffz);
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hy_offset, *_hz_offset;
        };
        
        
        
        
        class AreaYMoment: public MAST::FieldFunction<Real> {
        public:
            AreaYMoment(MAST::FieldFunction<Real>* hy,
                        MAST::FieldFunction<Real>* hz,
                        MAST::FieldFunction<Real>* hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _hy(hy),
            _hz(hz),
            _hz_offset(hz_offset) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
                _functions.insert(hz_offset->master());
            }
            
            AreaYMoment(const MAST::Solid1DSectionProperty::AreaYMoment &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hz_offset(f._hz_offset->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
                _functions.insert(_hz_offset->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionProperty::AreaYMoment(*this));
            }
            
            virtual ~AreaYMoment() {
                delete _hy;
                delete _hz;
                delete _hz_offset;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hz_offset)(p, t, off);
                
                m = hy*hz*off;
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)        (p, t, hy);         _hy->derivative(d, f, p, t, dhy);
                (*_hz)        (p, t, hz);         _hz->derivative(d, f, p, t, dhz);
                (*_hz_offset) (p, t, off); _hz_offset->derivative(d, f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hz_offset;
        };
        
        
        
        class AreaZMoment: public MAST::FieldFunction<Real> {
        public:
            AreaZMoment(MAST::FieldFunction<Real>* hy,
                        MAST::FieldFunction<Real>* hz,
                        MAST::FieldFunction<Real>* hy_offset):
            MAST::FieldFunction<Real>("AreaZMoment"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
                _functions.insert(hy_offset->master());
            }
            
            AreaZMoment(const MAST::Solid1DSectionProperty::AreaZMoment &f):
            MAST::FieldFunction<Real>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hy_offset(f._hy_offset->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
                _functions.insert(_hy_offset->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Solid1DSectionProperty::AreaZMoment(*this));
            }
            
            virtual ~AreaZMoment() {
                delete _hy;
                delete _hz;
                delete _hy_offset;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off;
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hy_offset)(p, t, off);
                
                m = hy*hz*off;
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                (*_hy)(p, t, hy); _hy->derivative(d, f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->derivative(d, f, p, t, dhz);
                (*_hy_offset)(p, t, off); _hy_offset->derivative(d, f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hy_offset;
        };
        
        
        
        
        class AreaInertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            AreaInertiaMatrix(MAST::FieldFunction<Real>* hy,
                              MAST::FieldFunction<Real>* hz,
                              MAST::FieldFunction<Real>* hy_offset,
                              MAST::FieldFunction<Real>* hz_offset):
            MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(hy->master());
                _functions.insert(hz->master());
                _functions.insert(hy_offset->master());
                _functions.insert(hz_offset->master());
            }
            
            AreaInertiaMatrix(const MAST::Solid1DSectionProperty::AreaInertiaMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _hy(f._hy->clone().release()),
            _hz(f._hz->clone().release()),
            _hy_offset(f._hy_offset->clone().release()),
            _hz_offset(f._hz_offset->clone().release()) {
                _functions.insert(_hy->master());
                _functions.insert(_hz->master());
                _functions.insert(_hy_offset->master());
                _functions.insert(_hz_offset->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::AreaInertiaMatrix(*this));
            }
            
            virtual ~AreaInertiaMatrix() {
                delete _hy;
                delete _hz;
                delete _hy_offset;
                delete _hz_offset;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real hy, hz, offy, offz;
                m = RealMatrixX::Zero(2,2);
                (*_hy)(p, t, hy);
                (*_hz)(p, t, hz);
                (*_hy_offset)(p, t, offy);
                (*_hz_offset)(p, t, offz);
                
                m(0,0) = hz*pow(hy,3)/12. + hy*hz*pow(offy,2) ; // Izz for v-bending
                m(0,1) = hy*hz*offy*offz;
                m(1,0) = m(0,1);
                m(1,1) = hy*pow(hz,3)/12. + hy*hz*pow(offz,2) ; // Iyy for w-bending
            }
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real hy, hz, offy, offz, dhy, dhz, doffy, doffz;
                m = RealMatrixX::Zero(2,2);
                (*_hy)(p, t, hy); _hy->derivative(d, f, p, t, dhy);
                (*_hz)(p, t, hz); _hz->derivative(d, f, p, t, dhz);
                (*_hy_offset)(p, t, offy); _hy_offset->derivative(d, f, p, t, doffy);
                (*_hz_offset)(p, t, offz); _hz_offset->derivative(d, f, p, t, doffz);
                
                
                m(0,0) = dhz*pow(hy,3)/12. + hz*pow(hy,2)/4.*dhy +
                dhy*hz*pow(offy,2) + hy*dhz*pow(offy,2) + 2.*hy*hz*offy*doffy ;
                m(0,1) = dhy*hz*offy*offz + hy*dhz*offy*offz +
                hy*hz*doffy*offz + hy*hz*offy*doffz;
                m(1,0) = m(0,1);
                m(1,1) = dhy*pow(hz,3)/12. + hy*pow(hz,2)/4.*dhz +
                dhy*hz*pow(offz,2) + hy*dhz*pow(offz,2) + 2.*hy*hz*offz*doffz ;
            }
            
        protected:
            
            MAST::FieldFunction<Real> *_hy, *_hz, *_hy_offset, *_hz_offset;
        };
        
        
        class ExtensionStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                     MAST::FieldFunction<Real>* A,
                                     MAST::FieldFunction<Real>* J);
            
            ExtensionStiffnessMatrix(const MAST::Solid1DSectionProperty::ExtensionStiffnessMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A(f._A->clone().release()),
            _J(f._J->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_A->master());
                _functions.insert(_J->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::ExtensionStiffnessMatrix(*this));
            }
            
            virtual ~ExtensionStiffnessMatrix() {
                delete _material_stiffness;
                delete _A;
                delete _J;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<Real> *_A, *_J;
        };
        
        
        
        class ExtensionBendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionBendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                            MAST::FieldFunction<Real>* A_y_moment,
                                            MAST::FieldFunction<Real>* A_z_moment);
            
            ExtensionBendingStiffnessMatrix(const MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_A_y_moment->master());
                _functions.insert(_A_z_moment->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix(*this));
            }
            
            virtual ~ExtensionBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<Real> *_A_y_moment, *_A_z_moment;
        };
        
        
        class BendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            BendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                   MAST::FieldFunction<RealMatrixX> *I);
            
            BendingStiffnessMatrix(const MAST::Solid1DSectionProperty::BendingStiffnessMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _I(f._I->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_I->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::BendingStiffnessMatrix(*this));
            }
            
            virtual ~BendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _I;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<RealMatrixX> *_I;
        };
        
        
        class TransverseStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                      MAST::FieldFunction<Real>* A):
            MAST::FieldFunction<RealMatrixX>("TransverseStiffnessMatrix1D"),
            _material_stiffness(mat),
            _A(A) {
                _functions.insert(mat->master());
                _functions.insert(A->master());
            }
            
            
            TransverseStiffnessMatrix(const MAST::Solid1DSectionProperty::TransverseStiffnessMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _A(f._A->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_A->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::TransverseStiffnessMatrix(*this));
            }
            
            virtual ~TransverseStiffnessMatrix() {
                delete _material_stiffness;
                delete _A;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real A;
                (*_A)(p, t, A);
                (*_material_stiffness)(p, t, m);
                m *= A;
            }
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                RealMatrixX dm;
                Real A, dA;
                (*_A)                  (p, t, A);                  _A->derivative(d, f, p, t, dA);
                (*_material_stiffness) (p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
                
                m *= dA;
                m += A*dm;
            }
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<Real> *_A;
        };
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(MAST::FieldFunction<Real>* rho,
                          MAST::FieldFunction<Real>* A,
                          MAST::FieldFunction<Real>* A_y_moment,
                          MAST::FieldFunction<Real>* A_z_moment,
                          MAST::FieldFunction<Real>* Ip,
                          MAST::FieldFunction<RealMatrixX>* I);
            
            InertiaMatrix(const MAST::Solid1DSectionProperty::InertiaMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _rho(f._rho->clone().release()),
            _A(f._A->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release()),
            _Ip(f._Ip->clone().release()),
            _I(f._I->clone().release()) {
                _functions.insert(_rho->master());
                _functions.insert(_A->master());
                _functions.insert(_A_y_moment->master());
                _functions.insert(_A_z_moment->master());
                _functions.insert(_Ip->master());
                _functions.insert(_I->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::InertiaMatrix(*this));
            }
            
            virtual ~InertiaMatrix() {
                delete _rho;
                delete _A;
                delete _A_y_moment;
                delete _A_z_moment;
                delete _Ip;
                delete _I;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real> *_rho, *_A, *_A_y_moment, *_A_z_moment, *_Ip;
            MAST::FieldFunction<RealMatrixX> *_I;
        };
        
        
        
        class ThermalExpansionAMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionAMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                                    MAST::FieldFunction<RealMatrixX> *mat_expansion,
                                    MAST::FieldFunction<Real> *A);
            
            ThermalExpansionAMatrix(const MAST::Solid1DSectionProperty::ThermalExpansionAMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _A(f._A->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_material_expansion->master());
                _functions.insert(_A->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::ThermalExpansionAMatrix(*this));
            }
            
            virtual ~ThermalExpansionAMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _A;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<RealMatrixX> *_material_expansion;
            MAST::FieldFunction<Real> *_A;
        };
        
        
        
        class ThermalExpansionBMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionBMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                                    MAST::FieldFunction<RealMatrixX> *mat_expansion,
                                    MAST::FieldFunction<Real> *A_y_moment,
                                    MAST::FieldFunction<Real> *A_z_moment);
            
            ThermalExpansionBMatrix(const MAST::Solid1DSectionProperty::ThermalExpansionBMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_material_expansion->master());
                _functions.insert(_A_y_moment->master());
                _functions.insert(_A_z_moment->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::ThermalExpansionBMatrix(*this));
            }
            
            virtual ~ThermalExpansionBMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<RealMatrixX> *_material_expansion;
            MAST::FieldFunction<Real> *_A_y_moment, *_A_z_moment;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                             MAST::FieldFunction<RealMatrixX> *T,
                             MAST::FieldFunction<Real> *A);
            
            PrestressAMatrix(const MAST::Solid1DSectionProperty::PrestressAMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _A(f._A->clone().release()) {
                _functions.insert(_prestress->master());
                _functions.insert(_T->master());
                _functions.insert(_A->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::PrestressAMatrix(*this));
            }
            
            virtual ~PrestressAMatrix() {
                delete _prestress;
                delete _T;
                delete _A;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            //virtual void convert_to_vector(const RealMatrixX& m, DenseRealVector& v) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_prestress, *_T;
            MAST::FieldFunction<Real> *_A;
        };
        
        
        
        class PrestressBMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressBMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                             MAST::FieldFunction<RealMatrixX> *T,
                             MAST::FieldFunction<Real> *A_y_moment,
                             MAST::FieldFunction<Real> *A_z_moment);
            
            PrestressBMatrix(const MAST::Solid1DSectionProperty::PrestressBMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _A_y_moment(f._A_y_moment->clone().release()),
            _A_z_moment(f._A_z_moment->clone().release()) {
                _functions.insert(_prestress->master());
                _functions.insert(_T->master());
                _functions.insert(_A_y_moment->master());
                _functions.insert(_A_z_moment->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid1DSectionProperty::PrestressBMatrix(*this));
            }
            
            virtual ~PrestressBMatrix() {
                delete _prestress;
                delete _T;
                delete _A_y_moment;
                delete _A_z_moment;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            //virtual void convert_to_vector(const RealMatrixX& m, DenseRealVector& v) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_prestress, *_T;
            MAST::FieldFunction<Real> *_A_y_moment, *_A_z_moment;
        };
        
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                                     MAST::FieldFunction<Real> *A);
            
            ThermalConductanceMatrix(const MAST::Solid1DSectionProperty::ThermalConductanceMatrix &f);
            
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const;
            
            virtual ~ThermalConductanceMatrix();
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX>* _mat_cond;
            
            MAST::FieldFunction<Real>* _A;
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                                     MAST::FieldFunction<Real> *h);
            
            ThermalCapacitanceMatrix(const MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix &f);
            
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const;
            
            virtual ~ThermalCapacitanceMatrix();
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX>* _mat_cap;
            
            MAST::FieldFunction<Real>* _h;
        };
    }
    
    
}


MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::
ExtensionStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                         MAST::FieldFunction<Real> *A,
                         MAST::FieldFunction<Real> *J):
MAST::FieldFunction<RealMatrixX> ("ExtensionStiffnessMatrix1D"),
_material_stiffness(mat),
_A(A),
_J(J) {
    _functions.insert(mat->master());
    _functions.insert(A->master());
    _functions.insert(J->master());
}





void
MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::operator() (const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    // [C]*h
    Real A, J;
    (*_A)(p, t, A);
    (*_J)(p, t, J);
    (*_material_stiffness)(p, t, m);
    m.row(0) *= A;
    m.row(1) *= J;
}




void
MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                      const MAST::FunctionBase& f,
                                      const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2);
    Real A, J, dA, dJ;
    (*_A)(p, t, A); _A->derivative(d, f, p, t, dA);
    (*_J)(p, t, J); _J->derivative(d, f, p, t, dJ);
    (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
    
    // [C]*dh
    m.row(0) *= dA;
    m.row(1) *= dJ;
    
    // += [dC]*h
    dm.row(0) *= A;
    dm.row(1) *= J;
    m         += dm;
}






MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix::
ExtensionBendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                MAST::FieldFunction<Real> *A_y_moment,
                                MAST::FieldFunction<Real> *A_z_moment):
MAST::FieldFunction<RealMatrixX> ("ExtensionBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(mat->master());
    _functions.insert(A_y_moment->master());
    _functions.insert(A_z_moment->master());
}



void
MAST::Solid1DSectionProperty::
ExtensionBendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    Real Ay, Az;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_material_stiffness)(p, t, m);
    
    m(0,1)    = m(0,0)*Ay;  // coupling of u and w bending (== theta_y)
    m(0,0)   *= Az;        // coupling of u and v bending (== theta_z)
    
    m.row(1) *= 0; // no coupling for torsion for symmetic sections
}





void
MAST::Solid1DSectionProperty::
ExtensionBendingStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                             const MAST::FunctionBase& f,
                                             const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    RealMatrixX dm;
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->derivative(d, f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->derivative(d, f, p, t, dAz);
    (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
    
    m(0,1)    = m(0,0)*dAy;  // coupling of u and w bending (== theta_y)
    m(0,0)   *= dAz;        // coupling of u and v bending (== theta_z)
    m.row(1) *= 0;     // no coupling for torsion for symmetic sections
    
    dm(0,1)   = dm(0,0)*Ay;
    dm(0,0)  *= Az;
    dm.row(1)*= 0.;
    m        += dm;
}




MAST::Solid1DSectionProperty::BendingStiffnessMatrix::
BendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                       MAST::FieldFunction<RealMatrixX> *I):
MAST::FieldFunction<RealMatrixX> ("BendingStiffnessMatrix1D"),
_material_stiffness(mat),
_I(I) {
    _functions.insert(mat->master());
    _functions.insert(I->master());
}



void
MAST::Solid1DSectionProperty::
BendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat;
    (*_I)(p, t, m);
    (*_material_stiffness)(p, t, mat);
    
    // E*I
    m *= mat(0,0); // scale the inertia matrix with modulus of elasticity
}








void
MAST::Solid1DSectionProperty::
BendingStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat, dmat, dm;
    (*_I)(p, t, m); _I->derivative(d, f, p, t, dm);
    (*_material_stiffness)(p, t, mat); _material_stiffness->derivative(d, f, p, t, dmat);
    
    // dE*I
    m *= dmat(0,0); // scale the inertia matrix with modulus of elasticity
    
    // E*dI
    m += mat(0,0)*dm; // scale the inertia matrix with modulus of elasticity
}





MAST::Solid1DSectionProperty::InertiaMatrix::
InertiaMatrix(MAST::FieldFunction<Real> *rho,
              MAST::FieldFunction<Real> *A,
              MAST::FieldFunction<Real> *A_y_moment,
              MAST::FieldFunction<Real> *A_z_moment,
              MAST::FieldFunction<Real> *Ip,
              MAST::FieldFunction<RealMatrixX> *I):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix1D"),
_rho(rho),
_A(A),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment),
_Ip(Ip),
_I(I) {
    _functions.insert(_rho->master());
    _functions.insert(_A->master());
    _functions.insert(_A_y_moment->master());
    _functions.insert(_A_z_moment->master());
    _functions.insert(_Ip->master());
    _functions.insert(_I->master());
}




void
MAST::Solid1DSectionProperty::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m = RealMatrixX::Zero(6, 6);
    RealMatrixX I;
    Real rho, A, Ay, Az, Ip;
    (*_rho)(p, t, rho);
    (*_A)(p, t, A);
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_Ip)(p, t, Ip);
    (*_I)(p, t, I);
    
    // translation velocities
    m(0,0) = A; m(1,1) = A; m(2,2) = A;
    
    // torsion
    m(3,3) = Ip;
    
    // rotational velocities
    m(0,4) = Ay;  m(4,0) = Ay;   // w-displacement
    m(0,5) = -Az; m(5,0) = -Az;  // v-displacement
    
    // bending rotation inertia
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(4+i,4+j) = I(i,j);
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<3; i++)
        m(i+3,i+3) *= 1.0e-16;
    
    m *= rho;
}








void
MAST::Solid1DSectionProperty::
InertiaMatrix::derivative (const MAST::DerivativeType d,
                           const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(6, 6); dm = RealMatrixX::Zero(6, 6);
    RealMatrixX I, dI;
    Real rho, A, Ay, Az, Ip, drho, dA, dAy, dAz, dIp;
    (*_rho)(p, t, rho); _rho->derivative(d, f, p, t, drho);
    (*_A)(p, t, A); _A->derivative(d, f, p, t, dA);
    (*_A_y_moment)(p, t, Ay); _A_y_moment->derivative(d, f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->derivative(d, f, p, t, dAz);
    (*_Ip)(p, t, Ip); _Ip->derivative(d, f, p, t, dIp);
    (*_I)(p, t, I); _I->derivative(d, f, p, t, dI);
    
    // translation velocities
    m(0,0) = A;  m(1,1) = A;  m(2,2) = A;
    dm(0,0) = dA; dm(1,1) = dA; dm(2,2) = dA;
    
    // torsion
    m(3,3) = Ip;
    dm(3,3) = dIp;
    
    // rotational velocities
    m(0,4) = Ay;  m(4,0) = Ay;   // w-displacement
    dm(0,4) = dAy;  dm(4,0) = dAy;   // w-displacement
    m(0,5) = -Az; m(5,0) = -Az;  // v-displacement
    dm(0,5) = -dAz; m(5,0) = -dAz;  // v-displacement
    
    // bending rotation inertia
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++) {
            m(4+i,4+j) = I(i,j);
            dm(4+i,4+j) = dI(i,j);
        }
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<3; i++) {
        m(i+3,i+3) *= 1.0e-16;
        dm(i+3,i+3) *= 1.0e-16;
    }
    
    m *= drho;
    m += rho*dm;
}




MAST::Solid1DSectionProperty::ThermalExpansionAMatrix::
ThermalExpansionAMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                        MAST::FieldFunction<RealMatrixX> *mat_expansion,
                        MAST::FieldFunction<Real> *A):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionAMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A(A) {
    _functions.insert(mat_stiff->master());
    _functions.insert(mat_expansion->master());
    _functions.insert(_A->master());
}




void
MAST::Solid1DSectionProperty::
ThermalExpansionAMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real A;
    RealMatrixX at;
    (*_A)(p, t, A);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m *= at;
    m *= A;
}





void
MAST::Solid1DSectionProperty::
ThermalExpansionAMatrix::derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real A, dA;
    RealMatrixX m1, at, dat, dm;
    (*_A)(p, t, A); _A->derivative(d, f, p, t, dA);
    (*_material_stiffness)(p, t, m1); _material_stiffness->derivative(d, f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->derivative(d, f, p, t, dat);
    
    m=m1;
    
    m *= at;
    m *= dA;
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m  += A*m1;
}




MAST::Solid1DSectionProperty::ThermalExpansionBMatrix::
ThermalExpansionBMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                        MAST::FieldFunction<RealMatrixX> *mat_expansion,
                        MAST::FieldFunction<Real> *A_y_moment,
                        MAST::FieldFunction<Real> *A_z_moment):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionBMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(mat_stiff->master());
    _functions.insert(mat_expansion->master());
    _functions.insert(_A_y_moment->master());
    _functions.insert(_A_z_moment->master());
}




void
MAST::Solid1DSectionProperty::
ThermalExpansionBMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real Ay, Az;
    RealMatrixX at;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m *= at;
    m(1,0)  = Ay * m(0,0);
    m(0,0) *= Az;
}






void
MAST::Solid1DSectionProperty::
ThermalExpansionBMatrix::derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real Ay, Az, dAy, dAz;
    RealMatrixX at, dat, m1, dm;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->derivative(d, f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->derivative(d, f, p, t, dAz);
    (*_material_stiffness)(p, t, m1); _material_stiffness->derivative(d, f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->derivative(d, f, p, t, dat);
    
    m = m1;
    m *= at;
    m(1,0)  = dAy * m(0,0);
    m(0,0) *= dAz;
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    m1(1,0)  = Ay * m1(0,0);
    m1(0,0) *= Az;
    
    m += m1;
}




MAST::Solid1DSectionProperty::PrestressAMatrix::
PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                 MAST::FieldFunction<RealMatrixX> *T,
                 MAST::FieldFunction<Real> *A):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix1D"),
_prestress(prestress),
_T(T),
_A(A) {
    _functions.insert(prestress->master());
    _functions.insert(T->master());
    _functions.insert(A->master());
}




void
MAST::Solid1DSectionProperty::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real A;
    (*_A)(p, t, A);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    m(0,0) = s(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}





void
MAST::Solid1DSectionProperty::
PrestressAMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real A, dA;
    (*_A)(p, t, A); _A->derivative(d, f, p, t, dA);
    (*_prestress)(p, t, s); _prestress->derivative(d, f, p, t, ds);
    (*_T)(p, t, T); _T->derivative(d, f, p, t, dT);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    RealMatrixX tmp;
    ds * T;
    ds = T.transpose()*ds;
    
    tmp  = T.transpose() * s * dT + dT.transpose() * s * T;
    ds  += tmp;
    
    m(0,0) = s(0,0)*dA + ds(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}



//void
//MAST::Solid1DSectionProperty::
//PrestressAMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v)  const {
//    libmesh_assert_equal_to(m.rows(), 2);
//    libmesh_assert_equal_to(m.cols(), 2);
//    v.resize(2);
//    v(0) = m(0,0);
//}




MAST::Solid1DSectionProperty::PrestressBMatrix::
PrestressBMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                 MAST::FieldFunction<RealMatrixX> *T,
                 MAST::FieldFunction<Real> *A_y_moment,
                 MAST::FieldFunction<Real> *A_z_moment):
MAST::FieldFunction<RealMatrixX>("PrestressBMatrix1D"),
_prestress(prestress),
_T(T),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(prestress->master());
    _functions.insert(T->master());
    _functions.insert(A_y_moment->master());
    _functions.insert(A_z_moment->master());
}




void
MAST::Solid1DSectionProperty::
PrestressBMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real Ay, Az;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s = T.transpose() * s * T;
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  s(0,0)*Az;
    m(0,1) =  s(0,0)*Ay;
}






void
MAST::Solid1DSectionProperty::
PrestressBMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->derivative(d, f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->derivative(d, f, p, t, dAz);
    (*_prestress)(p, t, s); _prestress->derivative(d, f, p, t, ds);
    (*_T)(p, t, T); _T->derivative(d, f, p, t, dT);
    
    // convert the stress to the local coordinate
    s = T.transpose() * s * T;
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    RealMatrixX tmp;
    ds = (T.transpose()  * ds * T +
          T.transpose()  *  s * dT +
          dT.transpose() *  s * T);
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  (s(0,0)*dAz + ds(0,0)*Az);
    m(0,1) =  s(0,0)*dAy + ds(0,0)*Ay;
}



//void
//MAST::Solid1DSectionProperty::
//PrestressBMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v)  const {
//    libmesh_assert_equal_to(m.rows(), 2);
//    libmesh_assert_equal_to(m.cols(), 2);
//    v.resize(2);
//    v(0) = m(0,0);
//    v(1) = m(0,1);
//}




MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                         MAST::FieldFunction<Real> *A):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_A(A) {
    _functions.insert(mat_cond->master());
    _functions.insert(A->master());
}


MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::Solid1DSectionProperty::ThermalConductanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cond(f._mat_cond->clone().release()),
_A(f._A->clone().release()) {
    _functions.insert(_mat_cond->master());
    _functions.insert(_A->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionProperty::ThermalConductanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalConductanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
~ThermalConductanceMatrix() {
    
    delete _mat_cond;
    delete _A;
}


void
MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m = RealMatrixX::Zero(1, 1);
    Real A;
    (*_mat_cond)(p, t, m);
    (*_A)(p, t, A);
    
    m *= A;
}



void
MAST::Solid1DSectionProperty::ThermalConductanceMatrix::derivative (const MAST::DerivativeType d,
                                                                    const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m = RealMatrixX::Zero(1, 1);
    RealMatrixX dm;
    Real A, dA;
    (*_mat_cond)(p, t, m);
    _mat_cond->derivative(d, f, p, t, dm);
    (*_A)(p, t, A);
    _A->derivative(d, f, p, t, dA);
    
    m *= dA;
    m += dm*A;
}




MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cap,
                         MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap),
_h(h) {
    _functions.insert(mat_cap->master());
    _functions.insert(h->master());
}


MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cap(f._mat_cap->clone().release()),
_h(f._h->clone().release()) {
    _functions.insert(_mat_cap->master());
    _functions.insert(_h->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() {
    
    delete _mat_cap;
    delete _h;
}


void
MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m = RealMatrixX::Zero(1, 1);
    Real h;
    (*_mat_cap)(p, t, m);
    (*_h)(p, t, h);
    
    m *= h;
}



void
MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::derivative (const MAST::DerivativeType d,
                                                                    const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m = RealMatrixX::Zero(1, 1);
    RealMatrixX dm;
    Real h, dh;
    (*_mat_cap)(p, t, m);
    _mat_cap->derivative(d, f, p, t, dm);
    (*_h)(p, t, h);
    _h->derivative(d, f, p, t, dh);
    
    m *= dh;
    m += dm*h;
}




bool
MAST::Solid1DSectionElementPropertyCard::
depends_on(const MAST::FunctionBase& f) const {
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}



const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::A() const {
    
    libmesh_assert(_initialized);
    return *_A;
}



const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::J() const {
    
    libmesh_assert(_initialized);
    return *_J;
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Ip() const {
    
    libmesh_assert(_initialized);
    return *_Ip;
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Ay() const {
    
    libmesh_assert(_initialized);
    return *_Ay;
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Az() const {
    
    libmesh_assert(_initialized);
    return *_Az;
}


const MAST::FieldFunction<RealMatrixX>&
MAST::Solid1DSectionElementPropertyCard::I() const {
    
    libmesh_assert(_initialized);
    return *_AI;
}




MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::A() {
    
    libmesh_assert(_initialized);
    return *_A;
}



MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::J() {
    
    libmesh_assert(_initialized);
    return *_J;
}


MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Ip() {
    
    libmesh_assert(_initialized);
    return *_Ip;
}


MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Ay() {
    
    libmesh_assert(_initialized);
    return *_Ay;
}


MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Az() {
    
    libmesh_assert(_initialized);
    return *_Az;
}


MAST::FieldFunction<RealMatrixX>&
MAST::Solid1DSectionElementPropertyCard::I() {
    
    libmesh_assert(_initialized);
    return *_AI;
}





void
MAST::Solid1DSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DSectionProperty::Area(hy.clone().release(),
                                                    hz.clone().release()));
    _Ay.reset(new MAST::Solid1DSectionProperty::AreaYMoment(hy.clone().release(),
                                                            hz.clone().release(),
                                                            hy_off.clone().release()));
    _Az.reset(new MAST::Solid1DSectionProperty::AreaZMoment(hy.clone().release(),
                                                            hz.clone().release(),
                                                            hz_off.clone().release()));
    _J.reset(new MAST::Solid1DSectionProperty::TorsionalConstant(hy.clone().release(),
                                                                 hz.clone().release()));
    _Ip.reset(new MAST::Solid1DSectionProperty::PolarInertia(hy.clone().release(),
                                                             hz.clone().release(),
                                                             hy_off.clone().release(),
                                                             hz_off.clone().release()));
    _AI.reset(new MAST::Solid1DSectionProperty::AreaInertiaMatrix(hy.clone().release(),
                                                                  hz.clone().release(),
                                                                  hy_off.clone().release(),
                                                                  hz_off.clone().release()));
    
    _initialized = true;
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(1).release(),
     _A->clone().release(),
     _J->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(1).release(),
     _Ay->clone().release(),
     _Az->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(1).release(),
     _AI->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_error();
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (NULL);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
inertia_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho").clone().release(),
     _A->clone().release(),
     _Ay->clone().release(),
     _Az->clone().release(),
     _Ip->clone().release(),
     _AI->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(1).release(),
     _material->thermal_expansion_matrix(1).release(),
     _A->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(1).release(),
     _material->thermal_expansion_matrix(1).release(),
     _Ay->clone().release(),
     _Az->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix().release(),
     _A->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_A_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::PrestressAMatrix
    (this->get<MAST::FieldFunction<RealMatrixX> >("prestress").clone().release(),
     e.local_elem().T_matrix_function().release(),
     _A->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_B_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::PrestressBMatrix
    (this->get<MAST::FieldFunction<RealMatrixX> >("prestress").clone().release(),
     e.local_elem().T_matrix_function().release(),
     _Ay->clone().release(),
     _Az->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(1).release(),
     _A->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(1).release(),
     _A->clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



