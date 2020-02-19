/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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



namespace MAST {
    namespace Solid1DSectionProperty {
        
        class Area: public MAST::FieldFunction<Real> {
        public:
            Area(const MAST::FieldFunction<Real>& hy,
                 const MAST::FieldFunction<Real>&  hz):
            MAST::FieldFunction<Real>("Area"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(&hy);
                _functions.insert(&hz);
            }
            
            virtual ~Area() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz;
                _hy(p, t, hy);
                _hz(p, t, hz);
                
                m = hy*hz;
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz;
                _hy(p, t, hy); _hy.derivative( f, p, t, dhy);
                _hz(p, t, hz); _hz.derivative( f, p, t, dhz);
                
                m = dhy*hz + hy*dhz;
            }
            
        protected:
            
            const MAST::FieldFunction<Real>& _hy, &_hz;
        };
        
        
        
        class TorsionalConstant: public MAST::FieldFunction<Real> {
        public:
            TorsionalConstant(const MAST::FieldFunction<Real>& hy,
                              const MAST::FieldFunction<Real>&  hz):
            MAST::FieldFunction<Real>("TorsionalConstant"),
            _hy(hy),
            _hz(hz) {
                _functions.insert(&hy);
                _functions.insert(&hz);
            }
            
            
            virtual ~TorsionalConstant() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, a, b;
                _hy(p, t, hy);
                _hz(p, t, hz);
                
                
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
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz, a, b, da, db;
                _hy(p, t, hy); _hy.derivative( f, p, t, dhy);
                _hz(p, t, hz); _hz.derivative( f, p, t, dhz);
                
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
            
            const MAST::FieldFunction<Real>& _hy, &_hz;
        };
        
        
        
        class PolarInertia: public MAST::FieldFunction<Real> {
        public:
            PolarInertia(const MAST::FieldFunction<Real>& hy,
                         const MAST::FieldFunction<Real>&  hz,
                         const MAST::FieldFunction<Real>&  hy_offset,
                         const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("PolarInertia"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&hy);
                _functions.insert(&hz);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~PolarInertia() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, offy, offz;
                _hy(p, t, hy);
                _hz(p, t, hz);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m = hy*hz*((pow(hy,2) + pow(hz,2))/12. +
                           pow(offy,2) + pow(offz,2));
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, dhy, dhz, offy, offz, doffy, doffz;
                _hy        (p, t, hy);           _hy.derivative( f, p, t, dhy);
                _hz        (p, t, hz);           _hz.derivative( f, p, t, dhz);
                _hy_offset (p, t, offy);  _hy_offset.derivative( f, p, t, doffy);
                _hz_offset (p, t, offz);  _hz_offset.derivative( f, p, t, doffz);
                
                
                m =
                (dhy*hz + hy*dhz) * ((pow(hy,2) + pow(hz,2))/12. +
                                     pow(offy,2) + pow(offz,2)) +
                2.*hy*hz*((hy*dhy + hz*dhz)/12. +
                          offy*doffy + offz*doffz);
            }
            
        protected:
            
            const MAST::FieldFunction<Real>& _hy, &_hz, &_hy_offset, &_hz_offset;
        };
        
        
        
        /*!
         *   calculates the area moment about the Y-axis due to an offset 
         *   along the Z-axis
         */
        class AreaYMoment: public MAST::FieldFunction<Real> {
        public:
            AreaYMoment(const MAST::FieldFunction<Real>&  hy,
                        const MAST::FieldFunction<Real>&  hz,
                        const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<Real>("AreaYMoment"),
            _hy(hy),
            _hz(hz),
            _hz_offset(hz_offset) {
                _functions.insert(&hy);
                _functions.insert(&hz);
                _functions.insert(&hz_offset);
            }
            
            
            virtual ~AreaYMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off;
                _hy(p, t, hy);
                _hz(p, t, hz);
                _hz_offset(p, t, off);
                
                m = hy*hz*off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                _hy        (p, t, hy);         _hy.derivative( f, p, t, dhy);
                _hz        (p, t, hz);         _hz.derivative( f, p, t, dhz);
                _hz_offset (p, t, off); _hz_offset.derivative( f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            const MAST::FieldFunction<Real>& _hy, &_hz, &_hz_offset;
        };
        
        
        
        /*!
         *   calculates the area moment about the Z-axis due to an offset
         *   along the Y-axis
         */
        class AreaZMoment: public MAST::FieldFunction<Real> {
        public:
            AreaZMoment(const MAST::FieldFunction<Real>&  hy,
                        const MAST::FieldFunction<Real>&  hz,
                        const MAST::FieldFunction<Real>&  hy_offset):
            MAST::FieldFunction<Real>("AreaZMoment"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset) {
                _functions.insert(&hy);
                _functions.insert(&hz);
                _functions.insert(&hy_offset);
            }
            
            virtual ~AreaZMoment() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off;
                _hy(p, t, hy);
                _hz(p, t, hz);
                _hy_offset(p, t, off);
                
                m = hy*hz*off;
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const {
                Real hy, hz, off, dhy, dhz, doff;
                _hy(p, t, hy); _hy.derivative( f, p, t, dhy);
                _hz(p, t, hz); _hz.derivative( f, p, t, dhz);
                _hy_offset(p, t, off); _hy_offset.derivative( f, p, t, doff);
                
                m = dhy*hz*off + hy*dhz*off + hy*hz*doff;
            }
            
        protected:
            
            const MAST::FieldFunction<Real>& _hy, &_hz, &_hy_offset;
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
            AreaInertiaMatrix(const MAST::FieldFunction<Real>&  hy,
                              const MAST::FieldFunction<Real>&  hz,
                              const MAST::FieldFunction<Real>&  hy_offset,
                              const MAST::FieldFunction<Real>&  hz_offset):
            MAST::FieldFunction<RealMatrixX>("AreaInertiaMatrix"),
            _hy(hy),
            _hz(hz),
            _hy_offset(hy_offset),
            _hz_offset(hz_offset) {
                _functions.insert(&hy);
                _functions.insert(&hz);
                _functions.insert(&hy_offset);
                _functions.insert(&hz_offset);
            }
            
            virtual ~AreaInertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real hy, hz, offy, offz;
                m = RealMatrixX::Zero(2,2);
                _hy(p, t, hy);
                _hz(p, t, hz);
                _hy_offset(p, t, offy);
                _hz_offset(p, t, offz);
                
                m(0,0) = hz*pow(hy,3)/12. + hy*hz*pow(offy,2) ; // Izz for v-bending
                m(0,1) = hy*hz*offy*offz;
                m(1,0) = m(0,1);
                m(1,1) = hy*pow(hz,3)/12. + hy*hz*pow(offz,2) ; // Iyy for w-bending
            }
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real hy, hz, offy, offz, dhy, dhz, doffy, doffz;
                m = RealMatrixX::Zero(2,2);
                _hy(p, t, hy); _hy.derivative( f, p, t, dhy);
                _hz(p, t, hz); _hz.derivative( f, p, t, dhz);
                _hy_offset(p, t, offy); _hy_offset.derivative( f, p, t, doffy);
                _hz_offset(p, t, offz); _hz_offset.derivative( f, p, t, doffz);
                
                
                m(0,0) = dhz*pow(hy,3)/12. + hz*pow(hy,2)/4.*dhy +
                dhy*hz*pow(offy,2) + hy*dhz*pow(offy,2) + 2.*hy*hz*offy*doffy ;
                m(0,1) = dhy*hz*offy*offz + hy*dhz*offy*offz +
                hy*hz*doffy*offz + hy*hz*offy*doffz;
                m(1,0) = m(0,1);
                m(1,1) = dhy*pow(hz,3)/12. + hy*pow(hz,2)/4.*dhz +
                dhy*hz*pow(offz,2) + hy*dhz*pow(offz,2) + 2.*hy*hz*offz*doffz ;
            }
            
        protected:
            
            const MAST::FieldFunction<Real>& _hy, &_hz, &_hy_offset, &_hz_offset;
        };
        
        
        class WarpingConstant: public MAST::FieldFunction<Real> 
        {
        public:
            WarpingConstant():
            MAST::FieldFunction<Real>("WarpingConstant")
            {
            }
            
            virtual ~WarpingConstant()
            {
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const
            {
                m = 0.0;
            }
            
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& m) const
            {
                m = 0.0;
            }
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
        
        
        class ExtensionStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                     const MAST::FieldFunction<Real>&  A,
                                     const MAST::FieldFunction<Real>&  J);
            
            virtual ~ExtensionStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<Real>& _A, &_J;
        };
        
        
        
        class ExtensionBendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionBendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                            const MAST::FieldFunction<Real>&  A_y_moment,
                                            const MAST::FieldFunction<Real>&  A_z_moment);
            
            virtual ~ExtensionBendingStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<Real>& _A_y_moment, &_A_z_moment;
        };
        
        
        class BendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            BendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                   const MAST::FieldFunction<RealMatrixX>& I);
            
            virtual ~BendingStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<RealMatrixX>& _I;
        };
        
        
        class TransverseStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                      const MAST::FieldFunction<Real>&  A,
                                      const MAST::FieldFunction<RealMatrixX>& Kappa):
            MAST::FieldFunction<RealMatrixX>("TransverseStiffnessMatrix1D"),
            _material_stiffness(mat),
            _A(A),
            _Kappa(Kappa)
            {
                _functions.insert(&mat);
                _functions.insert(&A);
                _functions.insert(&Kappa);
            }
            
            virtual ~TransverseStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real A;
                RealMatrixX Kappa;
                _A(p, t, A);
                _Kappa(p, t, Kappa);
                _material_stiffness(p, t, m);

                m *= A;
                
                m(0,0) *= Kappa(0,0);
                m(1,1) *= Kappa(1,1);
                m(0,1) *= Kappa(0,1);
                m(1,0) *= Kappa(1,0);
                
                /**
                 * m = [G*kappa*A,          0,
                 *              0,   G*kappa*A]
                 */
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                RealMatrixX G, dG, dm, Kappa, dKappa;
                Real A, dA;
                _A                  (p, t, A);                  _A.derivative( f, p, t, dA);
                _material_stiffness (p, t, G); _material_stiffness.derivative( f, p, t, dG);
                _Kappa              (p, t, Kappa);
                _Kappa.derivative(f, p, t, dKappa);
                
                m(0,0) = (dG(0,0) * Kappa(0,0) * A) + (G(0,0) * dKappa(0,0) * A) + (G(0,0) * Kappa(0,0) * dA);
                m(1,1) = (dG(1,1) * Kappa(1,1) * A) + (G(1,1) * dKappa(1,1) * A) + (G(1,1) * Kappa(1,1) * dA);
                m(0,1) = (dG(0,1) * Kappa(0,1) * A) + (G(0,1) * dKappa(0,1) * A) + (G(0,1) * Kappa(0,1) * dA);
                m(1,0) = (dG(1,0) * Kappa(1,0) * A) + (G(1,0) * dKappa(1,0) * A) + (G(1,0) * Kappa(1,0) * dA);
            }
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<Real>& _A;
            const MAST::FieldFunction<RealMatrixX>& _Kappa;
        };
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(const MAST::FieldFunction<Real>&  rho,
                          const MAST::FieldFunction<Real>&  A,
                          const MAST::FieldFunction<Real>&  A_y_moment,
                          const MAST::FieldFunction<Real>&  A_z_moment,
                          const MAST::FieldFunction<Real>&  Ip,
                          const MAST::FieldFunction<RealMatrixX>& I);
            
            virtual ~InertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _rho, &_A, &_A_y_moment, &_A_z_moment, &_Ip;
            const MAST::FieldFunction<RealMatrixX>& _I;
        };
        
        
        
        class ThermalExpansionAMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionAMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                                    const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                                    const MAST::FieldFunction<Real>& A);
            
            virtual ~ThermalExpansionAMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<RealMatrixX>& _material_expansion;
            const MAST::FieldFunction<Real>& _A;
        };
        
        
        
        class ThermalExpansionBMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionBMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                                    const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                                    const MAST::FieldFunction<Real>& A_y_moment,
                                    const MAST::FieldFunction<Real>& A_z_moment);
            
            virtual ~ThermalExpansionBMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<RealMatrixX>& _material_expansion;
            const MAST::FieldFunction<Real>& _A_y_moment, &_A_z_moment;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                             const MAST::FieldFunction<RealMatrixX>& T,
                             const MAST::FieldFunction<Real>& A);
            
            virtual ~PrestressAMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            //virtual void convert_to_vector(const RealMatrixX& m, DenseRealVector& v) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _prestress, &_T;
            const MAST::FieldFunction<Real>& _A;
        };
        
        
        
        class PrestressBMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressBMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                             const MAST::FieldFunction<RealMatrixX>& T,
                             const MAST::FieldFunction<Real>& A_y_moment,
                             const MAST::FieldFunction<Real>& A_z_moment);
            
            virtual ~PrestressBMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            //virtual void convert_to_vector(const RealMatrixX& m, DenseRealVector& v) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _prestress, &_T;
            const MAST::FieldFunction<Real>& _A_y_moment, &_A_z_moment;
        };
        
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                                     const MAST::FieldFunction<Real>& A);
            
            virtual ~ThermalConductanceMatrix();
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _mat_cond;
            
            const MAST::FieldFunction<Real>&  _A;
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                                     const MAST::FieldFunction<Real>& h);
            
            virtual ~ThermalCapacitanceMatrix();
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _mat_cap;
            
            const MAST::FieldFunction<Real>&  _h;
        };
        
        
        
    }
    
    
}


MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::
ExtensionStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                         const MAST::FieldFunction<Real>& A,
                         const MAST::FieldFunction<Real>& J):
MAST::FieldFunction<RealMatrixX> ("ExtensionStiffnessMatrix1D"),
_material_stiffness(mat),
_A(A),
_J(J) {
    _functions.insert(&mat);
    _functions.insert(&A);
    _functions.insert(&J);
}





void
MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::operator() (const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    // [C]*h
    Real A, J;
    _A(p, t, A);
    _J(p, t, J);
    _material_stiffness(p, t, m);
    m.row(0) *= A;
    m.row(1) *= J;
}




void
MAST::Solid1DSectionProperty::
ExtensionStiffnessMatrix::derivative (     const MAST::FunctionBase& f,
                                      const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2);
    Real A, J, dA, dJ;
    _A(p, t, A); _A.derivative( f, p, t, dA);
    _J(p, t, J); _J.derivative( f, p, t, dJ);
    _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
    
    // [C]*dh
    m.row(0) *= dA;
    m.row(1) *= dJ;
    
    // += [dC]*h
    dm.row(0) *= A;
    dm.row(1) *= J;
    m         += dm;
}






MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix::
ExtensionBendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                const MAST::FieldFunction<Real>& A_y_moment,
                                const MAST::FieldFunction<Real>& A_z_moment):
MAST::FieldFunction<RealMatrixX> ("ExtensionBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(&mat);
    _functions.insert(&A_y_moment);
    _functions.insert(&A_z_moment);
}



void
MAST::Solid1DSectionProperty::
ExtensionBendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    Real Ay, Az;
    _A_y_moment(p, t, Ay);
    _A_z_moment(p, t, Az);
    _material_stiffness(p, t, m);
    
    m(0,1)    = m(0,0)*Ay;  // coupling of u and w bending (== theta_y)
    m(0,0)   *= Az;        // coupling of u and v bending (== theta_z)
    
    m.row(1) *= 0; // no coupling for torsion for symmetic sections
}





void
MAST::Solid1DSectionProperty::
ExtensionBendingStiffnessMatrix::derivative (            const MAST::FunctionBase& f,
                                             const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    RealMatrixX dm;
    Real Ay, Az, dAy, dAz;
    _A_y_moment(p, t, Ay); _A_y_moment.derivative( f, p, t, dAy);
    _A_z_moment(p, t, Az); _A_z_moment.derivative( f, p, t, dAz);
    _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
    
    m(0,1)    = m(0,0)*dAy;  // coupling of u and w bending (== theta_y)
    m(0,0)   *= dAz;        // coupling of u and v bending (== theta_z)
    m.row(1) *= 0;     // no coupling for torsion for symmetic sections
    
    dm(0,1)   = dm(0,0)*Ay;
    dm(0,0)  *= Az;
    dm.row(1)*= 0.;
    m        += dm;
}




MAST::Solid1DSectionProperty::BendingStiffnessMatrix::
BendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                       const MAST::FieldFunction<RealMatrixX>& I):
MAST::FieldFunction<RealMatrixX> ("BendingStiffnessMatrix1D"),
_material_stiffness(mat),
_I(I) {
    _functions.insert(&mat);
    _functions.insert(&I);
}



void
MAST::Solid1DSectionProperty::
BendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat;
    _I(p, t, m);
    _material_stiffness(p, t, mat);
    
    // E*I
    m *= mat(0,0); // scale the inertia matrix with modulus of elasticity
}








void
MAST::Solid1DSectionProperty::
BendingStiffnessMatrix::derivative (   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat, dmat, dm;
    _I(p, t, m); _I.derivative( f, p, t, dm);
    _material_stiffness(p, t, mat); _material_stiffness.derivative( f, p, t, dmat);
    
    // dE*I
    m *= dmat(0,0); // scale the inertia matrix with modulus of elasticity
    
    // E*dI
    m += mat(0,0)*dm; // scale the inertia matrix with modulus of elasticity
}





MAST::Solid1DSectionProperty::InertiaMatrix::
InertiaMatrix(const MAST::FieldFunction<Real>& rho,
              const MAST::FieldFunction<Real>& A,
              const MAST::FieldFunction<Real>& A_y_moment,
              const MAST::FieldFunction<Real>& A_z_moment,
              const MAST::FieldFunction<Real>& Ip,
              const MAST::FieldFunction<RealMatrixX>& I):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix1D"),
_rho(rho),
_A(A),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment),
_Ip(Ip),
_I(I) {
    _functions.insert(&_rho);
    _functions.insert(&_A);
    _functions.insert(&_A_y_moment);
    _functions.insert(&_A_z_moment);
    _functions.insert(&_Ip);
    _functions.insert(&_I);
}




void
MAST::Solid1DSectionProperty::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m = RealMatrixX::Zero(6, 6);
    RealMatrixX I;
    Real rho, A, Ay, Az, Ip;
    _rho(p, t, rho);
    _A(p, t, A);
    _A_y_moment(p, t, Ay);
    _A_z_moment(p, t, Az);
    _Ip(p, t, Ip);
    _I(p, t, I);
    
    // translation velocities
    m(0,0) = A; m(1,1) = A; m(2,2) = A;
    
    // torsion
    m(3,3) = Ip;
    
    // rotational velocities
    // theta-y rotation
    m(0,4) =  Ay; m(4,0) =  Ay;
    m(4,4) =  I(1,1); // I11 is defined about y-y axis for theta-y

    // theta-z rotation
    m(0,5) =  Az; m(5,0) =  Az;
    m(5,5) =  I(0,0); // I00 is defined about z-z axis for theta-z

    m(4,5) =  m(5,4) = I(0,1);

    m *= rho;
}








void
MAST::Solid1DSectionProperty::
InertiaMatrix::derivative (               const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(6, 6); dm = RealMatrixX::Zero(6, 6);
    RealMatrixX I, dI;
    Real rho, A, Ay, Az, Ip, drho, dA, dAy, dAz, dIp;
    _rho(p, t, rho); _rho.derivative( f, p, t, drho);
    _A(p, t, A); _A.derivative( f, p, t, dA);
    _A_y_moment(p, t, Ay); _A_y_moment.derivative( f, p, t, dAy);
    _A_z_moment(p, t, Az); _A_z_moment.derivative( f, p, t, dAz);
    _Ip(p, t, Ip); _Ip.derivative( f, p, t, dIp);
    _I(p, t, I); _I.derivative( f, p, t, dI);
    
    // translation velocities
    m(0,0) = A;  m(1,1) = A;  m(2,2) = A;
    dm(0,0) = dA; dm(1,1) = dA; dm(2,2) = dA;
    
    // torsion
    m(3,3) = Ip;
    dm(3,3) = dIp;
    
    // rotational velocities
    // theta-y rotation
    m(0,4) =  Ay; m(4,0) =  Ay;
    m(4,4) =  I(1,1);
    dm(0,4) = dAy;  dm(4,0) = dAy;
    dm(4,4) =  dI(1,1);

    // theta-z rotation
    m(0,5) =  Az; m(5,0) =  Az;
    m(5,5) =  I(0,0);
    dm(0,5) = dAz; dm(5,0) = dAz;  // v-displacement
    dm(5,5) = dI(0,0);

    m(4,5) =  m(5,4) = I(0,1);
    dm(4,5) =  dm(5,4) = dI(0,1);

    m *= drho;
    m += rho*dm;
}




MAST::Solid1DSectionProperty::ThermalExpansionAMatrix::
ThermalExpansionAMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                        const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                        const MAST::FieldFunction<Real>& A):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionAMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A(A) {
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
    _functions.insert(&_A);
}




void
MAST::Solid1DSectionProperty::
ThermalExpansionAMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real A;
    RealMatrixX at;
    _A(p, t, A);
    _material_stiffness(p, t, m);
    _material_expansion(p, t, at);
    
    m *= at;
    m *= A;
}





void
MAST::Solid1DSectionProperty::
ThermalExpansionAMatrix::derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real A, dA;
    RealMatrixX m1, at, dat, dm;
    _A(p, t, A); _A.derivative( f, p, t, dA);
    _material_stiffness(p, t, m1); _material_stiffness.derivative( f, p, t, dm);
    _material_expansion(p, t, at); _material_expansion.derivative( f, p, t, dat);
    
    m=m1;
    
    m *= at;
    m *= dA;
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m  += A*m1;
}




MAST::Solid1DSectionProperty::ThermalExpansionBMatrix::
ThermalExpansionBMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                        const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                        const MAST::FieldFunction<Real>& A_y_moment,
                        const MAST::FieldFunction<Real>& A_z_moment):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionBMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
    _functions.insert(&_A_y_moment);
    _functions.insert(&_A_z_moment);
}




void
MAST::Solid1DSectionProperty::
ThermalExpansionBMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real Ay, Az;
    RealMatrixX at;
    _A_y_moment(p, t, Ay);
    _A_z_moment(p, t, Az);
    _material_stiffness(p, t, m);
    _material_expansion(p, t, at);
    
    m *= at;
    m(1,0)  = Ay * m(0,0); // for w-displacement, area moment about Y-axis
    m(0,0) *= Az;          // for v-displacement, area moment about Z-axis
}






void
MAST::Solid1DSectionProperty::
ThermalExpansionBMatrix::derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    Real Ay, Az, dAy, dAz;
    RealMatrixX at, dat, m1, dm;
    _A_y_moment(p, t, Ay); _A_y_moment.derivative( f, p, t, dAy);
    _A_z_moment(p, t, Az); _A_z_moment.derivative( f, p, t, dAz);
    _material_stiffness(p, t, m1); _material_stiffness.derivative( f, p, t, dm);
    _material_expansion(p, t, at); _material_expansion.derivative( f, p, t, dat);
    
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
PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                 const MAST::FieldFunction<RealMatrixX>& T,
                 const MAST::FieldFunction<Real>& A):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix1D"),
_prestress(prestress),
_T(T),
_A(A) {
    _functions.insert(&prestress);
    _functions.insert(&T);
    _functions.insert(&A);
}




void
MAST::Solid1DSectionProperty::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real A;
    _A(p, t, A);
    _prestress(p, t, s);
    _T(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    m(0,0) = s(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}





void
MAST::Solid1DSectionProperty::
PrestressAMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real A, dA;
    _A(p, t, A); _A.derivative( f, p, t, dA);
    _prestress(p, t, s); _prestress.derivative( f, p, t, ds);
    _T(p, t, T); _T.derivative( f, p, t, dT);
    
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




MAST::Solid1DSectionProperty::PrestressBMatrix::
PrestressBMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                 const MAST::FieldFunction<RealMatrixX>& T,
                 const MAST::FieldFunction<Real>& A_y_moment,
                 const MAST::FieldFunction<Real>& A_z_moment):
MAST::FieldFunction<RealMatrixX>("PrestressBMatrix1D"),
_prestress(prestress),
_T(T),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(&prestress);
    _functions.insert(&T);
    _functions.insert(&A_y_moment);
    _functions.insert(&A_z_moment);
}




void
MAST::Solid1DSectionProperty::
PrestressBMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real Ay, Az;
    _A_y_moment(p, t, Ay);
    _A_z_moment(p, t, Az);
    _prestress(p, t, s);
    _T(p, t, T);
    
    // convert the stress to the local coordinate
    s = T.transpose() * s * T;
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  s(0,0)*Az;
    m(0,1) =  s(0,0)*Ay;
}






void
MAST::Solid1DSectionProperty::
PrestressBMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real Ay, Az, dAy, dAz;
    _A_y_moment(p, t, Ay); _A_y_moment.derivative( f, p, t, dAy);
    _A_z_moment(p, t, Az); _A_z_moment.derivative( f, p, t, dAz);
    _prestress(p, t, s); _prestress.derivative( f, p, t, ds);
    _T(p, t, T); _T.derivative( f, p, t, dT);
    
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






MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                         const MAST::FieldFunction<Real>& A):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_A(A) {
    _functions.insert(&mat_cond);
    _functions.insert(&A);
}


MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
~ThermalConductanceMatrix() { }


void
MAST::Solid1DSectionProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m = RealMatrixX::Zero(1, 1);
    Real A;
    _mat_cond(p, t, m);
    _A(p, t, A);
    
    m *= A;
}



void
MAST::Solid1DSectionProperty::ThermalConductanceMatrix::derivative (                                   const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m = RealMatrixX::Zero(1, 1);
    RealMatrixX dm;
    Real A, dA;
    _mat_cond(p, t, m);
    _mat_cond.derivative( f, p, t, dm);
    _A(p, t, A);
    _A.derivative( f, p, t, dA);
    
    m *= dA;
    m += dm*A;
}




MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cap,
                         const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap),
_h(h) {
    _functions.insert(&mat_cap);
    _functions.insert(&h);
}




MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() { }


void
MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m = RealMatrixX::Zero(1, 1);
    Real h;
    _mat_cap(p, t, m);
    _h(p, t, h);
    
    m *= h;
}



void
MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix::derivative (                                   const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m = RealMatrixX::Zero(1, 1);
    RealMatrixX dm;
    Real h, dh;
    _mat_cap(p, t, m);
    _mat_cap.derivative( f, p, t, dm);
    _h(p, t, h);
    _h.derivative( f, p, t, dh);
    
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


const MAST::FieldFunction<RealMatrixX>&
MAST::Solid1DSectionElementPropertyCard::Kap() const {
    
    libmesh_assert(_initialized);
    return *_Kappa;
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Gam() const {
    
    libmesh_assert(_initialized);
    return *_Gamma;
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


MAST::FieldFunction<RealMatrixX>&
MAST::Solid1DSectionElementPropertyCard::Kap() {
    
    libmesh_assert(_initialized);
    return *_Kappa;
}


MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::Gam() {
    
    libmesh_assert(_initialized);
    return *_Gamma;
}


void
MAST::Solid1DSectionElementPropertyCard::clear() {

    _A.reset();
    _Ay.reset();
    _Az.reset();
    _J.reset();
    _Ip.reset();
    _AI.reset();
    _Gamma.reset();
    _Kappa.reset();
    
    _initialized = false;
}



void
MAST::Solid1DSectionElementPropertyCard::init() {
    
    libmesh_assert(!_initialized);
    
    MAST::FieldFunction<Real>
    &hy         =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz         =  this->get<MAST::FieldFunction<Real> >("hz"),
    &Kappazz    =  this->get<MAST::FieldFunction<Real> >("Kappazz"),
    &Kappayy    =  this->get<MAST::FieldFunction<Real> >("Kappayy"),
    &hy_off     =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off     =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    _A.reset(new MAST::Solid1DSectionProperty::Area(hy,
                                                    hz));
    _Ay.reset(new MAST::Solid1DSectionProperty::AreaYMoment(hy,
                                                            hz,
                                                            hz_off));
    _Az.reset(new MAST::Solid1DSectionProperty::AreaZMoment(hy,
                                                            hz,
                                                            hy_off));
    _J.reset(new MAST::Solid1DSectionProperty::TorsionalConstant(hy,
                                                                 hz));
    _Ip.reset(new MAST::Solid1DSectionProperty::PolarInertia(hy,
                                                             hz,
                                                             hy_off,
                                                             hz_off));
    _AI.reset(new MAST::Solid1DSectionProperty::AreaInertiaMatrix(hy,
                                                                  hz,
                                                                  hy_off,
                                                                  hz_off));
    
    _Kappa.reset(new MAST::Solid1DSectionProperty::ShearCoefficientMatrix(Kappazz, Kappayy));
    
    _Gamma.reset(new MAST::Solid1DSectionProperty::WarpingConstant());
    
    _initialized = true;
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(1), *_A, *_J);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_A_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(1), *_A, *_J);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(1), *_Ay, *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_B_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(1), *_Ay, *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(1), *_AI);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
stiffness_D_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(1), *_AI);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_error();
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
damping_matrix() const {
    
    libmesh_error();
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
inertia_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho"),
     *_A,
     *_Ay,
     *_Az,
     *_Ip,
     *_AI);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
inertia_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho"),
     *_A,
     *_Ay,
     *_Az,
     *_Ip,
     *_AI);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(1),
     _material->thermal_expansion_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_A_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(1),
     _material->thermal_expansion_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(1),
     _material->thermal_expansion_matrix(1),
     *_Ay,
     *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_expansion_B_matrix() const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(1),
     _material->thermal_expansion_matrix(1),
     *_Ay,
     *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix(),
     *_A, *_Kappa);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
transverse_shear_stiffness_matrix() const {
    
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix(),
     *_A, *_Kappa);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_A_matrix(MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid1DSectionProperty::PrestressAMatrix
    //(this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_A_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid1DSectionProperty::PrestressAMatrix
    //(this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_B_matrix(MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid1DSectionProperty::PrestressBMatrix
    //(this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // *_Ay,
    // *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
prestress_B_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);

    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid1DSectionProperty::PrestressBMatrix
    //(this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // *_Ay,
    // *_Az);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_conductance_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid1DSectionElementPropertyCard::
thermal_capacitance_matrix() const {
    
    // make sure that the init method has been called on the card
    libmesh_assert(_initialized);
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid1DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(1),
     *_A);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::
section(const MAST::ElementBase& e) const {
    
    return *_A;
}


const MAST::FieldFunction<Real>&
MAST::Solid1DSectionElementPropertyCard::
section() const {
    
    return *_A;
}


const libMesh::Point 
MAST::Solid1DSectionElementPropertyCard::
get_shear_center(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    hy(p, t, h_y);
    hz(p, t, h_z);
    hy_off(p, t, h_y_off);
    hz_off(p, t, h_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    
    libMesh::Point ps(0., 0.);
    return  ps+offset;
}


const libMesh::Point 
MAST::Solid1DSectionElementPropertyCard::
get_shear_center_derivative(MAST::FunctionBase& f, const libMesh::Point& p, 
                            const Real t)
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    Real dh_y, dh_z, dh_y_off, dh_z_off;
    
    hy(p, t, h_y);          hy.derivative(f, p, t, dh_y);
    hz(p, t, h_z);          hz.derivative(f, p, t, dh_z);
    hy_off(p, t, h_y_off);  hy_off.derivative(f, p, t, dh_y_off);
    hz_off(p, t, h_z_off);  hz_off.derivative(f, p, t, dh_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    libMesh::Point doffset(dh_z_off, dh_y_off);
    
    return  doffset;
}


const libMesh::Point 
MAST::Solid1DSectionElementPropertyCard::
get_centroid(const libMesh::Point& p, const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    hy(p, t, h_y);
    hz(p, t, h_z);
    hy_off(p, t, h_y_off);
    hz_off(p, t, h_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    
    libMesh::Point ps(0., 0.);
    return  ps+offset;
}


const libMesh::Point 
MAST::Solid1DSectionElementPropertyCard::
get_centroid_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, 
                        const Real t) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    Real dh_y, dh_z, dh_y_off, dh_z_off;
    
    hy(p, t, h_y);          hy.derivative(f, p, t, dh_y);
    hz(p, t, h_z);          hz.derivative(f, p, t, dh_z);
    hy_off(p, t, h_y_off);  hy_off.derivative(f, p, t, dh_y_off);
    hz_off(p, t, h_z_off);  hz_off.derivative(f, p, t, dh_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    libMesh::Point doffset(dh_z_off, dh_y_off);
    
    return  doffset;
}


const std::vector<libMesh::Point> 
MAST::Solid1DSectionElementPropertyCard::
get_geom_points(const libMesh::Point& p, const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    hy(p, t, h_y);
    hz(p, t, h_z);
    hy_off(p, t, h_y_off);
    hz_off(p, t, h_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    
    std::vector<libMesh::Point> geom_points = {
        libMesh::Point(-0.5*h_z, -0.5*h_y) + offset,
        libMesh::Point(0.5*h_z, -0.5*h_y) + offset,
        libMesh::Point(0.5*h_z, 0.5*h_y) + offset,
        libMesh::Point(-0.5*h_z, 0.5*h_y) + offset
    };
    
    return geom_points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DSectionElementPropertyCard::
get_geom_points_derivative(const MAST::FunctionBase& f, const libMesh::Point& p, 
                           const Real t, const uint n) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    Real dh_y, dh_z, dh_y_off, dh_z_off;
    
    hy(p, t, h_y);          hy.derivative(f, p, t, dh_y);
    hz(p, t, h_z);          hz.derivative(f, p, t, dh_z);
    hy_off(p, t, h_y_off);  hy_off.derivative(f, p, t, dh_y_off);
    hz_off(p, t, h_z_off);  hz_off.derivative(f, p, t, dh_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    libMesh::Point doffset(dh_z_off, dh_y_off);
    
    std::vector<libMesh::Point> geom_points = {
        libMesh::Point(-0.5*dh_z, -0.5*dh_y) + doffset,
        libMesh::Point(0.5*dh_z, -0.5*dh_y) + doffset,
        libMesh::Point(0.5*dh_z, 0.5*dh_y) + doffset,
        libMesh::Point(-0.5*dh_z, 0.5*dh_y) + doffset
    };
    
    return geom_points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DSectionElementPropertyCard::
get_stress_points(const libMesh::Point& p, const Real t, const libMesh::Point ps) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    hy(p, t, h_y);
    hz(p, t, h_z);
    hy_off(p, t, h_y_off);
    hz_off(p, t, h_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
                
    std::vector<libMesh::Point> stress_points = {
        libMesh::Point(0.5*h_z, 0.5*h_y) + offset - ps,
        libMesh::Point(0.5*h_z, -0.5*h_y) + offset - ps,
        libMesh::Point(-0.5*h_z, -0.5*h_y) + offset - ps,
        libMesh::Point(-0.5*h_z, 0.5*h_y) + offset - ps
    };
    
    return stress_points;
}


const std::vector<libMesh::Point> 
MAST::Solid1DSectionElementPropertyCard::
get_stress_points_derivative(const MAST::FunctionBase& f, 
                             const libMesh::Point& p, const Real t,
                             const libMesh::Point dps) const
{
    const MAST::FieldFunction<Real>
    &hy     =  this->get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  this->get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  this->get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  this->get<MAST::FieldFunction<Real> >("hz_off");
    
    Real h_y, h_z, h_y_off, h_z_off;
    Real dh_y, dh_z, dh_y_off, dh_z_off;
    
    hy(p, t, h_y);          hy.derivative(f, p, t, dh_y);
    hz(p, t, h_z);          hz.derivative(f, p, t, dh_z);
    hy_off(p, t, h_y_off);  hy_off.derivative(f, p, t, dh_y_off);
    hz_off(p, t, h_z_off);  hz_off.derivative(f, p, t, dh_z_off);
    
    libMesh::Point offset(h_z_off, h_y_off);
    libMesh::Point doffset(dh_z_off, dh_y_off);
                
    std::vector<libMesh::Point> stress_points = {
        libMesh::Point(0.5*dh_z, 0.5*dh_y) + doffset - dps,
        libMesh::Point(0.5*dh_z, -0.5*dh_y) + doffset - dps,
        libMesh::Point(-0.5*dh_z, -0.5*dh_y) + doffset - dps,
        libMesh::Point(-0.5*dh_z, 0.5*dh_y) + doffset - dps
    };
    
    return stress_points;
}
