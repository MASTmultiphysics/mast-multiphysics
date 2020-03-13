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
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"
#include "base/elem_base.h"


namespace MAST {
    namespace Solid2DSectionProperty {
        
        class ExtensionStiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                     const MAST::FieldFunction<Real>& h);
            
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
            const MAST::FieldFunction<Real>& _h;
        };
        
        
        
        class ExtensionBendingStiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionBendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                            const MAST::FieldFunction<Real>& h,
                                            const MAST::FieldFunction<Real>& off);
            
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
            const MAST::FieldFunction<Real>& _h, &_off;
        };
        
        
        class BendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            BendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                   const MAST::FieldFunction<Real>& h,
                                   const MAST::FieldFunction<Real>& off);
            
            
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
            const MAST::FieldFunction<Real>& _h, &_off;
        };
        
        
        
        
        class TransverseStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                      const MAST::FieldFunction<Real>& h,
                                      const MAST::FieldFunction<Real>& kappa):
            MAST::FieldFunction<RealMatrixX>("TransverseStiffnessMatrix2D"),
            _material_stiffness(mat),
            _h(h), _kappa(kappa) {
                _functions.insert(&mat);
                _functions.insert(&h);
                _functions.insert(&kappa);
            }
            
            
            virtual ~TransverseStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real h, kappa;
                _h(p, t, h);
                _kappa(p, t, kappa);
                _material_stiffness(p, t, m);
                //m = m*h*kappa;
                m *= h;
                m *= kappa;
            }
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                RealMatrixX dm;
                Real h, dh, kappa, dkappa;
                _h(p, t, h);        _h.derivative(f, p, t, dh);
                _kappa(p, t, kappa);    _kappa.derivative(f, p, t, dkappa); 
                _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
                
                // dm = dm*h*kappa + m*dh*kappa + m*h*dkappa
                m *= (dh*kappa + h*dkappa);
                m += dm*h*kappa;
            }
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
            const MAST::FieldFunction<Real>& _h;
            const MAST::FieldFunction<Real>& _kappa;
        };
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(const MAST::FieldFunction<Real>& rho,
                          const MAST::FieldFunction<Real>& h,
                          const MAST::FieldFunction<Real>& off);
            
            
            virtual ~InertiaMatrix() { }
            
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _rho, &_h, &_off;
        };
        
        
        
        class ThermalExpansionAMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionAMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                                    const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                                    const MAST::FieldFunction<Real>& h);
            
            
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
            const MAST::FieldFunction<Real>& _h;
        };
        
        
        
        class ThermalExpansionBMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionBMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                                    const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                                    const MAST::FieldFunction<Real>& h,
                                    const MAST::FieldFunction<Real>& off);
            
            
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
            const MAST::FieldFunction<Real>& _h, &_off;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                             const MAST::FieldFunction<RealMatrixX>& T,
                             const MAST::FieldFunction<Real>& h);
            
            
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
            const MAST::FieldFunction<Real>& _h;
        };
        
        
        
        class PrestressBMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressBMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                             const MAST::FieldFunction<RealMatrixX>& T,
                             const MAST::FieldFunction<Real>& h,
                             const MAST::FieldFunction<Real>& off);
            
            
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
            const MAST::FieldFunction<Real>& _h, &_off;
        };
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                                     const MAST::FieldFunction<Real>& h);
            
            
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
            
            const MAST::FieldFunction<Real>& _h;
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
            
            const MAST::FieldFunction<Real>& _h;
        };
    }
}




bool
MAST::Solid2DSectionElementPropertyCard::depends_on(const MAST::FunctionBase& f) const {
    
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}




MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix::
ExtensionStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                         const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX> ("ExtensionStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(&mat);
    _functions.insert(&h);
}



void
MAST::Solid2DSectionProperty::
ExtensionStiffnessMatrix::operator() (const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    // [C]*h
    Real h;
    _h(p, t, h);
    _material_stiffness(p, t, m);
    m *= h;
}




void
MAST::Solid2DSectionProperty::
ExtensionStiffnessMatrix::derivative (     const MAST::FunctionBase& f,
                                      const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    RealMatrixX dm;
    Real h, dhdf;
    _h(p, t, h); _h.derivative( f, p, t, dhdf);
    _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
    
    // [C]*dh
    m *= dhdf;
    
    // += [dC]*h
    m += h*dm;
}






MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix::
ExtensionBendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                                const MAST::FieldFunction<Real>& h,
                                const MAST::FieldFunction<Real>& off):
MAST::FieldFunction<RealMatrixX> ("ExtensionBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(&mat);
    _functions.insert(&h);
    _functions.insert(&off);
}



void
MAST::Solid2DSectionProperty::
ExtensionBendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    // [C]*h
    Real h, off;
    _h(p, t, h);
    _off(p, t, off);
    _material_stiffness(p, t, m);
    m *= h*off;
}




void
MAST::Solid2DSectionProperty::
ExtensionBendingStiffnessMatrix::derivative (            const MAST::FunctionBase& f,
                                             const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(3,3); dm = RealMatrixX::Zero(3, 3);
    Real h, off, dh, doff;
    
    _h(p, t, h); _h.derivative( f, p, t, dh);
    _off(p, t, off); _off.derivative( f, p, t, doff);
    _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
    m *= dh*off + h*doff;
    m += h*off*dm;
}




MAST::Solid2DSectionProperty::BendingStiffnessMatrix::
BendingStiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                       const MAST::FieldFunction<Real>& h,
                       const MAST::FieldFunction<Real>& off):
MAST::FieldFunction<RealMatrixX> ("BendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(&mat);
    _functions.insert(&h);
    _functions.insert(&off);
}



void
MAST::Solid2DSectionProperty::
BendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    // [C]*h
    Real h, off;
    _h(p, t, h);
    _off(p, t, off);
    _material_stiffness(p, t, m);
    m *= (pow(h,3)/12. + h*pow(off,2));
}




void
MAST::Solid2DSectionProperty::
BendingStiffnessMatrix::derivative (   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(3,3); dm = RealMatrixX::Zero(3, 3);
    Real h, dhdf, off, doff;
    _h(p, t, h); _h.derivative( f, p, t, dhdf);
    _off(p, t, off); _off.derivative( f, p, t, doff);
    _material_stiffness(p, t, m); _material_stiffness.derivative( f, p, t, dm);
    
    // [C]*dh
    m *= (pow(h,2)/4.*dhdf + dhdf*pow(off,2) + h*2.*off*doff);
    
    // += [dC]*h
    m += (pow(h,3)/12. + h*pow(off, 2))* dm;
}





MAST::Solid2DSectionProperty::InertiaMatrix::
InertiaMatrix(const MAST::FieldFunction<Real>& rho,
              const MAST::FieldFunction<Real>& h,
              const MAST::FieldFunction<Real>& off):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix2D"),
_rho(rho),
_h(h),
_off(off) {
    _functions.insert(&rho);
    _functions.insert(&h);
    _functions.insert(&off);
}




void
MAST::Solid2DSectionProperty::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m = RealMatrixX::Zero(6, 6);
    Real h, rho, off;
    _h(p, t, h);
    _off(p, t, off);
    _rho(p, t, rho);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = h;
    
    m(0,4) = off*h; m(4,0) = m(0,4);       // extension-bending coupling
    m(1,3) = -off*h; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(4,4) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(5,5) = pow(h,3)/12.*1.0e-6; // neglect the rotary inertia wrt theta_z
    // FIXME: The line above, a small value in m(5,5) can cause many artificial eigenvalues around the true one. Resulting in difficult to interpret modal
    m *= rho;
}





void
MAST::Solid2DSectionProperty::
InertiaMatrix::derivative (               const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m = RealMatrixX::Zero(6,6);
    Real h, dhdf, rho, drhodf, off, doff;
    _h(p, t, h); _h.derivative( f, p, t, dhdf);
    _off(p, t, off); _off.derivative( f, p, t, doff);
    _rho(p, t, rho); _rho.derivative( f, p, t, drhodf);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = drhodf*h + rho*dhdf;
    
    m(0,4) = doff*h+off*dhdf; m(4,0) = m(0,4);        // extension-bending coupling
    m(1,3) = -doff*h-off*dhdf; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(4,4) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(5,5) = (drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf)*1.0e-6; // neglect the rotary inertia wrt theta_z
    // FIXME: The line above, a small value in m(5,5) can cause many artificial eigenvalues around the true one. Resulting in difficult to interpret modal
}




MAST::Solid2DSectionProperty::ThermalExpansionAMatrix::
ThermalExpansionAMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                        const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                        const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionAMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h) {
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
    _functions.insert(&h);
}




void
MAST::Solid2DSectionProperty::
ThermalExpansionAMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX at;
    Real h;
    _h(p, t, h);
    _material_stiffness(p, t, m);
    _material_expansion(p, t, at);
    
    m *= at;
    m *= h;
}






void
MAST::Solid2DSectionProperty::
ThermalExpansionAMatrix::derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX m1, at, dm, dat;
    Real h, dh;
    _h(p, t, h); _h.derivative( f, p, t, dh);
    _material_stiffness(p, t, m1); _material_stiffness.derivative( f, p, t, dm);
    _material_expansion(p, t, at); _material_expansion.derivative( f, p, t, dat);
    
    m = m1;
    
    m *= at;
    m *= dh;
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m += h*m1;
}




MAST::Solid2DSectionProperty::ThermalExpansionBMatrix::
ThermalExpansionBMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                        const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                        const MAST::FieldFunction<Real>& h,
                        const MAST::FieldFunction<Real>& off):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionBMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h),
_off(off) {
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
    _functions.insert(&h);
    _functions.insert(&off);
}




void
MAST::Solid2DSectionProperty::
ThermalExpansionBMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX at;
    Real h, off;
    _h(p, t, h);
    _off(p, t, off);
    _material_stiffness(p, t, m);
    _material_expansion(p, t, at);
    
    m *= at;
    m *= h*off;
}





void
MAST::Solid2DSectionProperty::
ThermalExpansionBMatrix::derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX m1, at, dm, dat;
    Real h, dh, off, doff;
    _h(p, t, h); _h.derivative( f, p, t, dh);
    _off(p, t, off); _off.derivative( f, p, t, doff);
    _material_stiffness(p, t, m1); _material_stiffness.derivative( f, p, t, dm);
    _material_expansion(p, t, at); _material_expansion.derivative( f, p, t, dat);
    
    m = m1;
    
    m *= at;
    m *= (dh*off+h*doff);
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m += h*off*m1;
}




MAST::Solid2DSectionProperty::PrestressAMatrix::
PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                 const MAST::FieldFunction<RealMatrixX>& T,
                 const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix2D"),
_prestress(prestress),
_T(T),
_h(h) {
    _functions.insert(&prestress);
    _functions.insert(&T);
    _functions.insert(&h);
}




void
MAST::Solid2DSectionProperty::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real h;
    _h(p, t, h);
    _prestress(p, t, s);
    _T(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*h;
}






void
MAST::Solid2DSectionProperty::
PrestressAMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real h, dh;
    _h(p, t, h); _h.derivative( f, p, t, dh);
    _prestress(p, t, s); _prestress.derivative( f, p, t, ds);
    _T(p, t, T); _T.derivative( f, p, t, dT);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    RealMatrixX tmp;
    ds *= T;
    ds = T.transpose()*ds;
    
    tmp = s;
    tmp *= dT;
    tmp = T.transpose() * tmp;
    ds += tmp;
    
    tmp = s;
    tmp *= T;
    tmp = dT.transpose() * tmp;
    ds += tmp;
    
    
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*h + s(i,j)*dh;
}




//void
//MAST::Solid2DSectionProperty::
//PrestressAMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v) const {
//    libmesh_assert_equal_to(m.rows(), 2);
//    libmesh_assert_equal_to(m.cols(), 2);
//    v.resize(3);
//    v(0) = m(0,0);  // sigma x
//    v(2) = m(1,1);  // sigma y
//    v(2) = m(0,1);  // tau xy
//}




MAST::Solid2DSectionProperty::PrestressBMatrix::
PrestressBMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                 const MAST::FieldFunction<RealMatrixX>& T,
                 const MAST::FieldFunction<Real>& h,
                 const MAST::FieldFunction<Real>& off):
MAST::FieldFunction<RealMatrixX>("PrestressBMatrix2D"),
_prestress(prestress),
_T(T),
_h(h),
_off(off) {
    _functions.insert(&prestress);
    _functions.insert(&T);
    _functions.insert(&h);
    _functions.insert(&off);
}




void
MAST::Solid2DSectionProperty::
PrestressBMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m = RealMatrixX::Zero(2, 2);
    Real h, off;
    _h(p, t, h);
    _off(p, t, off);
    _prestress(p, t, s);
    _T(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*(h*off);
}





void
MAST::Solid2DSectionProperty::
PrestressBMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m = RealMatrixX::Zero(2, 2);
    Real h, dh, off, doff;
    _h(p, t, h); _h.derivative( f, p, t, dh);
    _off(p, t, off); _off.derivative( f, p, t, doff);
    _prestress(p, t, s); _prestress.derivative( f, p, t, ds);
    _T(p, t, T); _T.derivative( f, p, t, dT);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    RealMatrixX tmp;
    ds *= T;
    ds = T.transpose() * ds;
    
    tmp = s;
    tmp *= dT;
    tmp = dT.transpose() * tmp;
    ds += tmp;
    
    tmp = s;
    tmp *= T;
    tmp = T.transpose()*tmp;
    ds += tmp;
    
    
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*(h*off) + s(i,j)*(dh*off+h*doff);
}




MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                         const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_h(h) {
    _functions.insert(&mat_cond);
    _functions.insert(&h);
}



MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
~ThermalConductanceMatrix() { }


void
MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m = RealMatrixX::Zero(2, 2);
    Real h;
    _mat_cond(p, t, m);
    _h(p, t, h);
    
    m *= h;
}




void
MAST::Solid2DSectionProperty::ThermalConductanceMatrix::derivative (                                   const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m = RealMatrixX::Zero(2, 2);
    RealMatrixX dm;
    Real h, dh;
    _mat_cond(p, t, m);
    _mat_cond.derivative( f, p, t, dm);
    _h(p, t, h);
    _h.derivative( f, p, t, dh);
    
    m *= dh;
    m += dm*h;
}





MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cap,
                         const MAST::FieldFunction<Real>& h):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap),
_h(h) {
    _functions.insert(&mat_cap);
    _functions.insert(&h);
}



MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() { }


void
MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
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
MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
derivative (const MAST::FunctionBase& f,
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



//void
//MAST::Solid2DSectionProperty::
//PrestressBMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v) const {
//    // nothing to be done for a symmetric section
//    v.resize(3);
//}




std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<const FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_A_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<const FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_B_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_D_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_error();
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (nullptr);
}



std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
inertia_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho"),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
inertia_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho"),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     _material->thermal_expansion_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_A_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     _material->thermal_expansion_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     _material->thermal_expansion_matrix(2),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_B_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(2, _if_plane_stress),
     _material->thermal_expansion_matrix(2),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix(),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("kappa")
    );
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
transverse_shear_stiffness_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix(),
     this->get<FieldFunction<Real> >("h"),
     this->get<FieldFunction<Real> >("kappa")
    );
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
prestress_A_matrix( MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid2DSectionProperty::PrestressAMatrix
    // (this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
prestress_B_matrix( MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval;
    // TODO: figure out the interface for prestress and T matrix
    libmesh_assert(false);
    // = new MAST::Solid2DSectionProperty::PrestressBMatrix
    // (this->get<MAST::FieldFunction<RealMatrixX> >("prestress"),
    // e.local_elem().T_matrix(),
    // this->get<FieldFunction<Real> >("h"),
    // this->get<FieldFunction<Real> >("off"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_conductance_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_capacitance_matrix() const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(2),
     this->get<FieldFunction<Real> >("h"));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}

const MAST::FieldFunction<Real>*
MAST::Solid2DSectionElementPropertyCard::
section(const MAST::ElementBase& e) const {
    
    return &(this->get<FieldFunction<Real>>("h"));
}
