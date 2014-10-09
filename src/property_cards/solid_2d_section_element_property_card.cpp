/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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
#include "mesh/local_elem_base.h"


namespace MAST {
    namespace Solid2DSectionProperty {
        
        class ExtensionStiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                     MAST::FieldFunction<Real> *h);
            
            ExtensionStiffnessMatrix(const MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
            }
            
            
            virtual ~ExtensionStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix(*this));
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
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class ExtensionBendingStiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            ExtensionBendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                            MAST::FieldFunction<Real> *h,
                                            MAST::FieldFunction<Real> *off);
            
            ExtensionBendingStiffnessMatrix(const MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            virtual ~ExtensionBendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix(*this));
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
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        class BendingStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            BendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                   MAST::FieldFunction<Real> *h,
                                   MAST::FieldFunction<Real> *off);
            
            BendingStiffnessMatrix(const MAST::Solid2DSectionProperty::BendingStiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            
            virtual ~BendingStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::BendingStiffnessMatrix(*this));
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
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        
        
        class TransverseStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                      MAST::FieldFunction<Real>* h):
            MAST::FieldFunction<RealMatrixX>("TransverseStiffnessMatrix2D"),
            _material_stiffness(mat),
            _h(h) {
                _functions.insert(mat);
                _functions.insert(h);
            }
            
            
            TransverseStiffnessMatrix(const MAST::Solid2DSectionProperty::TransverseStiffnessMatrix &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_h);
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::TransverseStiffnessMatrix(*this));
            }
            
            virtual ~TransverseStiffnessMatrix() {
                delete _material_stiffness;
                delete _h;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                Real h;
                (*_h)(p, t, h);
                (*_material_stiffness)(p, t, m);
                m *= h;
            }
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
                RealMatrixX dm;
                Real h, dh;
                (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
                (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
                
                m *= dh;
                m += h*dm;
            }
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_stiffness;
            MAST::FieldFunction<Real> *_h;
        };
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(MAST::FieldFunction<Real> *rho,
                          MAST::FieldFunction<Real> *h,
                          MAST::FieldFunction<Real> *off);
            
            InertiaMatrix(const MAST::Solid2DSectionProperty::InertiaMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _rho(f._rho->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_rho);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            
            virtual ~InertiaMatrix() {
                delete _rho;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::InertiaMatrix(*this));
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
            
            MAST::FieldFunction<Real> *_rho, *_h, *_off;
        };
        
        
        
        class ThermalExpansionAMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionAMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                                    MAST::FieldFunction<RealMatrixX> *mat_expansion,
                                    MAST::FieldFunction<Real> *h);
            
            ThermalExpansionAMatrix(const MAST::Solid2DSectionProperty::ThermalExpansionAMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
                _functions.insert(_h);
            }
            
            
            
            virtual ~ThermalExpansionAMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
            }
            
            
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::ThermalExpansionAMatrix(*this));
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
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class ThermalExpansionBMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionBMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                                    MAST::FieldFunction<RealMatrixX> *mat_expansion,
                                    MAST::FieldFunction<Real> *h,
                                    MAST::FieldFunction<Real> *off);
            
            ThermalExpansionBMatrix(const MAST::Solid2DSectionProperty::ThermalExpansionBMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            
            virtual ~ThermalExpansionBMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::ThermalExpansionBMatrix(*this));
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
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                             MAST::FieldFunction<RealMatrixX> *T,
                             MAST::FieldFunction<Real> *h);
            
            PrestressAMatrix(const MAST::Solid2DSectionProperty::PrestressAMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release()) {
                _functions.insert(_prestress);
                _functions.insert(_T);
                _functions.insert(_h);
            }
            
            
            virtual ~PrestressAMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::PrestressAMatrix(*this));
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
            MAST::FieldFunction<Real> *_h;
        };
        
        
        
        class PrestressBMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressBMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                             MAST::FieldFunction<RealMatrixX> *T,
                             MAST::FieldFunction<Real> *h,
                             MAST::FieldFunction<Real> *off);
            
            PrestressBMatrix(const MAST::Solid2DSectionProperty::PrestressBMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()),
            _T(f._T->clone().release()),
            _h(f._h->clone().release()),
            _off(f._off->clone().release()) {
                _functions.insert(_prestress);
                _functions.insert(_T);
                _functions.insert(_h);
                _functions.insert(_off);
            }
            
            
            virtual ~PrestressBMatrix() {
                delete _prestress;
                delete _T;
                delete _h;
                delete _off;
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::Solid2DSectionProperty::PrestressBMatrix(*this));
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
            MAST::FieldFunction<Real> *_h, *_off;
        };
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                                     MAST::FieldFunction<Real> *h);
            
            ThermalConductanceMatrix(const MAST::Solid2DSectionProperty::ThermalConductanceMatrix &f);
            
            
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
            
            MAST::FieldFunction<Real>* _h;
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                                     MAST::FieldFunction<Real> *h);
            
            ThermalCapacitanceMatrix(const MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix &f);
            
            
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




bool
MAST::Solid2DSectionElementPropertyCard::depends_on(const MAST::FunctionBase& f) const {
    
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}




MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix::
ExtensionStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                         MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX> ("ExtensionStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat);
    _functions.insert(h);
}



void
MAST::Solid2DSectionProperty::
ExtensionStiffnessMatrix::operator() (const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    // [C]*h
    Real h;
    (*_h)(p, t, h);
    (*_material_stiffness)(p, t, m);
    m *= h;
}




void
MAST::Solid2DSectionProperty::
ExtensionStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                      const MAST::FunctionBase& f,
                                      const libMesh::Point& p,
                                      const Real t,
                                      RealMatrixX& m) const {
    RealMatrixX dm;
    Real h, dhdf;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dhdf);
    (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
    
    // [C]*dh
    m *= dhdf;
    
    // += [dC]*h
    m += h*dm;
}






MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix::
ExtensionBendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                                MAST::FieldFunction<Real> *h,
                                MAST::FieldFunction<Real> *off):
MAST::FieldFunction<RealMatrixX> ("ExtensionBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(mat);
    _functions.insert(h);
    _functions.insert(off);
}



void
MAST::Solid2DSectionProperty::
ExtensionBendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    // [C]*h
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    m *= h*off;
}




void
MAST::Solid2DSectionProperty::
ExtensionBendingStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                             const MAST::FunctionBase& f,
                                             const libMesh::Point& p,
                                             const Real t,
                                             RealMatrixX& m) const {
    RealMatrixX dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, off, dh, doff;
    
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
    (*_off)(p, t, off); _off->derivative(d, f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
    m *= dh*off + h*doff;
    m += h*off*dm;
}




MAST::Solid2DSectionProperty::BendingStiffnessMatrix::
BendingStiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                       MAST::FieldFunction<Real> *h,
                       MAST::FieldFunction<Real> *off):
MAST::FieldFunction<RealMatrixX> ("BendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(mat);
    _functions.insert(h);
    _functions.insert(off);
}



void
MAST::Solid2DSectionProperty::
BendingStiffnessMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    // [C]*h
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    m *= pow(h,3)/12. + h*pow(off,2);
}




void
MAST::Solid2DSectionProperty::
BendingStiffnessMatrix::derivative (const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, dhdf, off, doff;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dhdf);
    (*_off)(p, t, off); _h->derivative(d, f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->derivative(d, f, p, t, dm);
    
    // [C]*dh
    m *= pow(h,2)/4.*dhdf + dhdf*pow(off,2) + h*2.*off*doff;
    
    // += [dC]*h
    m += (pow(h,3)/12. + h*pow(off, 2))* dm;
}





MAST::Solid2DSectionProperty::InertiaMatrix::
InertiaMatrix(MAST::FieldFunction<Real> *rho,
              MAST::FieldFunction<Real> *h,
              MAST::FieldFunction<Real> *off):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix2D"),
_rho(rho),
_h(h),
_off(off) {
    _functions.insert(rho);
    _functions.insert(h);
    _functions.insert(off);
}




void
MAST::Solid2DSectionProperty::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m.resize(6, 6);
    Real h, rho, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_rho)(p, t, rho);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = h;
    
    m(0,4) = off*h; m(4,0) = m(0,4);       // extension-bending coupling
    m(1,3) = -off*h; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(4,4) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(5,5) = pow(h,3)/12.*1.0e-12; // neglect the rotary inertia wrt theta_z
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<2; i++)
        m(i+3,i+3) *= 1.0e-16;
    
    m *= rho;
}





void
MAST::Solid2DSectionProperty::
InertiaMatrix::derivative (const MAST::DerivativeType d,
                           const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    m.resize(6,6);
    Real h, dhdf, rho, drhodf, off, doff;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dhdf);
    (*_off)(p, t, off); _off->derivative(d, f, p, t, doff);
    (*_rho)(p, t, rho); _rho->derivative(d, f, p, t, drhodf);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = drhodf*h + rho*dhdf;
    
    m(0,4) = doff*h+off*dhdf; m(4,0) = m(0,4);        // extension-bending coupling
    m(1,3) = -doff*h-off*dhdf; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(4,4) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(5,5) = (drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf)*1.0e-12; // neglect the rotary inertia wrt theta_z
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<2; i++)
        m(i+3,i+3) *= 1.0e-16;
}




MAST::Solid2DSectionProperty::ThermalExpansionAMatrix::
ThermalExpansionAMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                        MAST::FieldFunction<RealMatrixX> *mat_expansion,
                        MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionAMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
    _functions.insert(h);
}




void
MAST::Solid2DSectionProperty::
ThermalExpansionAMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX at;
    Real h;
    (*_h)(p, t, h);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m *= at;
    m *= h;
}






void
MAST::Solid2DSectionProperty::
ThermalExpansionAMatrix::derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX m1, at, dm, dat;
    Real h, dh;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
    (*_material_stiffness)(p, t, m1); _material_stiffness->derivative(d, f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->derivative(d, f, p, t, dat);
    
    m = m1;
    
    m *= at;
    m *= dh;
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m += h*m1;
}




MAST::Solid2DSectionProperty::ThermalExpansionBMatrix::
ThermalExpansionBMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                        MAST::FieldFunction<RealMatrixX> *mat_expansion,
                        MAST::FieldFunction<Real> *h,
                        MAST::FieldFunction<Real> *off):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionBMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h),
_off(off) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
    _functions.insert(h);
    _functions.insert(off);
}




void
MAST::Solid2DSectionProperty::
ThermalExpansionBMatrix::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX at;
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m *= at;
    m *= h*off;
}





void
MAST::Solid2DSectionProperty::
ThermalExpansionBMatrix::derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const {
    RealMatrixX m1, at, dm, dat;
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
    (*_off)(p, t, off); _off->derivative(d, f, p, t, doff);
    (*_material_stiffness)(p, t, m1); _material_stiffness->derivative(d, f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->derivative(d, f, p, t, dat);
    
    m = m1;
    
    m *= at;
    m *= (dh*off+h*doff);
    
    m1 *= dat;
    dm *= at;
    m1 += dm;
    
    m += h*off*m1;
}




MAST::Solid2DSectionProperty::PrestressAMatrix::
PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                 MAST::FieldFunction<RealMatrixX> *T,
                 MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix2D"),
_prestress(prestress),
_T(T),
_h(h) {
    _functions.insert(prestress);
    _functions.insert(T);
    _functions.insert(h);
}




void
MAST::Solid2DSectionProperty::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m.resize(2, 2);
    Real h;
    (*_h)(p, t, h);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*h;
}






void
MAST::Solid2DSectionProperty::
PrestressAMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
    (*_prestress)(p, t, s); _prestress->derivative(d, f, p, t, ds);
    (*_T)(p, t, T); _T->derivative(d, f, p, t, dT);
    
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
PrestressBMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                 MAST::FieldFunction<RealMatrixX> *T,
                 MAST::FieldFunction<Real> *h,
                 MAST::FieldFunction<Real> *off):
MAST::FieldFunction<RealMatrixX>("PrestressBMatrix2D"),
_prestress(prestress),
_T(T),
_h(h),
_off(off) {
    _functions.insert(prestress);
    _functions.insert(T);
    _functions.insert(h);
    _functions.insert(off);
}




void
MAST::Solid2DSectionProperty::
PrestressBMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, T;
    m.resize(2, 2);
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s *= T;
    s = T.transpose() * s;
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*(h*off);
}





void
MAST::Solid2DSectionProperty::
PrestressBMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    RealMatrixX s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->derivative(d, f, p, t, dh);
    (*_off)(p, t, off); _off->derivative(d, f, p, t, doff);
    (*_prestress)(p, t, s); _prestress->derivative(d, f, p, t, ds);
    (*_T)(p, t, T); _T->derivative(d, f, p, t, dT);
    
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
ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                         MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_h(h) {
    _functions.insert(mat_cond);
    _functions.insert(h);
}


MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::Solid2DSectionProperty::ThermalConductanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cond(f._mat_cond),
_h(f._h) {
    _functions.insert(_mat_cond);
    _functions.insert(_h);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionProperty::ThermalConductanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalConductanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
~ThermalConductanceMatrix() {
    
    delete _mat_cond;
    delete _h;
}


void
MAST::Solid2DSectionProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m.resize(2, 2);
    Real h;
    (*_mat_cond)(p, t, m);
    (*_h)(p, t, h);
    
    m *= h;
}




void
MAST::Solid2DSectionProperty::ThermalConductanceMatrix::derivative (const MAST::DerivativeType d,
                                                                    const MAST::FunctionBase& f,
                                                                    const libMesh::Point& p,
                                                                    const Real t,
                                                                    RealMatrixX& m) const {
    m.resize(2, 2);
    RealMatrixX dm;
    Real h, dh;
    (*_mat_cond)(p, t, m);
    _mat_cond->derivative(d, f, p, t, dm);
    (*_h)(p, t, h);
    _h->derivative(d, f, p, t, dh);
    
    m *= dh;
    m += dm*h;
}





MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cap,
                         MAST::FieldFunction<Real> *h):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap),
_h(h) {
    _functions.insert(mat_cap);
    _functions.insert(h);
}


MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cap(f._mat_cap),
_h(f._h) {
    _functions.insert(_mat_cap);
    _functions.insert(_h);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() {
    
    delete _mat_cap;
    delete _h;
}


void
MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m.resize(1, 1);
    Real h;
    (*_mat_cap)(p, t, m);
    (*_h)(p, t, h);
    
    m *= h;
}




void
MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix::
derivative (const MAST::DerivativeType d,
            const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    m.resize(1, 1);
    RealMatrixX dm;
    Real h, dh;
    (*_mat_cap)(p, t, m);
    _mat_cap->derivative(d, f, p, t, dm);
    (*_h)(p, t, h);
    _h->derivative(d, f, p, t, dh);
    
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




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionStiffnessMatrix
    (_material->stiffness_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ExtensionBendingStiffnessMatrix
    (_material->stiffness_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release(),
     this->get<FieldFunction<Real> >("off").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::BendingStiffnessMatrix
    (_material->stiffness_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release(),
     this->get<FieldFunction<Real> >("off").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_error();
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (NULL);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
inertia_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::InertiaMatrix
    (_material->get<FieldFunction<Real> >("rho").clone().release(),
     this->get<FieldFunction<Real> >("h").clone().release(),
     this->get<FieldFunction<Real> >("off").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionAMatrix
    (_material->stiffness_matrix(2).release(),
     _material->thermal_expansion_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalExpansionBMatrix
    (_material->stiffness_matrix(2).release(),
     _material->thermal_expansion_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release(),
     this->get<FieldFunction<Real> >("off").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::TransverseStiffnessMatrix
    (_material->transverse_shear_stiffness_matrix().release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
prestress_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::PrestressAMatrix
    (this->get<MAST::FieldFunction<RealMatrixX> >("prestress").clone().release(),
     e.local_elem().T_matrix_function().release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
prestress_B_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::PrestressBMatrix
    (this->get<MAST::FieldFunction<RealMatrixX> >("prestress").clone().release(),
     e.local_elem().T_matrix_function().release(),
     this->get<FieldFunction<Real> >("h").clone().release(),
     this->get<FieldFunction<Real> >("off").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalConductanceMatrix
    (_material->conductance_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::Solid2DSectionElementPropertyCard::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::Solid2DSectionProperty::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(2).release(),
     this->get<FieldFunction<Real> >("h").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> > (rval);
}



