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
#include "property_cards/isotropic_material_property_card.h"
#include "base/field_function_base.h"


namespace MAST {
    namespace IsotropicMaterialProperty {
        
        
        class StiffnessMatrix1D:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            StiffnessMatrix1D( MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);
            
            StiffnessMatrix1D(const MAST::IsotropicMaterialProperty::StiffnessMatrix1D& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()){
                _functions.insert(_E->master());
                _functions.insert(_nu->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
            clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::StiffnessMatrix1D(*this));
            }
            
            virtual ~StiffnessMatrix1D() {
                delete _E;
                delete _nu;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseShearStiffnessMatrix( MAST::FieldFunction<Real>* E,
                                           MAST::FieldFunction<Real>* nu,
                                           MAST::FieldFunction<Real>* kappa);
            
            TransverseShearStiffnessMatrix(const MAST::IsotropicMaterialProperty::TransverseShearStiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _kappa(f._kappa->clone().release()) {
                _functions.insert(_E->master());
                _functions.insert(_nu->master());
                _functions.insert(_kappa->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const  {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::TransverseShearStiffnessMatrix(*this));
            }
            
            virtual ~TransverseShearStiffnessMatrix() {
                delete _E;
                delete _nu;
                delete _kappa;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            MAST::FieldFunction<Real>* _kappa;
        };
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix2D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu,
                              bool plane_stress);
            
            StiffnessMatrix2D(const MAST::IsotropicMaterialProperty::StiffnessMatrix2D& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()),
            _plane_stress(f._plane_stress) {
                _functions.insert(_E->master());
                _functions.insert(_nu->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::StiffnessMatrix2D(*this));
            }
            
            virtual ~StiffnessMatrix2D() {
                delete _E;
                delete _nu;
            }
            
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
            bool _plane_stress;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix3D(MAST::FieldFunction<Real>* E,
                              MAST::FieldFunction<Real>* nu);
            
            StiffnessMatrix3D(const MAST::IsotropicMaterialProperty::StiffnessMatrix3D &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _E(f._E->clone().release()),
            _nu(f._nu->clone().release()) {
                _functions.insert(_E->master());
                _functions.insert(_nu->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::StiffnessMatrix3D(*this));
            }
            
            virtual ~StiffnessMatrix3D() {
                delete _E;
                delete _nu;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real>* _E;
            MAST::FieldFunction<Real>* _nu;
        };
        
        
        
        class InertiaMatrix3D: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix3D(MAST::FieldFunction<Real>* rho);
            
            InertiaMatrix3D(const MAST::IsotropicMaterialProperty::InertiaMatrix3D &f):
            MAST::FieldFunction<RealMatrixX>(f),
            _rho(f._rho->clone().release()) {
                _functions.insert(_rho->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::InertiaMatrix3D(*this));
            }
            
            virtual ~InertiaMatrix3D() {
                delete _rho;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<Real> *_rho;
            
        };

        
        
        class ThermalExpansionMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalExpansionMatrix(unsigned int dim,
                                   MAST::FieldFunction<Real>* alpha):
            MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix"),
            _dim(dim),
            _alpha(alpha) {
                _functions.insert(_alpha->master());
            }
            
            ThermalExpansionMatrix(const MAST::IsotropicMaterialProperty::ThermalExpansionMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _dim(f._dim),
            _alpha(f._alpha->clone().release()) {
                _functions.insert(_alpha->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::ThermalExpansionMatrix(*this));
            }
            
            virtual ~ThermalExpansionMatrix() {
                delete _alpha;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            MAST::FieldFunction<Real>* _alpha;
        };
        
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalConductanceMatrix(unsigned int dim,
                              MAST::FieldFunction<Real>* k):
            MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
            _dim(dim),
            _k(k) {
                _functions.insert(_k->master());
            }
            
            ThermalConductanceMatrix(const MAST::IsotropicMaterialProperty::ThermalConductanceMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _dim(f._dim),
            _k(f._k->clone().release()) {
                _functions.insert(_k->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::ThermalConductanceMatrix(*this));
            }
            
            virtual ~ThermalConductanceMatrix() {
                delete _k;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            MAST::FieldFunction<Real>* _k;
        };

        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalCapacitanceMatrix(unsigned int dim,
                                     MAST::FieldFunction<Real>* rho,
                                     MAST::FieldFunction<Real>* cp):
            MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
            _dim(dim),
            _rho(rho),
            _cp(cp) {
                
                _functions.insert(_rho->master());
                _functions.insert(_cp->master());
            }
            
            ThermalCapacitanceMatrix(const MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _dim(f._dim),
            _rho(f._rho->clone().release()),
            _cp(f._cp->clone().release()) {
                
                _functions.insert(_rho->master());
                _functions.insert(_cp->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix(*this));
            }
            
            virtual ~ThermalCapacitanceMatrix() {
                
                delete _rho;
                delete _cp;
            }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            MAST::FieldFunction<Real>* _rho;

            MAST::FieldFunction<Real>* _cp;
        };

        
    }
}



MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::StiffnessMatrix1D(MAST::FieldFunction<Real>* E,
                                     MAST::FieldFunction<Real>* nu ):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix1D"),
_E(E),
_nu(nu)
{
    _functions.insert(E->master());
    _functions.insert(nu->master());
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m  = RealMatrixX::Zero(2,2);
    Real E, nu, G;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    G = E/2./(1.+nu);
    m(0,0) = E;
    m(1,1) = G;
}


void
MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::derivative(const MAST::DerivativeType d,
                              const MAST::FunctionBase &f,
                              const libMesh::Point &p,
                              const Real t,
                              RealMatrixX &m) const {
    
    
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2); dm = RealMatrixX::Zero(2,2);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E);     _E->derivative(d, f, p, t, dEdf);
    (*_nu)(p, t, nu);  _nu->derivative(d, f, p, t, dnudf);
    
    // parM/parE * parE/parf
    dm(0,0) = 1.;
    dm(1,1) = 1./2./(1.+nu);
    m += dEdf * dm;
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = 0.;
    dm(1,1) = -E/2./pow(1.+nu,2);
    m+= dnudf*dm;
}


MAST::IsotropicMaterialProperty::
TransverseShearStiffnessMatrix::TransverseShearStiffnessMatrix(MAST::FieldFunction<Real> * E,
                                                               MAST::FieldFunction<Real> * nu,
                                                               MAST::FieldFunction<Real> * kappa):
MAST::FieldFunction<RealMatrixX>("TransverseShearStiffnessMatrix"),
_E(E),
_nu(nu),
_kappa(kappa)
{
    _functions.insert(E->master());
    _functions.insert(nu->master());
    _functions.insert(kappa->master());
}



void
MAST::IsotropicMaterialProperty::
TransverseShearStiffnessMatrix::operator() (const libMesh::Point& p,
                                            const Real t,
                                            RealMatrixX& m) const {
    m = RealMatrixX::Zero(2,2);
    Real E, nu, kappa, G;
    (*_E)(p, t, E); (*_nu)(p, t, nu); (*_kappa)(p, t, kappa);
    G = E/2./(1.+nu);
    m(0,0) = G*kappa;
    m(1,1) = m(0,0);
}



void
MAST::IsotropicMaterialProperty::
TransverseShearStiffnessMatrix::derivative(const MAST::DerivativeType d,
                                           const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2); dm = RealMatrixX::Zero(2, 2);
    Real E, nu, kappa, dEdf, dnudf, dkappadf, G;
    (*_E)    (p, t, E);         _E->derivative(d, f, p, t, dEdf);
    (*_nu)   (p, t, nu);       _nu->derivative(d, f, p, t, dnudf);
    (*_kappa)(p, t, kappa); _kappa->derivative(d, f, p, t, dkappadf);
    G = E/2./(1.+nu);
    
    
    // parM/parE * parE/parf
    dm(0,0) = 1./2./(1.+nu)*kappa;
    dm(1,1) = dm(0,0);
    m += dEdf * dm;
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = -E/2./pow(1.+nu,2)*kappa;
    dm(1,1) = dm(0,0);
    m += dnudf*dm;
    
    // parM/parnu * parkappa/parf
    
    dm(0,0) = G; dm(1,1) = G;
    m += dkappadf*dm;
}




MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::StiffnessMatrix2D(MAST::FieldFunction<Real> * E,
                                     MAST::FieldFunction<Real> * nu ,
                                     bool plane_stress ):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix2D"),
_E(E),
_nu(nu),
_plane_stress(plane_stress)
{
    _functions.insert(E->master());
    _functions.insert(nu->master());
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    libmesh_assert(_plane_stress); // currently only implemented for plane stress
    m = RealMatrixX::Zero(3,3);
    Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                m(i,i) = E/(1.-nu*nu);
            else // offdiagonal: direct stress
                m(i,j) = E*nu/(1.-nu*nu);
    }
    m(2,2) = E/2./(1.+nu); // diagonal: shear stress
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::derivative (const MAST::DerivativeType d,
                               const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    libmesh_assert(_plane_stress); // currently only implemented for plane stress
    RealMatrixX dm;
    m = RealMatrixX::Zero(3,3); dm = RealMatrixX::Zero(3, 3);
    Real E, nu, dEdf, dnudf;
    (*_E)  (p, t, E);   _E->derivative(d, f, p, t, dEdf);
    (*_nu) (p, t, nu); _nu->derivative(d, f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = 1./(1.-nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = 1.*nu/(1.-nu*nu);
    }
    dm(2,2) = 1./2./(1.+nu); // diagonal: shear stress
    m += dEdf * dm;
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = E/pow(1.-nu*nu, 2)*2.*nu;
            else // offdiagonal: direct stress
                dm(i,j) = E/(1.-nu*nu) + E*nu/pow(1.-nu*nu,2)*2.*nu;
    }
    dm(2,2) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    m+= dnudf*dm;
}



MAST::IsotropicMaterialProperty::
StiffnessMatrix3D::StiffnessMatrix3D(MAST::FieldFunction<Real> * E,
                                     MAST::FieldFunction<Real> * nu):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix3D"),
_E(E),
_nu(nu)
{
    _functions.insert(E->master());
    _functions.insert(nu->master());
}





void
MAST::IsotropicMaterialProperty::
StiffnessMatrix3D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m = RealMatrixX::Zero(6,6);
    Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                m(i,i) = E*(1.-nu)/(1.-nu-2.*nu*nu);
            else // offdiagonal: direct stress
                m(i,j) = E*nu/(1.-nu-2.*nu*nu);
        m(i+3,i+3) = E/2./(1.+nu); // diagonal: shear stress
    }
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix3D::derivative (const MAST::DerivativeType d,
                               const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(6,6); dm = RealMatrixX::Zero(6,6);
    Real E, nu, dEdf, dnudf;
    (*_E)  (p, t, E);   _E->derivative(d, f, p, t, dEdf);
    (*_nu) (p, t, nu); _nu->derivative(d, f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = nu/(1.-nu-2.*nu*nu);
        dm(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
    }
    m += dEdf * dm;
    
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
            else // offdiagonal: direct stress
                dm(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
        dm(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    }
    m+= dnudf*dm;
}



MAST::IsotropicMaterialProperty::
InertiaMatrix3D::InertiaMatrix3D(MAST::FieldFunction<Real> * rho):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_rho(rho) {
    
    _functions.insert(rho->master());
}


void
MAST::IsotropicMaterialProperty::
InertiaMatrix3D::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    m  = RealMatrixX::Zero(3,3);
    Real rho;
    (*_rho)(p, t, rho);
    m(0,0) = rho;
    m(1,1) = rho;
    m(2,2) = rho;
}


void
MAST::IsotropicMaterialProperty::
InertiaMatrix3D::derivative(const MAST::DerivativeType d,
                            const MAST::FunctionBase &f,
                            const libMesh::Point &p,
                            const Real t,
                            RealMatrixX &m) const {
    
    
    m  = RealMatrixX::Zero(3,3);
    Real rho;
    _rho->derivative(d, f, p, t, rho);
    m(0,0) = rho;
    m(1,1) = rho;
    m(2,2) = rho;
}



void
MAST::IsotropicMaterialProperty::ThermalExpansionMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real alpha;
    (*_alpha)(p, t, alpha);
    switch (_dim) {
        case 1:
            m = RealMatrixX::Zero(2,1);
            break;
            
        case 2:
            m = RealMatrixX::Zero(3,1);
            break;
            
        case 3:
            m = RealMatrixX::Zero(6,1);
            break;
    }
    
    for (unsigned int i=0; i<_dim; i++)
        m(i,0) = alpha;
}





void
MAST::IsotropicMaterialProperty::ThermalExpansionMatrix::
derivative (const MAST::DerivativeType d,
            const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real alpha;
    _alpha->derivative(d, f, p, t, alpha);
    switch (_dim) {
        case 1:
            m = RealMatrixX::Zero(2,1);
            break;
            
        case 2:
            m = RealMatrixX::Zero(3,1);
            break;
            
        case 3:
            m = RealMatrixX::Zero(6,1);
            break;
    }
    
    for (unsigned int i=0; i<_dim; i++)
        m(i,0) = alpha;
    
}





void
MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real cp, rho;
    (*_cp)  (p, t, cp);
    (*_rho) (p, t, rho);
    
    m.setZero(1,1);
    
    m(0,0) = cp*rho;
}





void
MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix::
derivative (const MAST::DerivativeType d,
            const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real cp, dcp, rho, drho;
    (*_cp)  (p, t, cp);    _cp->derivative(d, f, p, t, dcp);
    (*_rho) (p, t, rho);  _rho->derivative(d, f, p, t, drho);
    
    m.setZero(1,1);
    
    m(0,0) = dcp*rho + cp*drho;
}




void
MAST::IsotropicMaterialProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real k;
    (*_k)  (p, t, k);
    
    m.setIdentity(_dim, _dim);
    
    m *= k;
}





void
MAST::IsotropicMaterialProperty::ThermalConductanceMatrix::
derivative (const MAST::DerivativeType d,
            const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real k;
    _k->derivative(d, f, p, t, k);
    
    m.setIdentity(_dim, _dim);
    
    m *= k;
}





std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::stiffness_matrix(const unsigned int dim,
                                                      const bool plane_stress) {
    
    MAST::FieldFunction<RealMatrixX> *rval = NULL;
    
    switch (dim) {
        case 1:
            rval = new MAST::IsotropicMaterialProperty::StiffnessMatrix1D
            (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
             this->get<MAST::FieldFunction<Real> >("nu").clone().release());
            break;
            
        case 2:
            rval = new MAST::IsotropicMaterialProperty::StiffnessMatrix2D
            (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
             this->get<MAST::FieldFunction<Real> >("nu").clone().release(),
             plane_stress);
            break;
            
        case 3:
            rval = new MAST::IsotropicMaterialProperty::StiffnessMatrix3D
            (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
             this->get<MAST::FieldFunction<Real> >("nu").clone().release());
            break;
            
        default:
            // should not get here
            libmesh_error();
    }
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::damping_matrix(const unsigned int dim) {
    
    MAST::FieldFunction<RealMatrixX> *rval = NULL;
    
    // make sure that this is not null
    libmesh_assert(rval);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::inertia_matrix(const unsigned int dim) {
    
    MAST::FieldFunction<RealMatrixX> *rval = NULL;
    
    switch (dim) {
        case 3:
            rval = new MAST::IsotropicMaterialProperty::InertiaMatrix3D
            (this->get<MAST::FieldFunction<Real> >("rho").clone().release());
            break;
            
        case 1:
        case 2:
        default:
            // implemented only for 3D since the 2D and 1D elements
            // calculate it themselves
            libmesh_error();
            
    }
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::thermal_expansion_matrix(const unsigned int dim) {
    
    MAST::FieldFunction<RealMatrixX> *rval =
    new MAST::IsotropicMaterialProperty::ThermalExpansionMatrix
    (dim,
     this->get<MAST::FieldFunction<Real> >("alpha_expansion").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::transverse_shear_stiffness_matrix() {
    
    MAST::FieldFunction<RealMatrixX> *rval =
    new MAST::IsotropicMaterialProperty::TransverseShearStiffnessMatrix
    (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
     this->get<MAST::FieldFunction<Real> >("nu").clone().release(),
     this->get<MAST::FieldFunction<Real> >("kappa").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::capacitance_matrix(const unsigned int dim) {
    
    MAST::FieldFunction<RealMatrixX> *rval =
    new MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix
    (dim,
     this->get<MAST::FieldFunction<Real> >("rho").clone().release(),
     this->get<MAST::FieldFunction<Real> >("cp").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}





std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicMaterialPropertyCard::conductance_matrix(const unsigned int dim) {
    
    MAST::FieldFunction<RealMatrixX> *rval =
    new MAST::IsotropicMaterialProperty::ThermalConductanceMatrix
    (dim,
     this->get<MAST::FieldFunction<Real> >("k_th").clone().release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




