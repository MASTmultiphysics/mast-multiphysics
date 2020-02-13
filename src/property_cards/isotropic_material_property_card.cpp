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
#include "property_cards/isotropic_material_property_card.h"
#include "base/field_function_base.h"


namespace MAST {
    namespace IsotropicMaterialProperty {
        
        
        class StiffnessMatrix1D:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            StiffnessMatrix1D(const MAST::FieldFunction<Real>& E,
                              const MAST::FieldFunction<Real>& nu);
            
            virtual ~StiffnessMatrix1D() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _E;
            const MAST::FieldFunction<Real>& _nu;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseShearStiffnessMatrix(const MAST::FieldFunction<Real>& E,
                                           const MAST::FieldFunction<Real>& nu);
            
            
            virtual ~TransverseShearStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _E;
            const MAST::FieldFunction<Real>& _nu;
        };
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix2D(const MAST::FieldFunction<Real>& E,
                              const MAST::FieldFunction<Real>& nu,
                              bool plane_stress);
            
            
            virtual ~StiffnessMatrix2D() { }
            
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _E;
            const MAST::FieldFunction<Real>& _nu;
            bool _plane_stress;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix3D(const MAST::FieldFunction<Real>& E,
                              const MAST::FieldFunction<Real>& nu);
            
            virtual ~StiffnessMatrix3D() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _E;
            const MAST::FieldFunction<Real>& _nu;
        };
        
        
        
        class InertiaMatrix3D: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix3D(const MAST::FieldFunction<Real>& rho);
            

            virtual ~InertiaMatrix3D() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _rho;
            
        };

        
        
        class ThermalExpansionMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalExpansionMatrix(unsigned int dim,
                                   const MAST::FieldFunction<Real>& alpha):
            MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix"),
            _dim(dim),
            _alpha(alpha) {
                _functions.insert(&_alpha);
            }
            

            virtual ~ThermalExpansionMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            const MAST::FieldFunction<Real>& _alpha;
        };
        
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalConductanceMatrix(unsigned int dim,
                              MAST::FieldFunction<Real>& k):
            MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
            _dim(dim),
            _k(k) {
                _functions.insert(&_k);
            }
            
            
            virtual ~ThermalConductanceMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            const MAST::FieldFunction<Real>& _k;
        };

        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalCapacitanceMatrix(unsigned int dim,
                                     const MAST::FieldFunction<Real>& rho,
                                     const MAST::FieldFunction<Real>& cp):
            MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
            _dim(dim),
            _rho(rho),
            _cp(cp) {
                
                _functions.insert(&_rho);
                _functions.insert(&_cp);
            }
            
            
            virtual ~ThermalCapacitanceMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const unsigned int _dim;
            
            const MAST::FieldFunction<Real>& _rho;

            const MAST::FieldFunction<Real>& _cp;
        };

        
    }
}



MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::StiffnessMatrix1D(const MAST::FieldFunction<Real>& E,
                                     const MAST::FieldFunction<Real>& nu ):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix1D"),
_E(E),
_nu(nu)
{
    _functions.insert(&E);
    _functions.insert(&nu);
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m  = RealMatrixX::Zero(2,2);
    Real E, nu, G;
    _E(p, t, E); _nu(p, t, nu);
    G = E/2./(1.+nu);
    m(0,0) = E;
    m(1,1) = G;
}


void
MAST::IsotropicMaterialProperty::
StiffnessMatrix1D::derivative( const MAST::FunctionBase &f,
                              const libMesh::Point &p,
                              const Real t,
                              RealMatrixX &m) const {
    
    
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2); dm = RealMatrixX::Zero(2,2);
    Real E, nu, dEdf, dnudf;
    _E(p, t, E);     _E.derivative( f, p, t, dEdf);
    _nu(p, t, nu);  _nu.derivative( f, p, t, dnudf);
    
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
TransverseShearStiffnessMatrix::
TransverseShearStiffnessMatrix(const MAST::FieldFunction<Real>& E,
                               const MAST::FieldFunction<Real>& nu):
MAST::FieldFunction<RealMatrixX>("TransverseShearStiffnessMatrix"),
_E(E),
_nu(nu)
{
    _functions.insert(&E);
    _functions.insert(&nu);
}



void
MAST::IsotropicMaterialProperty::
TransverseShearStiffnessMatrix::operator() (const libMesh::Point& p,
                                            const Real t,
                                            RealMatrixX& m) const {
    m = RealMatrixX::Zero(2,2);
    Real E, nu, G;
    _E(p, t, E); _nu(p, t, nu);
    G = E/2./(1.+nu);
    m(0,0) = G;
    m(1,1) = m(0,0);
}



void
MAST::IsotropicMaterialProperty::
TransverseShearStiffnessMatrix::derivative(          const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2); dm = RealMatrixX::Zero(2, 2);
    Real E, nu, dEdf, dnudf, G, dG;
    _E    (p, t, E);         _E.derivative( f, p, t, dEdf);
    _nu   (p, t, nu);       _nu.derivative( f, p, t, dnudf);
    G = E/2./(1.+nu);       
    dG = (1./2./(1.+nu) * dEdf) + (-E/2./pow(1.+nu,2) * dnudf);
    
    dm(0,0) = dm(1,1) = dG;
    
//     // parM/parE * parE/parf
//     dm(0,0) = 1./2./(1.+nu);
//     dm(1,1) = dm(0,0);
//     m += dEdf * dm;
//     
//     // parM/parnu * parnu/parf
//     dm(0,0) = -E/2./pow(1.+nu,2);
//     dm(1,1) = dm(0,0);
//     m += dnudf*dm;
}




MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::StiffnessMatrix2D(const MAST::FieldFunction<Real>& E,
                                     const MAST::FieldFunction<Real>& nu ,
                                     bool plane_stress ):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix2D"),
_E(E),
_nu(nu),
_plane_stress(plane_stress) {
    
    _functions.insert(&E);
    _functions.insert(&nu);
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    // libmesh_assert(_plane_stress); // currently only implemented for plane stress
    m = RealMatrixX::Zero(3,3);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
    if (_plane_stress)
    {
        for (unsigned int i=0; i<2; i++) {
            for (unsigned int j=0; j<2; j++)
                if (i == j) // diagonal: direct stress
                    m(i,i) = E/(1.-nu*nu);
                else // offdiagonal: direct stress
                    m(i,j) = E*nu/(1.-nu*nu);
        }
        m(2,2) = E/2./(1.+nu); // diagonal: shear stress
    }
    else // plane_strain
    {
        for (unsigned int i=0; i<2; i++) {
            for (unsigned int j=0; j<2; j++)
                if (i == j) // diagonal: direct stress
                    m(i,i) = E*(1.-nu)/((1.+nu)*(1.-2.*nu));
                else // offdiagonal: direct stress
                    m(i,j) = E*nu/((1.+nu)*(1.-2.*nu));
        }
        m(2,2) = E/2./(1.+nu); // diagonal: shear stress
    }
}




void
MAST::IsotropicMaterialProperty::
StiffnessMatrix2D::derivative (  const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    // libmesh_assert(_plane_stress); // currently only implemented for plane stress
    RealMatrixX dm;
    m = RealMatrixX::Zero(3,3); dm = RealMatrixX::Zero(3, 3);
    Real E, nu, dEdf, dnudf;
    _E  (p, t, E);   _E.derivative( f, p, t, dEdf);
    _nu (p, t, nu); _nu.derivative( f, p, t, dnudf);
    
    if (_plane_stress)
    {
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
    else // plane_strain
    {
        // parM/parE * parE/parf
        for (unsigned int i=0; i<2; i++) {
            for (unsigned int j=0; j<2; j++)
                if (i == j) // diagonal: direct stress
                    dm(i,i) = 1.*(1.-nu)/((1.+nu)*(1.-2.*nu));
                else // offdiagonal: direct stress
                    dm(i,j) = 1.*nu/((1.+nu)*(1.-2.*nu));
        }
        dm(2,2) = 1./2./(1.+nu); // diagonal: shear stress
        m += dEdf*dm;
        
        
        // parM/parnu * parnu/parf
        for (unsigned int i=0; i<2; i++) {
            for (unsigned int j=0; j<2; j++)
                if (i == j) // diagonal: direct stress
                    dm(i,i) = E*2.*nu*(2.-nu)/pow(2.*nu*nu+nu-1., 2);
                else // offdiagonal: direct stress
                    dm(i,j) = (2.*nu*nu*E+E)/pow(2*nu*nu+nu-1., 2);
        }
        dm(2,2) = -E/2./pow(1.+nu, 2); // diagonal: shear stress
        m += dnudf*dm;
    }
}



MAST::IsotropicMaterialProperty::
StiffnessMatrix3D::StiffnessMatrix3D(const MAST::FieldFunction<Real>& E,
                                     const MAST::FieldFunction<Real>& nu):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix3D"),
_E(E),
_nu(nu) {
    
    _functions.insert(&E);
    _functions.insert(&nu);
}





void
MAST::IsotropicMaterialProperty::
StiffnessMatrix3D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m = RealMatrixX::Zero(6,6);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
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
StiffnessMatrix3D::derivative (  const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(6,6); dm = RealMatrixX::Zero(6,6);
    Real E, nu, dEdf, dnudf;
    _E  (p, t, E);   _E.derivative( f, p, t, dEdf);
    _nu (p, t, nu); _nu.derivative( f, p, t, dnudf);
    
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
InertiaMatrix3D::InertiaMatrix3D(const MAST::FieldFunction<Real>& rho):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_rho(rho) {
    
    _functions.insert(&rho);
}


void
MAST::IsotropicMaterialProperty::
InertiaMatrix3D::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    m  = RealMatrixX::Zero(3,3);
    Real rho;
    _rho(p, t, rho);
    m(0,0) = rho;
    m(1,1) = rho;
    m(2,2) = rho;
}


void
MAST::IsotropicMaterialProperty::
InertiaMatrix3D::derivative(
                            const MAST::FunctionBase &f,
                            const libMesh::Point &p,
                            const Real t,
                            RealMatrixX &m) const {
    
    
    m  = RealMatrixX::Zero(3,3);
    Real rho;
    _rho.derivative( f, p, t, rho);
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
    _alpha(p, t, alpha);
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
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real alpha;
    _alpha.derivative( f, p, t, alpha);
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
    _cp  (p, t, cp);
    _rho (p, t, rho);
    
    m.setZero(1,1);
    
    m(0,0) = cp*rho;
}





void
MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix::
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real cp, dcp, rho, drho;
    _cp  (p, t, cp);    _cp.derivative( f, p, t, dcp);
    _rho (p, t, rho);  _rho.derivative( f, p, t, drho);
    
    m.setZero(1,1);
    
    m(0,0) = dcp*rho + cp*drho;
}




void
MAST::IsotropicMaterialProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real k;
    _k  (p, t, k);
    
    m.setIdentity(_dim, _dim);
    
    m *= k;
}





void
MAST::IsotropicMaterialProperty::ThermalConductanceMatrix::
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real k;
    _k.derivative( f, p, t, k);
    
    m.setIdentity(_dim, _dim);
    
    m *= k;
}





MAST::IsotropicMaterialPropertyCard::IsotropicMaterialPropertyCard():
MAST::MaterialPropertyCardBase (),
_stiff_mat_1d                  (nullptr),
_stiff_mat_2d                  (nullptr),
_stiff_mat_3d                  (nullptr),
_damp_mat                      (nullptr),
_inertia_mat_3d                (nullptr),
_thermal_exp_mat_1d            (nullptr),
_thermal_exp_mat_2d            (nullptr),
_thermal_exp_mat_3d            (nullptr),
_transverse_shear_mat          (nullptr),
_thermal_capacitance_mat_1d    (nullptr),
_thermal_capacitance_mat_2d    (nullptr),
_thermal_capacitance_mat_3d    (nullptr),
_thermal_conductance_mat_1d    (nullptr),
_thermal_conductance_mat_2d    (nullptr),
_thermal_conductance_mat_3d    (nullptr)
{ }



MAST::IsotropicMaterialPropertyCard::~IsotropicMaterialPropertyCard() {
    
    if (_stiff_mat_1d)                delete _stiff_mat_1d;
    if (_stiff_mat_2d)                delete _stiff_mat_2d;
    if (_stiff_mat_3d)                delete _stiff_mat_3d;
    
    if (_damp_mat)                    delete _damp_mat;
    if (_inertia_mat_3d)              delete _inertia_mat_3d;
    if (_thermal_exp_mat_1d)          delete _thermal_exp_mat_1d;
    if (_thermal_exp_mat_2d)          delete _thermal_exp_mat_2d;
    if (_thermal_exp_mat_3d)          delete _thermal_exp_mat_3d;
    if (_transverse_shear_mat)        delete _transverse_shear_mat;
    
    if (_thermal_capacitance_mat_1d)  delete _thermal_capacitance_mat_1d;
    if (_thermal_capacitance_mat_2d)  delete _thermal_capacitance_mat_2d;
    if (_thermal_capacitance_mat_3d)  delete _thermal_capacitance_mat_3d;
    
    if (_thermal_conductance_mat_1d)  delete _thermal_conductance_mat_1d;
    if (_thermal_conductance_mat_2d)  delete _thermal_conductance_mat_2d;
    if (_thermal_conductance_mat_3d)  delete _thermal_conductance_mat_3d;
}


const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::stiffness_matrix(const unsigned int dim,
                                                      const bool plane_stress) {
    
    MAST::FieldFunction<RealMatrixX> *rval = nullptr;
    
    switch (dim) {
        case 1: {
            
            if (!_stiff_mat_1d)
                _stiff_mat_1d = new MAST::IsotropicMaterialProperty::StiffnessMatrix1D
                (this->get<MAST::FieldFunction<Real> >("E"),
                 this->get<MAST::FieldFunction<Real> >("nu"));
            return *_stiff_mat_1d;
        }
            break;
            
        case 2: {
            
            if (!_stiff_mat_2d)
                _stiff_mat_2d = new MAST::IsotropicMaterialProperty::StiffnessMatrix2D
                (this->get<MAST::FieldFunction<Real> >("E"),
                 this->get<MAST::FieldFunction<Real> >("nu"),
             plane_stress);
            
            return *_stiff_mat_2d;
        }
            break;
            
        case 3: {
            
            if (!_stiff_mat_3d)
                _stiff_mat_3d = new MAST::IsotropicMaterialProperty::StiffnessMatrix3D
                (this->get<MAST::FieldFunction<Real> >("E"),
                 this->get<MAST::FieldFunction<Real> >("nu"));
            
            return *_stiff_mat_3d;
        }
            break;
            
        default:
            // should not get here
            libmesh_error_msg("Should not get here; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }
}




const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::damping_matrix(const unsigned int dim) {
    
    // make sure that this is not null
    libmesh_assert(false);
    
    return *_damp_mat;
}




const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::inertia_matrix(const unsigned int dim) {
    
    switch (dim) {
        case 3: {
            
            if (!_inertia_mat_3d)
                _inertia_mat_3d = new MAST::IsotropicMaterialProperty::InertiaMatrix3D
                (this->get<MAST::FieldFunction<Real> >("rho"));
            
            return *_inertia_mat_3d;
        }
            break;
            
        case 1:
        case 2:
        default:
            // implemented only for 3D since the 2D and 1D elements
            // calculate it themselves
            libmesh_error_msg("Implemented only for 3D since the 2D and 1D elements calculate it themselves; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }
}




const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::thermal_expansion_matrix(const unsigned int dim) {

    switch (dim) {
        case 3: {
            
            if (!_thermal_exp_mat_3d)
                _thermal_exp_mat_3d =
                new MAST::IsotropicMaterialProperty::ThermalExpansionMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("alpha_expansion"));
            
            return *_thermal_exp_mat_3d;
        }
            break;
            
            
        case 2: {
            
            if (!_thermal_exp_mat_2d)
                _thermal_exp_mat_2d =
                new MAST::IsotropicMaterialProperty::ThermalExpansionMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("alpha_expansion"));
            
            return *_thermal_exp_mat_2d;
        }
            break;

            
        case 1: {
            
            if (!_thermal_exp_mat_1d)
                _thermal_exp_mat_1d =
                new MAST::IsotropicMaterialProperty::ThermalExpansionMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("alpha_expansion"));
            
            return *_thermal_exp_mat_1d;
        }
            break;

 
        default:
            libmesh_error_msg("Should not get here; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }

    
}




const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::transverse_shear_stiffness_matrix() {
    
    if (!_transverse_shear_mat)
        _transverse_shear_mat =
        new MAST::IsotropicMaterialProperty::TransverseShearStiffnessMatrix
        (this->get<MAST::FieldFunction<Real> >("E"),
         this->get<MAST::FieldFunction<Real> >("nu"));
    
    return *_transverse_shear_mat;
}



const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::capacitance_matrix(const unsigned int dim) {
    
    switch (dim) {
            
        case 1: {
            
            if (!_thermal_capacitance_mat_1d)
                _thermal_capacitance_mat_1d =
                new MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_1d;
        }
            break;

            
        case 2: {
            
            if (!_thermal_capacitance_mat_2d)
                _thermal_capacitance_mat_2d =
                new MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_2d;
        }
            break;

            
        case 3: {
            
            if (!_thermal_capacitance_mat_3d)
                _thermal_capacitance_mat_3d =
                new MAST::IsotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_3d;
        }
            break;

            
        default:
            // should not get here
            libmesh_error_msg("Should not get here; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }
}





const MAST::FieldFunction<RealMatrixX>&
MAST::IsotropicMaterialPropertyCard::conductance_matrix(const unsigned int dim) {
    
    switch (dim) {
            
        case 1: {
            
            if (!_thermal_conductance_mat_1d)
                _thermal_conductance_mat_1d =
                new MAST::IsotropicMaterialProperty::ThermalConductanceMatrix
                (dim, this->get<MAST::FieldFunction<Real> >("k_th"));
            
            return *_thermal_conductance_mat_1d;
        }
            break;
            
            
        case 2: {
            
            if (!_thermal_conductance_mat_2d)
                _thermal_conductance_mat_2d =
                new MAST::IsotropicMaterialProperty::ThermalConductanceMatrix
                (dim, this->get<MAST::FieldFunction<Real> >("k_th"));
            
            return *_thermal_conductance_mat_2d;
        }
            break;
            
            
        case 3: {
            
            if (!_thermal_conductance_mat_3d)
                _thermal_conductance_mat_3d =
                new MAST::IsotropicMaterialProperty::ThermalConductanceMatrix
                (dim, this->get<MAST::FieldFunction<Real> >("k_th"));
            
            return *_thermal_conductance_mat_3d;
        }
            break;
            
            
        default:
            // should not get here
            libmesh_error_msg("Should not get here; " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line number " << __LINE__);
    }
}



