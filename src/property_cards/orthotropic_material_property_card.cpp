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
#include "property_cards/orthotropic_material_property_card.h"
#include "base/field_function_base.h"


namespace MAST {
    namespace OrthotropicMaterialProperty {
        
        
        class StiffnessMatrix1D: public MAST::FieldFunction<RealMatrixX> {
        public:
            
            StiffnessMatrix1D(const MAST::FieldFunction<Real>& E11,
                              const MAST::FieldFunction<Real>& G12);
            
            virtual ~StiffnessMatrix1D() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _E11;
            const MAST::FieldFunction<Real>& _G12;
        };
        
        
        
        class TransverseShearStiffnessMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            TransverseShearStiffnessMatrix(const MAST::FieldFunction<Real>& G12);
            
            
            virtual ~TransverseShearStiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real>& _G12;
        };
        
        
        class StiffnessMatrix2D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix2D(const MAST::FieldFunction<Real>& E11,
                              const MAST::FieldFunction<Real>& E22,
                              const MAST::FieldFunction<Real>& E33,
                              const MAST::FieldFunction<Real>& nu12,
                              const MAST::FieldFunction<Real>& nu23,
                              const MAST::FieldFunction<Real>& nu31,
                              const MAST::FieldFunction<Real>& G12,
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
            
            const MAST::FieldFunction<Real>& _E11;
            const MAST::FieldFunction<Real>& _E22;
            const MAST::FieldFunction<Real>& _E33;
            const MAST::FieldFunction<Real>& _nu12;
            const MAST::FieldFunction<Real>& _nu23;
            const MAST::FieldFunction<Real>& _nu31;
            const MAST::FieldFunction<Real>& _G12;
            bool _plane_stress;
        };
        
        
        
        class StiffnessMatrix3D: public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix3D(const MAST::FieldFunction<Real>& E11,
                              const MAST::FieldFunction<Real>& E22,
                              const MAST::FieldFunction<Real>& E33,
                              const MAST::FieldFunction<Real>& nu12,
                              const MAST::FieldFunction<Real>& nu31,
                              const MAST::FieldFunction<Real>& nu23,
                              const MAST::FieldFunction<Real>& G12,
                              const MAST::FieldFunction<Real>& G13,
                              const MAST::FieldFunction<Real>& G23);
            
            
            virtual ~StiffnessMatrix3D() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative(   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<Real> &_E11,  &_E22,  &_E33;
            const MAST::FieldFunction<Real> &_nu12, &_nu31, &_nu23;
            const MAST::FieldFunction<Real> &_G12,  &_G13,  &_G23;

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
            
            ThermalExpansionMatrix(const MAST::FieldFunction<Real>& alpha11):
            MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix"),
            _dim(1),
            _alpha(_dim) {

                _alpha[0] = &alpha11;
                _functions.insert(_alpha[0]);
            }

            
            ThermalExpansionMatrix(const MAST::FieldFunction<Real>& alpha11,
                                   const MAST::FieldFunction<Real>& alpha22):
            MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix"),
            _dim(2),
            _alpha(_dim) {
                
                _alpha[0] = &alpha11;
                _alpha[1] = &alpha22;
                for (unsigned int i=0; i<_dim; i++)
                    _functions.insert(_alpha[i]);
            }

            
            ThermalExpansionMatrix(const MAST::FieldFunction<Real>& alpha11,
                                   const MAST::FieldFunction<Real>& alpha22,
                                   const MAST::FieldFunction<Real>& alpha33):
            MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix"),
            _dim(3),
            _alpha(_dim) {
                
                _alpha[0] = &alpha11;
                _alpha[1] = &alpha22;
                _alpha[2] = &alpha33;
                
                for (unsigned int i=0; i<_dim; i++)
                    _functions.insert(_alpha[i]);
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
            
            std::vector<const MAST::FieldFunction<Real>*> _alpha;
        };
        
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            
            ThermalConductanceMatrix(const MAST::FieldFunction<Real>& k11):
            MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
            _dim(1),
            _k(_dim) {

                _k[0] = &k11;

                _functions.insert(_k[0]);
            }

            ThermalConductanceMatrix(const MAST::FieldFunction<Real>& k11,
                                     const MAST::FieldFunction<Real>& k22):
            MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
            _dim(2),
            _k(_dim) {
                
                _k[0] = &k11;
                _k[1] = &k22;
                
                for (unsigned int i=0; i<_dim; i++)
                    _functions.insert(_k[i]);
            }

            ThermalConductanceMatrix(const MAST::FieldFunction<Real>& k11,
                                     const MAST::FieldFunction<Real>& k22,
                                     const MAST::FieldFunction<Real>& k33):
            MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
            _dim(3),
            _k(_dim) {
                
                _k[0] = &k11;
                _k[1] = &k22;
                _k[2] = &k33;
                
                for (unsigned int i=0; i<_dim; i++)
                    _functions.insert(_k[i]);
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
            
            std::vector<const MAST::FieldFunction<Real>*> _k;
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



MAST::OrthotropicMaterialProperty::
StiffnessMatrix1D::StiffnessMatrix1D(const MAST::FieldFunction<Real>& E11,
                                     const MAST::FieldFunction<Real>& G12):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix1D"),
_E11(E11),
_G12(G12)
{
    _functions.insert(&E11);
    _functions.insert(&G12);
}




void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix1D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m  = RealMatrixX::Zero(2,2);
    Real E11, G12;
    _E11(p, t, E11); 
    _G12(p, t, G12);
    m(0,0) = E11;
    m(1,1) = G12;
}


void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix1D::derivative( const MAST::FunctionBase &f,
                              const libMesh::Point &p,
                              const Real t,
                              RealMatrixX &m) const {
    
    
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2); dm = RealMatrixX::Zero(2,2);
    Real E11, G12, dE11df,dG12df;
    _E11(p, t, E11);     _E11.derivative(f, p, t, dE11df);
    _G12(p, t, G12);    _G12.derivative(f, p, t, dG12df);
    
    m(0,0) = dE11df;
    m(1,1) = dG12df;
}


MAST::OrthotropicMaterialProperty::
TransverseShearStiffnessMatrix::TransverseShearStiffnessMatrix(const MAST::FieldFunction<Real>& G12):
MAST::FieldFunction<RealMatrixX>("TransverseShearStiffnessMatrix"),
_G12(G12)
{
    _functions.insert(&G12);
}



void
MAST::OrthotropicMaterialProperty::
TransverseShearStiffnessMatrix::operator() (const libMesh::Point& p,
                                            const Real t,
                                            RealMatrixX& m) const {
    // TODO: Allow for use of G1z and G2z for compatibility with Nastran input.
    /*!
     *  m = [G1z,  0.
     *        0., G2z];
     * 
     * If G1z and G2z are not available from test data, then G12 is a valid
     * approximate value for these.
     * 
     * Reference: Nastran PSHELL and MAT8 documentation.
     */
    m = RealMatrixX::Zero(2,2);
    Real G12;
    _G12(p, t, G12);
    m(0,0) = G12;
    m(1,1) = G12;
}



void
MAST::OrthotropicMaterialProperty::
TransverseShearStiffnessMatrix::derivative(const MAST::FunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(2,2);      dm = RealMatrixX::Zero(2, 2);
    
    Real G12, dG12;
    _G12(p, t, G12);        _G12.derivative(f, p, t, dG12);
    m(0,0) = dG12;
    m(1,1) = dG12;
    
    
//     // parM/parE * parE/parf
//     dm(0,0) = 1./2./(1.+nu)*kappa;
//     dm(1,1) = dm(0,0);
//     m += dEdf * dm;
//     
//     // parM/parnu * parnu/parf
//     dm(0,0) = -E/2./pow(1.+nu,2)*kappa;
//     dm(1,1) = dm(0,0);
//     m+= dnudf*dm;
//     
//     // parM/parnu * parkappa/parf
//     
//     dm(0,0) = G; dm(1,1) = G;
//     dm += dkappadf*dm;
}




MAST::OrthotropicMaterialProperty::
StiffnessMatrix2D::StiffnessMatrix2D(const MAST::FieldFunction<Real>& E11,
                                     const MAST::FieldFunction<Real>& E22,
                                     const MAST::FieldFunction<Real>& E33,
                                     const MAST::FieldFunction<Real>& nu12,
                                     const MAST::FieldFunction<Real>& nu23,
                                     const MAST::FieldFunction<Real>& nu31,
                                     const MAST::FieldFunction<Real>& G12,
                                     bool plane_stress ):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix2D"),
_E11(E11),
_E22(E22),
_E33(E33),
_nu12(nu12),
_nu23(nu23),
_nu31(nu31),
_G12(G12),
_plane_stress(plane_stress)
{
    _functions.insert(&E11);
    _functions.insert(&E22);
    _functions.insert(&nu12);
    _functions.insert(&G12);
}




void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix2D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    // https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
    // Applying plane strain and plane stress assumptions to the 3D matrix at
    // the website above, results in the same matrix for both. So, I think
    // plane stress, and plane strain, may be equivalent for orthotropic 
    // materials.
    // libmesh_assert(_plane_stress); // currently only implemented for plane stress
    m = RealMatrixX::Zero(3,3);
    Real E11, E22, E33, nu12, nu23, nu31, nu21, G12, D;
    _E11  (p, t,  E11);
    _E22  (p, t,  E22);
    _E33  (p, t,  E33);
    _nu12 (p, t, nu12);
    _nu23 (p, t, nu23);
    _nu31 (p, t, nu31);
    _G12  (p, t,  G12);
    nu21 = nu12 * E22/E11;
    
    if (_plane_stress) // plane stress
    {
        D = (1. - nu12*nu21);

        m(0,0)  =      E11;
        m(0,1)  = nu21*E11;
        
        m(1,0)  = nu12*E22;
        m(1,1)  =      E22;
        
        m.topLeftCorner(2, 2) /= D;
        
        m(2,2)  = G12;
    }
    else // plane strain
    {
        m(0,0) = -((E11*E11)*E33*(E22-E33*(nu23*nu23)))/((E22*E22)*E33
                *(nu12*nu12)+E11*(E33*E33)*(nu23*nu23)+(E11*E11)*E22
                *(nu31*nu31)-E11*E22*E33+E11*E22*E33*nu12*nu23*nu31*2.0);
        
        m(1,1) = -(E11*(E22*E22)*(E33-E11*(nu31*nu31)))/((E22*E22)*E33
                *(nu12*nu12)+E11*(E33*E33)*(nu23*nu23)+(E11*E11)*E22
                *(nu31*nu31)-E11*E22*E33+E11*E22*E33*nu12*nu23*nu31*2.0);
        
        m(2,2) = G12;
        
        m(0,1) = m(1,0) = -(E11*E22*E33*(E22*nu12+E11*nu23*nu31))
                        /((E22*E22)*E33*(nu12*nu12)+E11*(E33*E33)*(nu23*nu23)
                        +(E11*E11)*E22*(nu31*nu31)-E11*E22*E33+E11*E22*E33
                        *nu12*nu23*nu31*2.0);
    }
}




void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix2D::derivative (  const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    // https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
    // Applying plane strain and plane stress assumptions to the 3D matrix at
    // the website above, results in the same matrix for both. So, I think
    // plane stress, and plane strain, may be equivalent for orthotropic 
    // materials.
    // libmesh_assert(_plane_stress); // currently only implemented for plane stress
    RealMatrixX dm;
    m = RealMatrixX::Zero(3,3); dm = RealMatrixX::Zero(3, 3);
    Real E11, E22, E33, nu12, nu23, nu31, nu21, D;
    Real dE11df, dE22df, dE33df, dnu12df, dnu23df, dnu31df, dnu21df, dDdf, dG12df;
    _E11   (p, t,  E11);   _E11.derivative( f, p, t, dE11df);
    _E22   (p, t,  E22);   _E22.derivative( f, p, t, dE22df);
    _E33   (p, t,  E33);   _E33.derivative( f, p, t, dE33df);
    _nu12  (p, t, nu12);  _nu12.derivative( f, p, t, dnu12df);
    _nu23  (p, t, nu12);  _nu23.derivative( f, p, t, dnu12df);
    _nu31  (p, t, nu12);  _nu31.derivative( f, p, t, dnu12df);
    _G12.derivative( f, p, t, dG12df);
    
    if (_plane_stress)
    {
        nu21    = nu12 * E22/E11;
        dnu21df = dnu12df * E22/E11 + nu12 * dE22df/E11 - nu12 * E22/E11/E11*dE11df;
        D    = (1. - nu12*nu21);
        dDdf = (- dnu12df*nu21 - nu12*dnu21df);
        
        m(0,0)  =      E11;
        m(0,1)  = nu21*E11;
        
        m(1,0)  = nu12*E22;
        m(1,1)  =      E22;
        
        m.topLeftCorner(2, 2) *= -dDdf/D/D;
        m(2,2)  = dG12df;
        
        dm(0,0)  =      dE11df;
        dm(0,1)  = nu21*dE11df + dnu21df*E11;
        
        dm(1,0)  = nu12*dE22df + dnu12df*E22;
        dm(1,1)  =      dE22df;

        m.topLeftCorner(2, 2) += dm.topLeftCorner(3, 3)/D;
    }
    else
    {
        m(0,0) = 1.0/pow((E22*E22)*E33*(nu12*nu12)+E11*(E33*E33)*(nu23*nu23) 
            + (E11*E11)*E22*(nu31*nu31)-E11*E22*E33+E11*E22*E33*nu12*nu23*nu31*2.0,2.0)
            *(E11*(E11*(E22*E22)*(E33*E33)*dE11df+E11*(E33*E33*E33*E33)*dE11df
            *(nu23*nu23*nu23*nu23)-E11*E22*(E33*E33*E33)*dE11df*(nu23*nu23)*2.0)
            +nu31*(E11*((E11*E11*E11)*(E22*E22)*E33*dnu31df*2.0+(E11*E11)*(E22*E22)
            *(E33*E33)*dnu12df*nu23*2.0-(E11*E11)*E22*(E33*E33*E33)*dnu12df
            *(nu23*nu23*nu23)*2.0-(E11*E11*E11)*E22*(E33*E33)*dnu31df
            *(nu23*nu23)*2.0)+E11*nu12*((E11*E11)*(E22*E22)*(E33*E33)*dnu23df*2.0
            -(E11*E11)*(E33*E33*E33)*dE22df*(nu23*nu23*nu23)*2.0+(E11*E11)*E22*(E33*E33)
            *dE33df*(nu23*nu23*nu23)*2.0+(E11*E11)*E22*(E33*E33*E33)*dnu23df
            *(nu23*nu23)*2.0-E11*(E22*E22)*(E33*E33)*dE11df*nu23*2.0+E11*E22
            *(E33*E33*E33)*dE11df*(nu23*nu23*nu23)*2.0))+E11*nu12*(E11*(E22*E22*E22)
            *(E33*E33)*dnu12df*2.0-E11*(E22*E22)*(E33*E33*E33)*dnu12df*(nu23*nu23)
            *2.0+(E11*E11)*(E22*E22)*(E33*E33)*dnu31df*nu23*2.0-(E11*E11)*E22*(E33*E33*E33)
            *dnu31df*(nu23*nu23*nu23)*2.0)-E11*(nu31*nu31)*((E11*E11*E11)*(E22*E22)
            *dE33df+(E11*E11*E11)*(E33*E33)*dE22df*(nu23*nu23)-(E11*E11*E11)*E22*E33
            *dE33df*(nu23*nu23)*2.0-(E11*E11*E11)*E22*(E33*E33)*dnu23df*nu23*2.0)+
            E11*(nu12*nu12)*((E22*E22*E22)*(E33*E33)*dE11df*-2.0+E11*(E22*E22)*(E33*E33)
            *dE22df+(E22*E22)*(E33*E33*E33)*dE11df*(nu23*nu23)*2.0+E11*(E22*E22)
            *(E33*E33)*dE33df*(nu23*nu23)-E11*E22*(E33*E33*E33)*dE22df*(nu23*nu23)
            *2.0+E11*(E22*E22)*(E33*E33*E33)*dnu23df*nu23*2.0));
    
        m(1,1) = 1.0/pow((E22*E22)*E33*(nu12*nu12)+E11*(E33*E33)*(nu23*nu23)+(E11*E11)
            *E22*(nu31*nu31)-E11*E22*E33+E11*E22*E33*nu12*nu23*nu31*2.0,2.0)
            *(E22*((E11*E11)*E22*(E33*E33)*dE22df+(E11*E11*E11*E11)*E22*dE22df
            *(nu31*nu31*nu31*nu31)-(E11*E11*E11)*E22*E33*dE22df*(nu31*nu31)*2.0)
            +nu23*(E22*((E11*E11)*E22*(E33*E33*E33)*dnu23df*2.0+(E11*E11)*(E22*E22)
            *(E33*E33)*dnu12df*nu31*2.0-(E11*E11*E11)*(E22*E22)*E33*dnu12df
            *(nu31*nu31*nu31)*2.0-(E11*E11*E11)*E22*(E33*E33)*dnu23df*(nu31*nu31)
            *2.0)+E22*nu12*((E11*E11)*(E22*E22)*(E33*E33)*dnu31df*2.0-(E11*E11*E11)
            *(E22*E22)*dE33df*(nu31*nu31*nu31)*2.0+(E11*E11)*(E22*E22)*E33*dE11df
            *(nu31*nu31*nu31)*2.0+(E11*E11*E11)*(E22*E22)*E33*dnu31df*(nu31*nu31)*2.0
            -(E11*E11)*E22*(E33*E33)*dE22df*nu31*2.0+(E11*E11*E11)*E22*E33*dE22df
            *(nu31*nu31*nu31)*2.0))+E22*nu12*(E11*(E22*E22*E22)*(E33*E33)*dnu12df*2.0
            -(E11*E11)*(E22*E22*E22)*E33*dnu12df*(nu31*nu31)*2.0+(E11*E11)*(E22*E22)
            *(E33*E33)*dnu23df*nu31*2.0-(E11*E11*E11)*(E22*E22)*E33*dnu23df
            *(nu31*nu31*nu31)*2.0)-E22*(nu12*nu12)*((E22*E22*E22)*(E33*E33)*dE11df
            +(E11*E11)*(E22*E22*E22)*dE33df*(nu31*nu31)-E11*(E22*E22*E22)*E33*dE11df
            *(nu31*nu31)*2.0-(E11*E11)*(E22*E22*E22)*E33*dnu31df*nu31*2.0)+E22
            *(nu23*nu23)*((E11*E11)*(E33*E33*E33)*dE22df*-2.0+(E11*E11)*E22*(E33*E33)
            *dE33df+(E11*E11*E11)*(E33*E33)*dE22df*(nu31*nu31)*2.0+(E11*E11)*E22
            *(E33*E33)*dE11df*(nu31*nu31)-(E11*E11*E11)*E22*E33*dE33df*(nu31*nu31)
            *2.0+(E11*E11*E11)*E22*(E33*E33)*dnu31df*nu31*2.0));
        
        m(2,2) = dG12df;
        
        m(0,1) = m(1,0) = -1.0/pow((E22*E22)*E33*(nu12*nu12)+E11*(E33*E33)
            *(nu23*nu23)+(E11*E11)*E22*(nu31*nu31)-E11*E22*E33+E11*E22*E33*nu12*nu23
            *nu31*2.0,2.0)*(-nu31*((E33*E33)*((E11*E11*E11)*(E22*E22)*dnu23df+(E11*E11)
            *(E22*E22)*dE11df*nu23+(E11*E11*E11)*E22*dE33df*(nu23*nu23*nu23))
            -(E33*E33*E33)*((E11*E11*E11)*dE22df*(nu23*nu23*nu23)+(E11*E11)*E22*dE11df
            *(nu23*nu23*nu23)-(E11*E11*E11)*E22*dnu23df*(nu23*nu23)))
            +(nu31*nu31*nu31)*((E11*E11*E11*E11)*(E22*E22)*E33*dnu23df+(E11*E11*E11*E11)
            *(E22*E22)*dE33df*nu23)-(nu12*nu12)*((E33*E33)*(E11*(E22*E22*E22*E22)
            *dnu12df+(E11*E11)*(E22*E22*E22)*dnu31df*nu23)-(E33*E33)*nu31*(-(E11*E11)
            *(E22*E22*E22)*dnu23df+E11*(E22*E22*E22)*dE11df*nu23*2.0+(E11*E11)*(E22*E22)
            *dE22df*nu23))+(nu31*nu31)*(E33*((E11*E11*E11)*(E22*E22*E22)*dnu12df
            -(E11*E11*E11*E11)*(E22*E22)*dnu31df*nu23)-(E11*E11*E11)*(E22*E22)*(E33*E33)
            *dnu12df*(nu23*nu23)*2.0)+nu12*((nu31*nu31)*(-E33*((E11*E11)
            *(E22*E22*E22)*dE11df-(E11*E11*E11)*(E22*E22)*dE22df)+(E11*E11*E11)*(E22*E22*E22)
            *dE33df+(E11*E11)*(E22*E22)*(E33*E33)*dE11df*(nu23*nu23)*2.0)-(E33*E33)
            *((E11*E11)*(E22*E22)*dE22df+(E11*E11)*(E22*E22)*dE33df*(nu23*nu23))
            +(E33*E33*E33)*((E11*E11)*E22*dE22df*(nu23*nu23)*2.0-(E11*E11)*(E22*E22)
            *dnu23df*nu23*2.0)-nu31*((E11*E11*E11)*(E22*E22*E22)*E33*dnu31df*2.0
            +(E11*E11)*(E22*E22*E22)*(E33*E33)*dnu12df*nu23*2.0))+(E33*E33*E33)
            *((E11*E11*E11)*E22*dnu31df*(nu23*nu23*nu23)+(E11*E11)*(E22*E22)*dnu12df
            *(nu23*nu23))-(E33*E33)*((E11*E11)*(E22*E22*E22)*dnu12df+(E11*E11*E11)
            *(E22*E22)*dnu31df*nu23)+(E22*E22*E22*E22)*(E33*E33)*dE11df
            *(nu12*nu12*nu12));
    }
}



MAST::OrthotropicMaterialProperty::
StiffnessMatrix3D::StiffnessMatrix3D(const MAST::FieldFunction<Real>& E11,
                                     const MAST::FieldFunction<Real>& E22,
                                     const MAST::FieldFunction<Real>& E33,
                                     const MAST::FieldFunction<Real>& nu12,
                                     const MAST::FieldFunction<Real>& nu31,
                                     const MAST::FieldFunction<Real>& nu23,
                                     const MAST::FieldFunction<Real>& G12,
                                     const MAST::FieldFunction<Real>& G13,
                                     const MAST::FieldFunction<Real>& G23):
MAST::FieldFunction<RealMatrixX>("StiffnessMatrix3D"),
_E11(E11),
_E22(E22),
_E33(E33),
_nu12(nu12),
_nu31(nu31),
_nu23(nu23),
_G12(G12),
_G13(G13),
_G23(G23)
{
    _functions.insert(&E11);
    _functions.insert(&E22);
    _functions.insert(&E33);
    _functions.insert(&nu12);
    _functions.insert(&nu31);
    _functions.insert(&nu23);
    _functions.insert(&G12);
    _functions.insert(&G13);
    _functions.insert(&G23);
}





void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix3D::operator() (const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    m = RealMatrixX::Zero(6,6);
    Real E11, E22, E33, G12, G13, G23, nu12, nu13, nu23, nu21, nu31, nu32, D;
    _E11  (p, t,  E11);
    _E22  (p, t,  E22);
    _E33  (p, t,  E33);
    _nu12(p, t, nu12);
    _nu31(p, t, nu31);
    _nu23(p, t, nu23);
    _G12  (p, t,  G12);
    _G13  (p, t,  G13);
    _G23  (p, t,  G23);
    nu21 = nu12 * E22/E11;
    nu13 = nu31 * E11/E33;
    nu32 = nu23 * E33/E22;
    D = 1.- nu12*nu21 - nu13*nu31 - nu23*nu32 - nu12*nu23*nu31 - nu13*nu21*nu32;
    

    m(0,0) = (  1. - nu23*nu32)*E11;
    m(0,1) = (nu21 + nu23*nu31)*E11;
    m(0,2) = (nu31 + nu21*nu32)*E11;

    m(1,0) = (nu12 + nu13*nu32)*E22;
    m(1,1) = (  1. - nu13*nu31)*E22;
    m(1,2) = (nu32 + nu12*nu31)*E22;

    m(2,0) = (nu13 + nu12*nu23)*E33;
    m(2,1) = (nu23 + nu13*nu21)*E33;
    m(2,2) = (  1. - nu12*nu21)*E33;

    m.topLeftCorner(3, 3) /= D;
    m(3,3) = G12;
    m(4,4) = G23;
    m(5,5) = G13;
}




void
MAST::OrthotropicMaterialProperty::
StiffnessMatrix3D::derivative (  const MAST::FunctionBase& f,
                               const libMesh::Point& p,
                               const Real t,
                               RealMatrixX& m) const {
    RealMatrixX dm;
    m = RealMatrixX::Zero(6,6); dm = RealMatrixX::Zero(6,6);
    Real
    E11, dE11df,
    E22, dE22df,
    E33, dE33df,
    dG12df,
    dG13df,
    dG23df,
    nu12, dnu12df,
    nu13, dnu13df,
    nu23, dnu23df,
    nu21, dnu21df,
    nu31, dnu31df,
    nu32, dnu32df,
    D, dDdf;
    _E11   (p, t,  E11);   _E11.derivative( f, p, t, dE11df);
    _E22   (p, t,  E22);   _E22.derivative( f, p, t, dE22df);
    _E33   (p, t,  E33);   _E33.derivative( f, p, t, dE33df);
    _nu12  (p, t, nu12); _nu12.derivative( f, p, t, dnu12df);
    _nu31  (p, t, nu31); _nu31.derivative( f, p, t, dnu31df);
    _nu23  (p, t, nu23); _nu23.derivative( f, p, t, dnu23df);
    _G12.derivative( f, p, t,  dG12df);
    _G13.derivative( f, p, t,  dG13df);
    _G23.derivative( f, p, t,  dG23df);
    nu21    = nu12 * E22/E11;
    dnu21df = dnu12df * E22/E11 + nu12 * dE22df/E11 - nu12 * E22/E11/E11*dE11df;
    nu13    = nu31 * E11/E33;
    dnu13df = dnu31df * E11/E33 + nu31 * dE11df/E33 - nu31 * E11/E33/E33*dE33df;
    nu32    = nu23 * E33/E22;
    dnu32df = dnu23df * E33/E22 + nu23 * dE33df/E22 - nu23 * E33/E22/E22*dE22df;
    D    = 1.- nu12*nu21 - nu13*nu31 - nu23*nu32 - nu12*nu23*nu31 - nu13*nu21*nu32;
    dDdf =
    - dnu12df*nu21 - nu12*dnu21df
    - dnu13df*nu31 - nu13*dnu31df
    - dnu23df*nu32 - nu23*dnu32df
    - dnu12df*nu23*nu31
    - nu12*dnu23df*nu31
    - nu12*nu23*dnu31df
    - dnu13df*nu21*nu32
    - nu13*dnu21df*nu32
    - nu13*nu21*dnu32df;
    
    m(0,0) = (  1. - nu23*nu32)*E11;
    m(0,1) = (nu21 + nu23*nu31)*E11;
    m(0,2) = (nu31 + nu21*nu32)*E11;
    
    m(1,0) = (nu12 + nu13*nu32)*E22;
    m(1,1) = (  1. - nu13*nu31)*E22;
    m(1,2) = (nu32 + nu12*nu31)*E22;
    
    m(2,0) = (nu13 + nu12*nu23)*E33;
    m(2,1) = (nu23 + nu13*nu21)*E33;
    m(2,2) = (  1. - nu12*nu21)*E33;
    m *= -dDdf/D/D;

    m(3,3) = dG12df;
    m(4,4) = dG23df;
    m(5,5) = dG13df;


    dm(0,0) = (- dnu23df*nu32 - nu23*dnu32df)*E11 + (  1. - nu23*nu32)*dE11df;
    dm(0,1) = (dnu21df + dnu23df*nu31 + nu23*dnu31df)*E11 + (nu21 + nu23*nu31)*dE11df;
    dm(0,2) = (dnu31df + dnu21df*nu32 + nu21*dnu32df)*E11 + (nu31 + nu21*nu32)*dE11df;
    
    dm(1,0) = (dnu12df + dnu13df*nu32 + nu13*dnu32df)*E22 + (nu12 + nu13*nu32)*dE22df;
    dm(1,1) = (- dnu13df*nu31 - nu13*dnu31df)*E22 + (  1. - nu13*nu31)*dE22df;
    dm(1,2) = (dnu32df + dnu12df*nu31 + nu12*dnu31df)*E22 + (nu32 + nu12*nu31)*dE22df;
    
    dm(2,0) = (dnu13df + dnu12df*nu23 + nu12*dnu23df)*E33 + (nu13 + nu12*nu23)*dE33df;
    dm(2,1) = (dnu23df + dnu13df*nu21 + nu13*dnu21df)*E33 + (nu23 + nu13*nu21)*dE33df;
    dm(2,2) = (- dnu12df*nu21 - nu12*dnu21df)*E33 + (  1. - nu12*nu21)*dE33df;

    m += dm/D;
}



MAST::OrthotropicMaterialProperty::
InertiaMatrix3D::InertiaMatrix3D(const MAST::FieldFunction<Real>& rho):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_rho(rho) {
    
    _functions.insert(&rho);
}



void
MAST::OrthotropicMaterialProperty::
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
MAST::OrthotropicMaterialProperty::
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
MAST::OrthotropicMaterialProperty::ThermalExpansionMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    

    Real alpha;
    switch (_dim) {
        case 1:
            m.setZero(2, 1);
            break;
            
        case 2:
            m.setZero(3, 1);
            break;
            
        case 3:
            m.setZero(6, 1);
            break;
    }
    
    for (unsigned int i=0; i<_dim; i++) {
        (*_alpha[i])  (p, t, alpha);
        m(i,0) = alpha;
    }
}





void
MAST::OrthotropicMaterialProperty::ThermalExpansionMatrix::
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    
    Real alpha;
    switch (_dim) {
        case 1:
            m.setZero(2, 1);
            break;
            
        case 2:
            m.setZero(3, 1);
            break;
            
        case 3:
            m.setZero(6, 1);
            break;
    }
    
    for (unsigned int i=0; i<_dim; i++) {
        _alpha[i]->derivative( f, p, t, alpha);
        m(i,0) = alpha;
    }
}





void
MAST::OrthotropicMaterialProperty::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real cp, rho;
    _cp   (p, t, cp);
    _rho  (p, t, rho);
    
    m.setZero(1,1);
    
    m(0,0) = cp*rho;
}





void
MAST::OrthotropicMaterialProperty::ThermalCapacitanceMatrix::
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    
    Real cp, dcp, rho, drho;
    _cp   (p, t, cp);    _cp.derivative( f, p, t, dcp);
    _rho  (p, t, rho);  _rho.derivative( f, p, t, drho);
    
    m.setZero(1,1);
    
    m(0,0) = dcp*rho + cp*drho;
}




void
MAST::OrthotropicMaterialProperty::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real k;
    m.setZero(_dim, _dim);
    for (unsigned int i=0; i<_dim; i++) {
        (*_k[i])  (p, t, k);
        m(i,i) = k;
    }
    
}





void
MAST::OrthotropicMaterialProperty::ThermalConductanceMatrix::
derivative (const MAST::FunctionBase& f,
            const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    Real k;
    m.setZero(_dim, _dim);
    for (unsigned int i=0; i<_dim; i++) {
        _k[i]->derivative( f, p, t, k);
        m(i,i) = k;
    }
}



MAST::OrthotropicMaterialPropertyCard::OrthotropicMaterialPropertyCard():
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



MAST::OrthotropicMaterialPropertyCard::~OrthotropicMaterialPropertyCard() {
    
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
MAST::OrthotropicMaterialPropertyCard::stiffness_matrix(const unsigned int dim,
                                                      const bool plane_stress) {
    
    switch (dim) {
            
        case 1: {
            
            if (!_stiff_mat_1d)
                _stiff_mat_1d = new MAST::OrthotropicMaterialProperty::StiffnessMatrix1D
                (this->get<MAST::FieldFunction<Real> >("E11"),
                 this->get<MAST::FieldFunction<Real> >("G12"));
            
            return *_stiff_mat_1d;
        }
            break;
            
        case 2: {
            
            if (!_stiff_mat_2d)
                _stiff_mat_2d = new MAST::OrthotropicMaterialProperty::StiffnessMatrix2D
                (this->get<MAST::FieldFunction<Real> >("E11"),
                 this->get<MAST::FieldFunction<Real> >("E22"),
                 this->get<MAST::FieldFunction<Real> >("E33"),
                 this->get<MAST::FieldFunction<Real> >("nu12"),
                 this->get<MAST::FieldFunction<Real> >("nu23"),
                 this->get<MAST::FieldFunction<Real> >("nu31"),
                 this->get<MAST::FieldFunction<Real> >("G12"),
                 plane_stress);
            
            return *_stiff_mat_2d;
        }
            break;
            
        case 3: {
            
            if (!_stiff_mat_3d)
                _stiff_mat_3d = new MAST::OrthotropicMaterialProperty::StiffnessMatrix3D
                (this->get<MAST::FieldFunction<Real> >("E11"),
                 this->get<MAST::FieldFunction<Real> >("E22"),
                 this->get<MAST::FieldFunction<Real> >("E33"),
                 this->get<MAST::FieldFunction<Real> >("nu12"),
                 this->get<MAST::FieldFunction<Real> >("nu31"),
                 this->get<MAST::FieldFunction<Real> >("nu23"),
                 this->get<MAST::FieldFunction<Real> >("G12"),
                 this->get<MAST::FieldFunction<Real> >("G13"),
                 this->get<MAST::FieldFunction<Real> >("G23"));
            
            return *_stiff_mat_3d;
        }
            break;
            
        default:
            // should not get here
            libmesh_error();
    }
}




const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::damping_matrix(const unsigned int dim) {
    
    // make sure that this is not null
    libmesh_assert(false);
    
    return *_damp_mat;
}




const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::inertia_matrix(const unsigned int dim) {
    
    switch (dim) {
        case 3: {
            
            if (!_inertia_mat_3d)
                _inertia_mat_3d = new MAST::OrthotropicMaterialProperty::InertiaMatrix3D
                (this->get<MAST::FieldFunction<Real> >("rho"));
            
            return *_inertia_mat_3d;
        }
            break;
            
        case 1:
        case 2:
        default:
            // implemented only for 3D since the 2D and 1D elements
            // calculate it themselves
            libmesh_error();
            
    }
}




const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::thermal_expansion_matrix(const unsigned int dim) {

    switch (dim) {
        case 1: {
            
            if (!_thermal_exp_mat_1d)
                _thermal_exp_mat_1d =
                new MAST::OrthotropicMaterialProperty::ThermalExpansionMatrix
                (this->get<MAST::FieldFunction<Real> >("alpha11_expansion"));
            
            return *_thermal_exp_mat_1d;
        }
            break;
            
            
        case 2: {
            
            if (!_thermal_exp_mat_2d)
                _thermal_exp_mat_2d =
                new MAST::OrthotropicMaterialProperty::ThermalExpansionMatrix
                (this->get<MAST::FieldFunction<Real> >("alpha11_expansion"),
                 this->get<MAST::FieldFunction<Real> >("alpha22_expansion"));
            
            return *_thermal_exp_mat_2d;
        }
            break;
            
            
        case 3: {
            
            if (!_thermal_exp_mat_3d)
                _thermal_exp_mat_3d =
                new MAST::OrthotropicMaterialProperty::ThermalExpansionMatrix
                (this->get<MAST::FieldFunction<Real> >("alpha11_expansion"),
                 this->get<MAST::FieldFunction<Real> >("alpha22_expansion"),
                 this->get<MAST::FieldFunction<Real> >("alpha33_expansion"));
            
            return *_thermal_exp_mat_3d;
        }
            break;
            
            
        default:
            libmesh_error();
    }
}




const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::transverse_shear_stiffness_matrix() {
    
    if (!_transverse_shear_mat)
        _transverse_shear_mat =
        new MAST::OrthotropicMaterialProperty::TransverseShearStiffnessMatrix
        (this->get<MAST::FieldFunction<Real> >("G12"));
    
    return *_transverse_shear_mat;

}



const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::capacitance_matrix(const unsigned int dim) {
    
    switch (dim) {
            
        case 1: {
            
            if (!_thermal_capacitance_mat_1d)
                _thermal_capacitance_mat_1d =
                new MAST::OrthotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_1d;
        }
            break;
            
            
        case 2: {
            
            if (!_thermal_capacitance_mat_2d)
                _thermal_capacitance_mat_2d =
                new MAST::OrthotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_2d;
        }
            break;
            
            
        case 3: {
            
            if (!_thermal_capacitance_mat_3d)
                _thermal_capacitance_mat_3d =
                new MAST::OrthotropicMaterialProperty::ThermalCapacitanceMatrix
                (dim,
                 this->get<MAST::FieldFunction<Real> >("rho"),
                 this->get<MAST::FieldFunction<Real> >("cp"));
            
            return *_thermal_capacitance_mat_3d;
        }
            break;
            
            
        default:
            // should not get here
            libmesh_error();
    }
}





const MAST::FieldFunction<RealMatrixX>&
MAST::OrthotropicMaterialPropertyCard::conductance_matrix(const unsigned int dim) {
    
    switch (dim) {
            
        case 1: {
            
            if (!_thermal_conductance_mat_1d)
                _thermal_conductance_mat_1d =
                new MAST::OrthotropicMaterialProperty::ThermalConductanceMatrix
                (this->get<MAST::FieldFunction<Real> >("k11_th"));
            
            return *_thermal_conductance_mat_1d;
        }
            break;
            
            
        case 2: {
            
            if (!_thermal_conductance_mat_2d)
                _thermal_conductance_mat_2d =
                new MAST::OrthotropicMaterialProperty::ThermalConductanceMatrix
                (this->get<MAST::FieldFunction<Real> >("k11_th"),
                 this->get<MAST::FieldFunction<Real> >("k22_th"));
            
            return *_thermal_conductance_mat_2d;
        }
            break;
            
            
        case 3: {
            
            if (!_thermal_conductance_mat_3d)
                _thermal_conductance_mat_3d =
                new MAST::OrthotropicMaterialProperty::ThermalConductanceMatrix
                (this->get<MAST::FieldFunction<Real> >("k11_th"),
                 this->get<MAST::FieldFunction<Real> >("k22_th"),
                 this->get<MAST::FieldFunction<Real> >("k33_th"));
            
            return *_thermal_conductance_mat_3d;
        }
            break;
            
            
        default:
            // should not get here
            libmesh_error();
    }
}



