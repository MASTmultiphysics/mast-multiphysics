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
#include "property_cards/orthotropic_element_property_card_3D.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"
#include "coordinates/coordinate_base.h"


namespace MAST {
    
    namespace OrthotropicProperty3D {
        
        
        class StiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                            MAST::CoordinateBase *orient);
            
            StiffnessMatrix(const MAST::OrthotropicProperty3D::StiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _orient(dynamic_cast<MAST::CoordinateBase*>(f._orient->clone().release())) {
                
                _functions.insert(_material_stiffness->master());
                _functions.insert(_orient->master());
            }
            
            /*!
             *   @returns a clone of the functiofn
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::OrthotropicProperty3D::StiffnessMatrix(*this));
            }
            
            virtual ~StiffnessMatrix() {
                delete _material_stiffness;
                delete _orient;
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
            
            MAST::CoordinateBase *_orient;
        };
        
        
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(MAST::FieldFunction<RealMatrixX> *mat);
            
            InertiaMatrix(const MAST::OrthotropicProperty3D::InertiaMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_inertia(f._material_inertia->clone().release()) {
                _functions.insert(_material_inertia->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::OrthotropicProperty3D::InertiaMatrix(*this));
            }
            
            virtual ~InertiaMatrix() { delete _material_inertia;}
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative (const MAST::DerivativeType d,
                                     const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            MAST::FieldFunction<RealMatrixX> *_material_inertia;
        };
        
        
        
        
        class ThermalExpansionMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                                   MAST::FieldFunction<RealMatrixX> *mat_expansion,
                                   MAST::CoordinateBase *orient);
            
            ThermalExpansionMatrix(const MAST::OrthotropicProperty3D::
                                   ThermalExpansionMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()),
            _orient(dynamic_cast<MAST::CoordinateBase*>(f._orient->clone().release())) {
                
                _functions.insert(_material_stiffness->master());
                _functions.insert(_material_expansion->master());
                _functions.insert(_orient->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::OrthotropicProperty3D::ThermalExpansionMatrix(*this));
            }
            
            virtual ~ThermalExpansionMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
                delete _orient;
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
            MAST::CoordinateBase *_orient;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                             MAST::CoordinateBase *orient);
            
            PrestressAMatrix(const MAST::OrthotropicProperty3D::PrestressAMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()),
            _orient(dynamic_cast<MAST::CoordinateBase*>(f._orient->clone().release())){
                _functions.insert(_prestress->master());
                _functions.insert(_orient->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::OrthotropicProperty3D::PrestressAMatrix(*this));
            }
            
            virtual ~PrestressAMatrix() {
                delete _prestress;
                delete _orient;
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
            
            MAST::FieldFunction<RealMatrixX> *_prestress;
            MAST::CoordinateBase *_orient;
        };
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                                     MAST::CoordinateBase *orient);
            
            ThermalConductanceMatrix(const MAST::OrthotropicProperty3D::ThermalConductanceMatrix &f);
            
            
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
            
            MAST::CoordinateBase* _orient;
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond);
            
            ThermalCapacitanceMatrix(const MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix &f);
            
            
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
        };
        
    }
    
}


bool
MAST::OrthotropicElementPropertyCard3D::depends_on(const MAST::FunctionBase& f) const {
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}



MAST::OrthotropicProperty3D::StiffnessMatrix::
StiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat,
                MAST::CoordinateBase *orient):
MAST::FieldFunction<RealMatrixX> ("StiffnessMatrix3D"),
_material_stiffness(mat),
_orient(orient) {
    
    _functions.insert(mat->master());
    _functions.insert(orient->master());
}


void
MAST::OrthotropicProperty3D::
StiffnessMatrix::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    RealMatrixX
    A    = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6);
    
    (*_material_stiffness) (p, t, m);
    (*_orient)             (p, t, A);
    _orient->stress_strain_transformation_matrix(A.transpose(), Tinv);
    
    //    vk'  = vj ej.ek' = vj Ajk = A^T v
    //    v = A vk'
    //    sij' = skl ek.ei' el.ej' = skl Aki Alj = A^T s A
    //    s' = T s
    //    s = Tinv s'
    //    s' = C' e'
    //    T s = C' Rinv T R s
    //    s = Tinv C Rinv T R e
    //    C = Tinv C Rinv T R
    //    T R scales last three columns by 1/2
    //    Rinv T scales last three rows by 2
    //    therefore, Rinv T R scales top right 3x3 block by 1/2,
    //    and bottom left 3x3 block by 2.
    //    Also, Rinv T R = T^{-T}
    m =  Tinv * m * Tinv.transpose();
}



void
MAST::OrthotropicProperty3D::
StiffnessMatrix::derivative (const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {

    RealMatrixX
    dm   = RealMatrixX::Zero(6, 6),
    A    = RealMatrixX::Zero(3, 3),
    dA   = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6),
    dTinv= RealMatrixX::Zero(6, 6);
    
    (*_material_stiffness) (p, t, m);
    _material_stiffness->derivative(d, f, p, t, dm);
    (*_orient)             (p, t, A);
    _orient->stress_strain_transformation_matrix(A.transpose(), Tinv);
    _orient->derivative    (d, f, p, t, dA);
    _orient->stress_strain_transformation_matrix_sens(A.transpose(),
                                                      dA.transpose(),
                                                      dTinv);
    
    m =
    dTinv *  m *  Tinv.transpose() +
    Tinv  * dm *  Tinv.transpose() +
    Tinv  *  m * dTinv.transpose();
}




MAST::OrthotropicProperty3D::
InertiaMatrix::InertiaMatrix(MAST::FieldFunction<RealMatrixX> *mat):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(mat->master());
}



void
MAST::OrthotropicProperty3D::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    RealMatrixX mat;
    m = RealMatrixX::Zero(6, 6);
    // this only returns the material inertia
    (*_material_inertia)(p, t, mat);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            m(i,j) = mat(i,j);
        }
        m(i+3,i+3) = mat(i,i) * 1.0e-12;
    }
}




void
MAST::OrthotropicProperty3D::
InertiaMatrix::derivative (const MAST::DerivativeType d,
                           const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    RealMatrixX mat;
    m  = RealMatrixX::Zero(6, 6);
    // this only returns the material inertia
    // sensitivity of rotary inertia is assumed to be zero
    _material_inertia->derivative(d, f, p, t, mat);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            m(i,j) = mat(i,j);
}




MAST::OrthotropicProperty3D::ThermalExpansionMatrix::
ThermalExpansionMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                       MAST::FieldFunction<RealMatrixX> *mat_expansion,
                       MAST::CoordinateBase *orient):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_orient(orient) {
    
    _functions.insert(mat_stiff->master());
    _functions.insert(mat_expansion->master());
    _functions.insert(orient->master());
}




void
MAST::OrthotropicProperty3D::
ThermalExpansionMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX
    mat  = RealMatrixX::Zero(6,1),
    A    = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6);
    
    (*_material_stiffness) (p, t, m);
    (*_material_expansion)(p, t, mat);
    (*_orient)             (p, t, A);
    _orient->stress_strain_transformation_matrix(A.transpose(), Tinv);
    
    //    epsilon' = T^{-T} epsilon
    //    epsilon  = T^T epsilon'
    //    C = Tinv C' T^{-T}
    //    epsilon_delta_T = C epsilon
    //                    = Tinv C' T^{-T} T^T epsilon'
    //                    = Tinv C' alpha' delta_temp
    //    hence,
    //    mat = Tinv C' alpha'
    m =  Tinv * m * mat;
}






void
MAST::OrthotropicProperty3D::
ThermalExpansionMatrix::derivative (const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    
    RealMatrixX
    dm   = RealMatrixX::Zero(6, 6),
    mat  = RealMatrixX::Zero(6, 1),
    dmat = RealMatrixX::Zero(6, 1),
    A    = RealMatrixX::Zero(3, 3),
    dA   = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6),
    dTinv= RealMatrixX::Zero(6, 6);
    
    (*_material_stiffness) (p, t, m);
    _material_stiffness->derivative(d, f, p, t, dm);
    (*_material_expansion)(p, t, mat);
    _material_expansion->derivative(d, f, p, t, dmat);
    (*_orient)             (p, t, A);
    _orient->stress_strain_transformation_matrix(A.transpose(), Tinv);
    _orient->derivative    (d, f, p, t, dA);
    _orient->stress_strain_transformation_matrix_sens(A.transpose(),
                                                      dA.transpose(),
                                                      dTinv);
    
    m =
    dTinv *  m *  mat +
    Tinv  * dm *  mat +
    Tinv  *  m * dmat;
}




MAST::OrthotropicProperty3D::PrestressAMatrix::
PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress,
                 MAST::CoordinateBase *orient):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix3D"),
_prestress(prestress),
_orient(orient) {
    _functions.insert(prestress->master());
    _functions.insert(orient->master());
}




void
MAST::OrthotropicProperty3D::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    (*_prestress)(p, t, m);
    libmesh_error();
}






void
MAST::OrthotropicProperty3D::
PrestressAMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress->derivative(d, f, p, t, m);
}



MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond,
                         MAST::CoordinateBase *orient):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_orient(orient) {
    
    _functions.insert(mat_cond->master());
    _functions.insert(orient->master());
}


MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::OrthotropicProperty3D::ThermalConductanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cond(f._mat_cond->clone().release()),
_orient(dynamic_cast<MAST::CoordinateBase*>(f._orient->clone().release())) {
    
    _functions.insert(_mat_cond->master());
    _functions.insert(_orient->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicProperty3D::ThermalConductanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalConductanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
~ThermalConductanceMatrix() {
    
    delete _mat_cond;
    delete _orient;
}


void
MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    (*_mat_cond)(p, t, m);
    libmesh_error();
}




void
MAST::OrthotropicProperty3D::ThermalConductanceMatrix::derivative (const MAST::DerivativeType d,
                                                                          const MAST::FunctionBase& f,
                                                                          const libMesh::Point& p,
                                                                          const Real t,
                                                                          RealMatrixX& m) const {
    _mat_cond->derivative(d, f, p, t, m);
}






MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cap):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap) {
    
    _functions.insert(mat_cap->master());
}


MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cap(f._mat_cap->clone().release()) {
    
    _functions.insert(_mat_cap->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() {
    
    delete _mat_cap;
}


void
MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    (*_mat_cap)(p, t, m);
}



void
MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::derivative (const MAST::DerivativeType d,
                                                                          const MAST::FunctionBase& f,
                                                                          const libMesh::Point& p,
                                                                          const Real t,
                                                                          RealMatrixX& m) const {
    _mat_cap->derivative(d, f, p, t, m);
}



MAST::OrthotropicElementPropertyCard3D::
~OrthotropicElementPropertyCard3D() {
    
    if (_orient) delete _orient;
}



void
MAST::OrthotropicElementPropertyCard3D::
set_orientation(const MAST::CoordinateBase& orient) {
    
    _orient = dynamic_cast<MAST::CoordinateBase*>(orient.clone().release());
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::StiffnessMatrix
    (_material->stiffness_matrix(3).release(),
     dynamic_cast<MAST::CoordinateBase*>(_orient->clone().release()));
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
inertia_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::InertiaMatrix
    (_material->inertia_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalExpansionMatrix
    (_material->stiffness_matrix(3).release(),
     _material->thermal_expansion_matrix(3).release(),
     dynamic_cast<MAST::CoordinateBase*>(_orient->clone().release()));
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    // for 3D elements, there is no difference between the A and B matrices
    return this->thermal_expansion_A_matrix(e);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
prestress_A_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
prestress_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalConductanceMatrix
    (_material->conductance_matrix(3).release(),
     dynamic_cast<MAST::CoordinateBase*>(_orient->clone().release()));
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



