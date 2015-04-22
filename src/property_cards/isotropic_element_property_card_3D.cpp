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
#include "property_cards/isotropic_element_property_card_3D.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"



namespace MAST {
    
    namespace IsotropicElementProperty3D {
        
        
        class StiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat);
            
            StiffnessMatrix(const MAST::IsotropicElementProperty3D::StiffnessMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()) {
                _functions.insert(_material_stiffness->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicElementProperty3D::StiffnessMatrix(*this));
            }
            
            virtual ~StiffnessMatrix() { delete _material_stiffness;}
            
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
        };
        
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(MAST::FieldFunction<RealMatrixX> *mat);
            
            InertiaMatrix(const MAST::IsotropicElementProperty3D::InertiaMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_inertia(f._material_inertia->clone().release()) {
                _functions.insert(_material_inertia->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicElementProperty3D::InertiaMatrix(*this));
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
                                   MAST::FieldFunction<RealMatrixX> *mat_expansion);
            
            ThermalExpansionMatrix(const MAST::IsotropicElementProperty3D::
                                   ThermalExpansionMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()) {
                _functions.insert(_material_stiffness->master());
                _functions.insert(_material_expansion->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicElementProperty3D::ThermalExpansionMatrix(*this));
            }
            
            virtual ~ThermalExpansionMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
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
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress);
            
            PrestressAMatrix(const MAST::IsotropicElementProperty3D::PrestressAMatrix& f):
            MAST::FieldFunction<RealMatrixX>(f),
            _prestress(f._prestress->clone().release()) {
                _functions.insert(_prestress->master());
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<RealMatrixX> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
                (new MAST::IsotropicElementProperty3D::PrestressAMatrix(*this));
            }
            
            virtual ~PrestressAMatrix() { delete _prestress;}
            
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
        };
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond);
            
            ThermalConductanceMatrix(const MAST::IsotropicElementProperty3D::ThermalConductanceMatrix &f);
            
            
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
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond);
            
            ThermalCapacitanceMatrix(const MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix &f);
            
            
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
MAST::IsotropicElementPropertyCard3D::depends_on(const MAST::FunctionBase& f) const {
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}



MAST::IsotropicElementProperty3D::StiffnessMatrix::
StiffnessMatrix(MAST::FieldFunction<RealMatrixX> *mat):
MAST::FieldFunction<RealMatrixX> ("StiffnessMatrix3D"),
_material_stiffness(mat) {
    _functions.insert(mat->master());
}


void
MAST::IsotropicElementProperty3D::
StiffnessMatrix::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    // this only returns the material stiffness
    (*_material_stiffness)
    (p, t, m);
}



void
MAST::IsotropicElementProperty3D::
StiffnessMatrix::derivative (const MAST::DerivativeType d,
                             const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    // this only returns the material stiffness
    _material_stiffness->derivative(d, f, p, t, m);
}




MAST::IsotropicElementProperty3D::
InertiaMatrix::InertiaMatrix(MAST::FieldFunction<RealMatrixX> *mat):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(mat->master());
}



void
MAST::IsotropicElementProperty3D::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {

    (*_material_inertia)(p, t, m);
}




void
MAST::IsotropicElementProperty3D::
InertiaMatrix::derivative (const MAST::DerivativeType d,
                           const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {

    _material_inertia->derivative(d, f, p, t, m);
}




MAST::IsotropicElementProperty3D::ThermalExpansionMatrix::
ThermalExpansionMatrix(MAST::FieldFunction<RealMatrixX> *mat_stiff,
                       MAST::FieldFunction<RealMatrixX> *mat_expansion):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion) {
    _functions.insert(mat_stiff->master());
    _functions.insert(mat_expansion->master());
}




void
MAST::IsotropicElementProperty3D::
ThermalExpansionMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat;
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, mat);
    m *= mat;
}






void
MAST::IsotropicElementProperty3D::
ThermalExpansionMatrix::derivative (const MAST::DerivativeType d,
                                    const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    libmesh_error();
}




MAST::IsotropicElementProperty3D::PrestressAMatrix::
PrestressAMatrix(MAST::FieldFunction<RealMatrixX> *prestress):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix3D"),
_prestress(prestress){
    _functions.insert(prestress->master());
}




void
MAST::IsotropicElementProperty3D::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    (*_prestress)(p, t, m);
}






void
MAST::IsotropicElementProperty3D::
PrestressAMatrix::derivative (const MAST::DerivativeType d,
                              const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress->derivative(d, f, p, t, m);
}



MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cond):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond) {
    
    _functions.insert(mat_cond->master());
}


MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::IsotropicElementProperty3D::ThermalConductanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cond(f._mat_cond->clone().release()) {
    
    _functions.insert(_mat_cond->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalConductanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
~ThermalConductanceMatrix() {
    
    delete _mat_cond;
}


void
MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    (*_mat_cond)(p, t, m);
}




void
MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::derivative (const MAST::DerivativeType d,
                                                                        const MAST::FunctionBase& f,
                                                                        const libMesh::Point& p,
                                                                        const Real t,
                                                                        RealMatrixX& m) const {
    _mat_cond->derivative(d, f, p, t, m);
}






MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(MAST::FieldFunction<RealMatrixX> *mat_cap):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap) {
    
    _functions.insert(mat_cap->master());
}


MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix &f):
MAST::FieldFunction<RealMatrixX>(f),
_mat_cap(f._mat_cap->clone().release()) {
    
    _functions.insert(_mat_cap->master());
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::clone() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix(*this);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}




MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() {
    
    delete _mat_cap;
}


void
MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    (*_mat_cap)(p, t, m);
}



void
MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::derivative (const MAST::DerivativeType d,
                                                                        const MAST::FunctionBase& f,
                                                                        const libMesh::Point& p,
                                                                        const Real t,
                                                                        RealMatrixX& m) const {
    _mat_cap->derivative(d, f, p, t, m);
}



//void
//MAST::IsotropicElementProperty3D::
//PrestressAMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v) const {
//    libmesh_error();
//}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::StiffnessMatrix
    (_material->stiffness_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
inertia_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::InertiaMatrix
    (_material->inertia_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalExpansionMatrix
    (_material->stiffness_matrix(3).release(),
     _material->thermal_expansion_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    // for 3D elements, there is no difference between the A and B matrices
    return this->thermal_expansion_A_matrix(e);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
prestress_A_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
prestress_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(NULL);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalConductanceMatrix
    (_material->conductance_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(3).release());
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


