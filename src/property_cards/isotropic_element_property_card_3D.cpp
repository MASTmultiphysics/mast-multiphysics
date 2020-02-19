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
#include "property_cards/isotropic_element_property_card_3D.h"
#include "property_cards/material_property_card_base.h"
#include "base/field_function_base.h"



namespace MAST {
    
    namespace IsotropicElementProperty3D {
        
        
        class StiffnessMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            StiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat);
            
            virtual ~StiffnessMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_stiffness;
        };
        
        
        
        class InertiaMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            InertiaMatrix(const MAST::FieldFunction<RealMatrixX>& mat);
            
            virtual ~InertiaMatrix() { }
            
            virtual void operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
            virtual void derivative (    const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealMatrixX& m) const;
            
        protected:
            
            const MAST::FieldFunction<RealMatrixX>& _material_inertia;
        };
        
        
        
        class ThermalExpansionMatrix: public MAST::FieldFunction<RealMatrixX> {
        public:
            ThermalExpansionMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                                   const MAST::FieldFunction<RealMatrixX>& mat_expansion);
            
            virtual ~ThermalExpansionMatrix() { }
            
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
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress);
            
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
            
            const MAST::FieldFunction<RealMatrixX>& _prestress;
        };
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond);
            
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
        };
        
        
        
        
        class ThermalCapacitanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond);
            
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
        };
        
    }
    
}


bool
MAST::IsotropicElementPropertyCard3D::depends_on(const MAST::FunctionBase& f) const {
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}



MAST::IsotropicElementProperty3D::StiffnessMatrix::
StiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat):
MAST::FieldFunction<RealMatrixX> ("StiffnessMatrix3D"),
_material_stiffness(mat) {
    
    _functions.insert(&mat);
}


void
MAST::IsotropicElementProperty3D::
StiffnessMatrix::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    // this only returns the material stiffness
    _material_stiffness(p, t, m);
}



void
MAST::IsotropicElementProperty3D::
StiffnessMatrix::derivative (const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    // this only returns the material stiffness
    _material_stiffness.derivative( f, p, t, m);
}




MAST::IsotropicElementProperty3D::
InertiaMatrix::InertiaMatrix(const MAST::FieldFunction<RealMatrixX>& mat):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(&mat);
}



void
MAST::IsotropicElementProperty3D::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {

    _material_inertia(p, t, m);
}




void
MAST::IsotropicElementProperty3D::
InertiaMatrix::derivative (               const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {

    _material_inertia.derivative( f, p, t, m);
}




MAST::IsotropicElementProperty3D::ThermalExpansionMatrix::
ThermalExpansionMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                       const MAST::FieldFunction<RealMatrixX>& mat_expansion):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion) {
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
}




void
MAST::IsotropicElementProperty3D::
ThermalExpansionMatrix::operator() (const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    RealMatrixX mat;
    _material_stiffness(p, t, m);
    _material_expansion(p, t, mat);
    m *= mat;
}






void
MAST::IsotropicElementProperty3D::
ThermalExpansionMatrix::derivative (   const MAST::FunctionBase& f,
                                    const libMesh::Point& p,
                                    const Real t,
                                    RealMatrixX& m) const {
    libmesh_error();
}




MAST::IsotropicElementProperty3D::PrestressAMatrix::
PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix3D"),
_prestress(prestress){
    _functions.insert(&prestress);
}




void
MAST::IsotropicElementProperty3D::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress(p, t, m);
}






void
MAST::IsotropicElementProperty3D::
PrestressAMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress.derivative( f, p, t, m);
}



MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond) {
    
    _functions.insert(&mat_cond);
}





MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
~ThermalConductanceMatrix() { }


void
MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    _mat_cond(p, t, m);
}




void
MAST::IsotropicElementProperty3D::ThermalConductanceMatrix::derivative (                                       const MAST::FunctionBase& f,
                                                                        const libMesh::Point& p,
                                                                        const Real t,
                                                                        RealMatrixX& m) const {
    _mat_cond.derivative( f, p, t, m);
}






MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cap):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap) {
    
    _functions.insert(&mat_cap);
}




MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() { }


void
MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    _mat_cap(p, t, m);
}



void
MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix::derivative (                                       const MAST::FunctionBase& f,
                                                                        const libMesh::Point& p,
                                                                        const Real t,
                                                                        RealMatrixX& m) const {
    _mat_cap.derivative( f, p, t, m);
}



//void
//MAST::IsotropicElementProperty3D::
//PrestressAMatrix::convert_to_vector(const RealMatrixX &m,
//                                                     RealVectorX &v) const {
//    libmesh_error();
//}




std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::StiffnessMatrix
    (_material->stiffness_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_A_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::StiffnessMatrix
    (_material->stiffness_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_B_matrix() const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
stiffness_D_matrix() const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
inertia_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::InertiaMatrix
    (_material->inertia_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
inertia_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::InertiaMatrix
    (_material->inertia_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalExpansionMatrix
    (_material->stiffness_matrix(3),
     _material->thermal_expansion_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_A_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalExpansionMatrix
    (_material->stiffness_matrix(3),
     _material->thermal_expansion_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_B_matrix(const MAST::ElementBase& e) const {
    
    // for 3D elements, there is no difference between the A and B matrices
    return this->thermal_expansion_A_matrix(e);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_expansion_B_matrix() const {
    
    // for 3D elements, there is no difference between the A and B matrices
    return this->thermal_expansion_A_matrix();
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
transverse_shear_stiffness_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
prestress_A_matrix( MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
prestress_B_matrix( MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalConductanceMatrix
    (_material->conductance_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_conductance_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalConductanceMatrix
    (_material->conductance_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::IsotropicElementPropertyCard3D::
thermal_capacitance_matrix() const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::IsotropicElementProperty3D::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(3));
    
    return std::unique_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}

