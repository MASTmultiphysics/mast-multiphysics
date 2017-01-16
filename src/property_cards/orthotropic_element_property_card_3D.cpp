/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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
            StiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                            const MAST::CoordinateBase& orient);
            
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
            
            const MAST::CoordinateBase& _orient;
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
                                   const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                                   const MAST::CoordinateBase& orient);
            
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
            const MAST::CoordinateBase& _orient;
        };
        
        
        
        
        class PrestressAMatrix:
        public MAST::FieldFunction<RealMatrixX> {
        public:
            PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                             const MAST::CoordinateBase& orient);
            
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
            const MAST::CoordinateBase& _orient;
        };
        
        
        
        class ThermalConductanceMatrix:
        public MAST::FieldFunction<RealMatrixX> {
            
        public:
            
            ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                                     const MAST::CoordinateBase& orient);
            
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
            
            const MAST::CoordinateBase& _orient;
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
MAST::OrthotropicElementPropertyCard3D::depends_on(const MAST::FunctionBase& f) const {
    return _material->depends_on(f) ||            // check if the material property depends on the function
    MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
}



MAST::OrthotropicProperty3D::StiffnessMatrix::
StiffnessMatrix(const MAST::FieldFunction<RealMatrixX>& mat,
                const MAST::CoordinateBase& orient):
MAST::FieldFunction<RealMatrixX> ("StiffnessMatrix3D"),
_material_stiffness(mat),
_orient(orient) {
    
    _functions.insert(&mat);
    _functions.insert(&orient);
}


void
MAST::OrthotropicProperty3D::
StiffnessMatrix::operator() (const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {
    RealMatrixX
    A    = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6);
    
    _material_stiffness (p, t, m);
    _orient             (p, t, A);
    _orient.stress_strain_transformation_matrix(A.transpose(), Tinv);
    
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
StiffnessMatrix::derivative (const MAST::FunctionBase& f,
                             const libMesh::Point& p,
                             const Real t,
                             RealMatrixX& m) const {

    RealMatrixX
    dm   = RealMatrixX::Zero(6, 6),
    A    = RealMatrixX::Zero(3, 3),
    dA   = RealMatrixX::Zero(3, 3),
    Tinv = RealMatrixX::Zero(6, 6),
    dTinv= RealMatrixX::Zero(6, 6);
    
    _material_stiffness (p, t, m);
    _material_stiffness.derivative( f, p, t, dm);
    _orient             (p, t, A);
    _orient.stress_strain_transformation_matrix(A.transpose(), Tinv);
    _orient.derivative    (f, p, t, dA);
    _orient.stress_strain_transformation_matrix_sens(A.transpose(),
                                                     dA.transpose(),
                                                     dTinv);
    
    m =
    dTinv *  m *  Tinv.transpose() +
    Tinv  * dm *  Tinv.transpose() +
    Tinv  *  m * dTinv.transpose();
}




MAST::OrthotropicProperty3D::
InertiaMatrix::InertiaMatrix(const MAST::FieldFunction<RealMatrixX>& mat):
MAST::FieldFunction<RealMatrixX>("InertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(&mat);
}



void
MAST::OrthotropicProperty3D::
InertiaMatrix::operator() (const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    // this only returns the material inertia
    _material_inertia(p, t, m);
}




void
MAST::OrthotropicProperty3D::
InertiaMatrix::derivative (               const MAST::FunctionBase& f,
                           const libMesh::Point& p,
                           const Real t,
                           RealMatrixX& m) const {
    
    _material_inertia.derivative( f, p, t, m);
}




MAST::OrthotropicProperty3D::ThermalExpansionMatrix::
ThermalExpansionMatrix(const MAST::FieldFunction<RealMatrixX>& mat_stiff,
                       const MAST::FieldFunction<RealMatrixX>& mat_expansion,
                       const MAST::CoordinateBase& orient):
MAST::FieldFunction<RealMatrixX>("ThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_orient(orient) {
    
    _functions.insert(&mat_stiff);
    _functions.insert(&mat_expansion);
    _functions.insert(&orient);
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
    
    _material_stiffness  (p, t, m);
    _material_expansion  (p, t, mat);
    _orient              (p, t, A);
    _orient.stress_strain_transformation_matrix(A.transpose(), Tinv);
    
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
ThermalExpansionMatrix::derivative (   const MAST::FunctionBase& f,
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
    
    _material_stiffness  (p, t, m);
    _material_stiffness.derivative( f, p, t, dm);
    _material_expansion (p, t, mat);
    _material_expansion.derivative( f, p, t, dmat);
    _orient              (p, t, A);
    _orient.stress_strain_transformation_matrix(A.transpose(), Tinv);
    _orient.derivative    (f, p, t, dA);
    _orient.stress_strain_transformation_matrix_sens(A.transpose(),
                                                      dA.transpose(),
                                                      dTinv);
    
    m =
    dTinv *  m *  mat +
    Tinv  * dm *  mat +
    Tinv  *  m * dmat;
}




MAST::OrthotropicProperty3D::PrestressAMatrix::
PrestressAMatrix(const MAST::FieldFunction<RealMatrixX>& prestress,
                 const MAST::CoordinateBase& orient):
MAST::FieldFunction<RealMatrixX>("PrestressAMatrix3D"),
_prestress(prestress),
_orient(orient) {
    _functions.insert(&prestress);
    _functions.insert(&orient);
}




void
MAST::OrthotropicProperty3D::
PrestressAMatrix::operator() (const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress (p, t, m);
    libmesh_error();
}






void
MAST::OrthotropicProperty3D::
PrestressAMatrix::derivative ( const MAST::FunctionBase& f,
                              const libMesh::Point& p,
                              const Real t,
                              RealMatrixX& m) const {
    _prestress.derivative( f, p, t, m);
}



MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
ThermalConductanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cond,
                         const MAST::CoordinateBase& orient):
MAST::FieldFunction<RealMatrixX>("ThermalConductanceMatrix"),
_mat_cond(mat_cond),
_orient(orient) {
    
    _functions.insert(&mat_cond);
    _functions.insert(&orient);
}



MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
~ThermalConductanceMatrix() { }


void
MAST::OrthotropicProperty3D::ThermalConductanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    RealMatrixX
    A    = RealMatrixX::Zero(3, 3);
    
    _mat_cond (p, t, m);
    _orient   (p, t, A);

    m = A.transpose() * m * A;
}




void
MAST::OrthotropicProperty3D::ThermalConductanceMatrix::derivative (                                         const MAST::FunctionBase& f,
                                                                          const libMesh::Point& p,
                                                                          const Real t,
                                                                          RealMatrixX& m) const {
    RealMatrixX
    dm    = RealMatrixX::Zero(3, 3),
    A     = RealMatrixX::Zero(3, 3),
    dA    = RealMatrixX::Zero(3, 3);
    
    _mat_cond    (p, t, m);
    _orient      (p, t, A);
    _mat_cond.derivative( f, p, t, dm);
    _orient.derivative  ( f, p, t, dA);
    
    m =
    dA.transpose() *  m * A +
    A.transpose()  * dm * A +
    A.transpose()  *  m * dA;
}






MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
ThermalCapacitanceMatrix(const MAST::FieldFunction<RealMatrixX>& mat_cap):
MAST::FieldFunction<RealMatrixX>("ThermalCapacitanceMatrix"),
_mat_cap(mat_cap) {
    
    _functions.insert(&mat_cap);
}



MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
~ThermalCapacitanceMatrix() { }


void
MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::
operator() (const libMesh::Point& p,
            const Real t,
            RealMatrixX& m) const {
    
    _mat_cap(p, t, m);
}



void
MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix::derivative (                                         const MAST::FunctionBase& f,
                                                                          const libMesh::Point& p,
                                                                          const Real t,
                                                                          RealMatrixX& m) const {
    _mat_cap.derivative( f, p, t, m);
}



MAST::OrthotropicElementPropertyCard3D::
~OrthotropicElementPropertyCard3D() { }



void
MAST::OrthotropicElementPropertyCard3D::
set_orientation(const MAST::CoordinateBase& orient) {
    
    _orient = &orient;
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::StiffnessMatrix
    (_material->stiffness_matrix(3), *_orient);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_B_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}




std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
stiffness_D_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
damping_matrix(const MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
inertia_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::InertiaMatrix
    (_material->inertia_matrix(3));
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_expansion_A_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalExpansionMatrix
    (_material->stiffness_matrix(3),
     _material->thermal_expansion_matrix(3),
     *_orient);
    
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
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
prestress_A_matrix( MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
prestress_B_matrix( MAST::ElementBase& e) const {
    
    libmesh_assert(false);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(nullptr);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_conductance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalConductanceMatrix
    (_material->conductance_matrix(3), *_orient);
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}


std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
MAST::OrthotropicElementPropertyCard3D::
thermal_capacitance_matrix(const MAST::ElementBase& e) const {
    
    MAST::FieldFunction<RealMatrixX>* rval =
    new MAST::OrthotropicProperty3D::ThermalCapacitanceMatrix
    (_material->capacitance_matrix(3));
    
    return std::auto_ptr<MAST::FieldFunction<RealMatrixX> >(rval);
}



