/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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


// BOOST includes
#include <boost/test/unit_test.hpp>


// MAST includes
#include "tests/structural/build_structural_elem_1D.h"
#include "tests/structural/build_structural_elem_2D.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/dof_map.h"


extern
void
set_deformation(const unsigned int dim,
                const unsigned int case_num,
                const libMesh::ElemType e,
                RealVectorX& vec);



template <typename ValType>
void check_internal_force_jacobian (ValType& v,
                                    const RealVectorX& sol) {

    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;

    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(v._mesh->local_elements_begin());
    
    // now create the structural element
    std::unique_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*v._structural_sys,
                                     elem,
                                     *v._p_card).release());
    
    
    // number of dofs in this element
    const libMesh::DofMap& dofmap = v._sys->get_dof_map();
    std::vector<unsigned int> dof_ids;
    dofmap.dof_indices(&elem, dof_ids);
    
    const unsigned int ndofs = (unsigned int)dof_ids.size();

    
    // make sure that the number of dofs are appropriately set
    libmesh_assert_equal_to(sol.size(), ndofs);
    
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x0          = RealVectorX::Zero(ndofs),
    x           = RealVectorX::Zero(ndofs),
    res0        = RealVectorX::Zero(ndofs),
    res         = RealVectorX::Zero(ndofs);
    
    RealMatrixX
    jac_x       = RealMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = RealMatrixX::Zero(ndofs, ndofs),
    dummy;
    

    // set the solution about which to evaluate the nonlinearity
    x0          = sol;
    x           = sol;
    
    // tell the element about the solution and velocity
    e->set_solution(x);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    e->internal_residual(true, res0, jac_x);
    
    
    for (unsigned int i=0; i<ndofs; i++) {
        
        
        // first the Jacobian due to x
        x      = x0;
        // perturb the i^th element of the solution
        x(i)  += delta;
        e->set_solution(x);
        
        // get the new residual
        res.setZero();
        e->internal_residual(false, res, dummy);
        
        // set the i^th column of the finite-differenced Jacobian
        jac_x_fd.col(i) = (res-res0)/delta;
    }
    
    // now compare the matrices
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd,  jac_x,  tol));
}



BOOST_FIXTURE_TEST_SUITE  (Structural1DJacobianEvaluation, MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear1DIndependentOffset) {

    this->init(false, false);
    
    RealVectorX v;
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear1DIndependentOffset) {
    
    this->init(false, true);
    
    RealVectorX v;
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear1DWithConstantOffset) {
    
    this->init(false, false);

    RealVectorX v;
    
    // set a finite value of the offset
    (*_hy_off)()  = 0.5*(*_thy)();
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear1DWithConstantOffset) {
    
    this->init(false, true);
    
    RealVectorX v;
    
    // set a finite value of the offset
    (*_hy_off)()  = 0.5*(*_thy)();
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear1DIndependentOffsetDependentOffset) {
    
    this->init(true, false);
    
    RealVectorX v;
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear1DIndependentOffsetDependentOffset) {
    
    this->init(true, true);
    
    RealVectorX v;
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}




BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DJacobianEvaluation,
                           MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DIndependentOffsetQUAD4) {
    
    this->init(false, false, libMesh::QUAD4);

    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DIndependentOffsetQUAD4) {
    
    this->init(false, true, libMesh::QUAD4);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DWithConstantOffsetQUAD4) {

    this->init(false, false, libMesh::QUAD4);
    
    // set a finite value of the offset
    (*_hzoff)()  = 0.5*(*_thz)();

    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DWithConstantOffsetQUAD4) {
    
    this->init(false, true, libMesh::QUAD4);
    
    // set a finite value of the offset
    (*_hzoff)()  = 0.5*(*_thz)();
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DDependentOffsetQUAD4) {
    
    this->init(true, false, libMesh::QUAD4);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DDependentOffsetQUAD4) {
    
    this->init(true, true, libMesh::QUAD4);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DIndependentOffsetTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DIndependentOffsetTRI3) {
    
    this->init(false, true, libMesh::TRI3);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DWithConstantOffsetTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    
    // set a finite value of the offset
    (*_hzoff)()  = 0.5*(*_thz)();
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DWithConstantOffsetTRI3) {
    
    this->init(false, true, libMesh::TRI3);
    
    // set a finite value of the offset
    (*_hzoff)()  = 0.5*(*_thz)();
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinear2DDependentOffsetTRI3) {
    
    this->init(true, false, libMesh::TRI3);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinear2DDependentOffsetTRI3) {
    
    this->init(true, true, libMesh::TRI3);
    
    RealVectorX v;
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_SUITE_END()



