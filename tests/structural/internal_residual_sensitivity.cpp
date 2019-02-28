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

extern void
set_deformation(const unsigned int dim,
                const unsigned int case_num,
                const libMesh::ElemType e_type,
                RealVectorX& vec);


template <typename ValType>
void check_internal_force_and_jacobian_sensitivity (ValType& v,
                                                    const RealVectorX& x) {
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

    // make sure that the input dof vector is properly sized
    libmesh_assert(ndofs == x.size());
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    res0        = RealVectorX::Zero(ndofs),
    dresdp      = RealVectorX::Zero(ndofs),
    dresdp_fd   = RealVectorX::Zero(ndofs);
    
    RealMatrixX
    jac0       = RealMatrixX::Zero(ndofs, ndofs),
    djacdp     = RealMatrixX::Zero(ndofs, ndofs),
    djacdp_fd  = RealMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    Real
    p0      = 0.,
    dp      = 0.;
    
    
    // tell the element about the solution
    e->set_solution(x);
    e->internal_residual(true, res0, jac0);
    
    
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {

        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        // set the sensitivity of solution to be zero
        e->sensitivity_param  = &f;

        // get the base residual vector and the Jacobians for numerical comparisons
        // later.
        dresdp.setZero();
        djacdp.setZero();
        e->internal_residual_sensitivity(true, dresdp, djacdp);
        
        // reset the sensitivity parameter
        e->sensitivity_param  = nullptr;
        
        // now calculate the finite difference sensitivity
        
        // identify the perturbation in the parameter
        p0           = f();
        (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        dresdp_fd.setZero();
        djacdp_fd.setZero();
        e->internal_residual(true, dresdp_fd, djacdp_fd);
        
        // reset the parameter value
        f()        = p0;
        
        // calculate the finite-difference quantities
        dresdp_fd -= res0;
        dresdp_fd /= dp;
        djacdp_fd -= jac0;
        djacdp_fd /= dp;
        
        // now compare the matrices
        BOOST_TEST_MESSAGE("  ** dres/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(  dresdp_fd,  dresdp,    tol));
        BOOST_TEST_MESSAGE("  ** djac/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_matrix(   djacdp_fd,  djacdp,   tol));
    }
}



BOOST_FIXTURE_TEST_SUITE  (Structural1DInternalForceSensitivity,
                           MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity1DIndependentOffset) {
    
    this->init(false, false);
    
    RealVectorX v;

    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);

    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity1DIndependentOffset) {
    
    this->init(false, true);
    
    RealVectorX v;
    
    // large deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity1DDependentOffset) {
    
    this->init(true, false);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity1DDependentOffset) {
    
    this->init(true, true);
    
    RealVectorX v;
    
    // large deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DInternalForceSensitivity,
                           MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity2DIndependentOffsetQUAD4) {

    this->init(false, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);

    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);

    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}




BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity2DIndependentOffsetQUAD4) {
    
    this->init(false, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity2DIndependentOffsetTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity2DIndependentOffsetTRI3) {
    
    this->init(false, true, libMesh::TRI3);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}




BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity2DDependentOffsetQUAD4) {
    
    this->init(true, false, libMesh::QUAD4);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity2DDependentOffsetQUAD4) {
    
    this->init(true, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianLinearSensitivity2DDependentOffsetTRI3) {
    
    this->init(true, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}




BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonlinearSensitivity2DDependentOffsetTRI3) {
    
    this->init(true, true, libMesh::TRI3);
    
    RealVectorX v;
    
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
}



BOOST_AUTO_TEST_SUITE_END()

