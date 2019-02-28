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
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "elasticity/piston_theory_boundary_condition.h"
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
void
check_thermal_residual_force_jacobian (ValType& v,
                                       const RealVectorX& sol) {

    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;

    // tell the discipline about the section property and the piston theory
    // boundary condition
    v._discipline->add_volume_load(0, *v._thermal_load);
    
    
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
    
    // make sure that the solution vector has been appropriately dimensioned
    libmesh_assert_equal_to(sol.size(), ndofs);
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x0          = RealVectorX::Zero(ndofs),
    xdot0       = RealVectorX::Zero(ndofs),
    x           = RealVectorX::Zero(ndofs),
    xdot        = RealVectorX::Zero(ndofs),
    res0        = RealVectorX::Zero(ndofs),
    res         = RealVectorX::Zero(ndofs);
    
    RealMatrixX
    jac_x       = RealMatrixX::Zero(ndofs, ndofs),
    jac_xdot    = RealMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = RealMatrixX::Zero(ndofs, ndofs),
    jac_xdot_fd = RealMatrixX::Zero(ndofs, ndofs),
    dummy;
    

    // set the solution
    x0          = sol;
    x           = sol;
    
    
    // tell the element about the solution and velocity
    e->set_solution(x);
    e->set_velocity(xdot);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    e->volume_external_residual(true,
                                res0,
                                jac_xdot,
                                jac_x,
                                v._discipline->volume_loads());
    
    for (unsigned int i=0; i<ndofs; i++) {
        
        // first the Jacobian due to x
        x      = x0;
        xdot   = xdot0;
        // perturb the i^th element of the solution
        x(i)  += delta;
        e->set_solution(x);
        e->set_velocity(xdot);
        
        // get the new residual
        res.setZero();
        e->volume_external_residual(false,
                                    res,
                                    dummy,
                                    dummy,
                                    v._discipline->volume_loads());
        
        // set the i^th column of the finite-differenced Jacobian
        jac_x_fd.col(i) = (res-res0)/delta;
        
        
        
        // do the same for the Jacobian due to x_dot
        x        = x0;
        xdot     = xdot0;
        // perturb the i^th element of the velocity
        xdot(i) += delta;
        e->set_solution(x);
        e->set_velocity(xdot);
        
        // get the new residual
        res.setZero();
        e->volume_external_residual(false,
                                    res,
                                    dummy,
                                    dummy,
                                    v._discipline->volume_loads());
        
        // set the i^th column of the finite-differenced Jacobian
        jac_xdot_fd.col(i) = (res-res0)/delta;
    }
    
    
    // now compare the matrices
    BOOST_CHECK(MAST::compare_matrix( jac_x_fd,   jac_x,    tol));
    BOOST_CHECK(MAST::compare_matrix(jac_xdot_fd, jac_xdot, tol));
    
    // now clear the laod
    v._discipline->clear_volume_load(0, *v._thermal_load);
}




template <typename ValType>
void check_thermal_force_and_jacobian_sensitivity (ValType& v,
                                                   const RealVectorX& x) {
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;

    // tell the discipline about the section property and the piston theory
    // boundary condition
    v._discipline->add_volume_load(0, *v._thermal_load);
    
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
    libmesh_assert_equal_to(x.size(), ndofs);
    
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
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    e->volume_external_residual(true,
                                res0,
                                dummy,
                                jac0,
                                v._discipline->volume_loads());

    
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {
        
        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        // set the sensitivity of solution to be zero
        e->sensitivity_param  = &f;
        
        // get the base residual vector and the Jacobians for numerical comparisons
        // later.
        dresdp.setZero();
        djacdp.setZero();
        e->volume_external_residual_sensitivity(true,
                                                dresdp,
                                                dummy,
                                                djacdp,
                                                v._discipline->volume_loads());
        
        // reset the sensitivity parameter
        e->sensitivity_param  = nullptr;
        
        // now calculate the finite difference sensitivity
        
        // identify the perturbation in the parameter
        p0           = f();
        (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        dresdp_fd.setZero();
        djacdp_fd.setZero();
        e->volume_external_residual(true,
                                    dresdp_fd,
                                    dummy,
                                    djacdp_fd,
                                    v._discipline->volume_loads());
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

    // now clear the laod
    v._discipline->clear_volume_load(0, *v._thermal_load);
}




BOOST_FIXTURE_TEST_SUITE  (Structural1DJacobianEvaluation,
                           MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (ThermalResidualLinear1DIndependentOffset) {
    
    RealVectorX v;
    
    this->init(false, false);
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinear1DIndependentOffset) {

    RealVectorX v;
    
    this->init(false, true);
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualLinear1DDependentOffset) {
    
    RealVectorX v;
    
    this->init(true, false);
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinear1DDependentOffset) {
    
    RealVectorX v;

    this->init(true, true);
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural1DElem>(*this, v);
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural1DThermalForceSensitivity,
                           MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivity1DIndependentOffset) {

    this->init(false, false);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivity1DIndependentOffset) {
    
    this->init(false, true);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivity1DDependentOffset) {
    
    this->init(true, false);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(1, 0, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(1, 1, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(1, 2, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}




BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivity1DDependentOffset) {
    
    this->init(true, true);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(1, 3, libMesh::INVALID_ELEM, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DJacobianEvaluation,
                           MAST::BuildStructural2DElem)


BOOST_AUTO_TEST_CASE   (ThermalResidualLinearQUAD42DIndependentOffset) {
    
    RealVectorX v;
    
    this->init(false, false, libMesh::QUAD4);
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualLinearTRI32DIndependentOffset) {
    
    RealVectorX v;
    
    this->init(false, false, libMesh::TRI3);
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinearQUAD42DIndependentOffset) {
    
    RealVectorX v;
    
    this->init(false, true, libMesh::QUAD4);
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinearTRI32DIndependentOffset) {
    
    RealVectorX v;
    
    this->init(false, true, libMesh::TRI3);
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualLinearQUAD42DDependentOffset) {
    
    RealVectorX v;
    
    this->init(true, false, libMesh::QUAD4);
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualLinearTRI32DDependentOffset) {
    
    RealVectorX v;
    
    this->init(true, false, libMesh::TRI3);
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinearQUAD42DDependentOffset) {
    
    RealVectorX v;
    
    this->init(true, true, libMesh::QUAD4);
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}



BOOST_AUTO_TEST_CASE   (ThermalResidualNonlinearTRI32DDependentOffset) {
    
    RealVectorX v;
    
    this->init(true, true, libMesh::TRI3);
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_residual_force_jacobian<MAST::BuildStructural2DElem>(*this, v);
}


BOOST_AUTO_TEST_SUITE_END()




BOOST_FIXTURE_TEST_SUITE  (Structural2DThermalForceSensitivity,
                           MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivityQUAD42DIndependentOffset) {
    
    this->init(false, false, libMesh::QUAD4);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivityTRI32DIndependentOffset) {
    
    this->init(false, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}




BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivityQUAD42DIndependentOffset) {
    
    this->init(false, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivityTRI32DIndependentOffset) {
    
    this->init(false, true, libMesh::TRI3);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}




BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivityQUAD42DDependentOffset) {
    
    this->init(true, false, libMesh::QUAD4);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (ThermalForceLinearSensitivityTRI32DDependentOffset) {
    
    this->init(true, false, libMesh::TRI3);
    
    RealVectorX v;
    
    // pure axial deformation
    BOOST_TEST_MESSAGE("**** Pure Extension Deformation **");
    set_deformation(2, 0, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // pure bending deformation
    BOOST_TEST_MESSAGE("**** Pure Bending Deformation **");
    set_deformation(2, 1, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Deformation **");
    set_deformation(2, 2, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}




BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivityQUAD42DDependentOffset) {
    
    this->init(true, true, libMesh::QUAD4);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::QUAD4, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_CASE   (ThermalForceNonlinearSensitivityTRI32DDependentOffset) {
    
    this->init(true, true, libMesh::TRI3);
    
    RealVectorX v;
    
    // combination of axial and bending deformation
    BOOST_TEST_MESSAGE("**** Combined Extension-Bending Large Deformation **");
    set_deformation(2, 3, libMesh::TRI3, v);
    check_thermal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>
    (*this, v);
    
}



BOOST_AUTO_TEST_SUITE_END()
