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
#include "tests/fluid/build_conservative_fluid_elem.h"
#include "tests/base/test_comparisons.h"
#include "examples/base/rigid_surface_motion.h"
#include "fluid/frequency_domain_linearized_conservative_fluid_elem.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/dof_map.h"


BOOST_FIXTURE_TEST_SUITE  (FreqDomainJacobianEvaluation, MAST::BuildConservativeFluidElem)

void
check_internal_force_jacobian(MAST::BuildConservativeFluidElem& v,
                              Real omega) {

    // set the frequency value
    *v._omega = omega;
    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(v._mesh->n_elem(), 1);
    libMesh::Elem& e = **v._mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(v._discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*v._fluid_sys, e, p));
    elem->freq   = v._freq_function;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    const Complex
    iota(0., 1.);
    
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res         = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = v._base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->internal_residual(true, res0, jac_x);
    
    
    for (unsigned int i=0; i<ndofs; i++) {
        
        
        // first the Jacobian due to x
        x      = ComplexVectorX::Zero(ndofs);
        // perturb the i^th element of the solution
        x(i)  += delta;
        elem->set_complex_solution(x);
        
        // get the new residual
        res.setZero();
        elem->internal_residual(false, res, dummy);
        
        // set the i^th column of the finite-differenced Jacobian
        jac_x_fd.col(i) = (res-res0)/delta;
    }
    
    // now compare the matrices
    BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Real Perturbation **");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Real Perturbation  **");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));

    if (fabs(omega) > 0.) {
        // now perturb the imaginary part of the solution
        for (unsigned int i=0; i<ndofs; i++) {
            
            
            // first the Jacobian due to x
            x      = ComplexVectorX::Zero(ndofs);
            // perturb the i^th element of the solution
            x(i)  += delta * iota;
            elem->set_complex_solution(x);
            
            // get the new residual
            res.setZero();
            elem->internal_residual(false, res, dummy);
            
            // set the i^th column of the finite-differenced Jacobian
            jac_x_fd.col(i) = (res-res0)/delta/iota;
        }
        
        // now compare the matrices
        BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Imaginary Perturbation **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
        
        BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Imaginary Perturbation  **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));
    }
}


void
check_far_field_jacobian(MAST::BuildConservativeFluidElem& v,
                         Real omega) {
    
    // set the frequency value
    *v._omega = omega;
    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(v._mesh->n_elem(), 1);
    libMesh::Elem& e = **v._mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(v._discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*v._fluid_sys, e, p));
    elem->freq   = v._freq_function;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    
    const Complex
    iota(0., 1.);
    
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res         = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    MAST::SideBCMapType bc_map;
    bc_map.insert(std::pair<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>
                  (1, v._far_field));
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = v._base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->side_external_residual(true, res0, jac_x, bc_map);
    
    
    for (unsigned int i=0; i<ndofs; i++) {
        
        
        // first the Jacobian due to x
        x      = ComplexVectorX::Zero(ndofs);
        // perturb the i^th element of the solution
        x(i)  += delta;
        elem->set_complex_solution(x);
        
        // get the new residual
        res.setZero();
        elem->side_external_residual(false, res, dummy, bc_map);
        
        // set the i^th column of the finite-differenced Jacobian
        jac_x_fd.col(i) = (res-res0)/delta;
    }
    
    // now compare the matrices
    BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Real Perturbation **");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Real Perturbation **");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));

    
    
    if (fabs(omega) > 0.) {
        for (unsigned int i=0; i<ndofs; i++) {
            
            
            // first the Jacobian due to x
            x      = ComplexVectorX::Zero(ndofs);
            // perturb the i^th element of the solution
            x(i)  += delta*iota;
            elem->set_complex_solution(x);
            
            // get the new residual
            res.setZero();
            elem->side_external_residual(false, res, dummy, bc_map);
            
            // set the i^th column of the finite-differenced Jacobian
            jac_x_fd.col(i) = (res-res0)/delta/iota;
        }
        
        // now compare the matrices
        BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Imaginary Perturbation **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
        
        BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Imaginary Perturbation **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));
    }
}


void
check_moving_wall_jacobian(MAST::BuildConservativeFluidElem& v,
                           Real omega) {
    
    // set the frequency value
    *v._omega = omega;
    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(v._mesh->n_elem(), 1);
    libMesh::Elem& e = **v._mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(v._discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*v._fluid_sys, e, p));
    elem->freq   = v._freq_function;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    const Complex
    iota(0., 1.);
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res         = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    MAST::SideBCMapType bc_map;
    bc_map.insert(std::pair<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>
                  (0, v._slip_wall));
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = v._base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->side_external_residual(true, res0, jac_x, bc_map);
    
    
    for (unsigned int i=0; i<ndofs; i++) {
        
        
        // first the Jacobian due to x
        x      = ComplexVectorX::Zero(ndofs);
        // perturb the i^th element of the solution
        x(i)  += delta;
        elem->set_complex_solution(x);
        
        // get the new residual
        res.setZero();
        elem->side_external_residual(false, res, dummy, bc_map);
        
        // set the i^th column of the finite-differenced Jacobian
        jac_x_fd.col(i) = (res-res0)/delta;
    }
    
    // now compare the matrices
    BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Real Perturbation**");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Real Perturbation**");
    BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));

    
    
    if (fabs(omega) > 0.) {
        
        for (unsigned int i=0; i<ndofs; i++) {
            
            
            // first the Jacobian due to x
            x      = ComplexVectorX::Zero(ndofs);
            // perturb the i^th element of the solution
            x(i)  += delta * iota;
            elem->set_complex_solution(x);
            
            // get the new residual
            res.setZero();
            elem->side_external_residual(false, res, dummy, bc_map);
            
            // set the i^th column of the finite-differenced Jacobian
            jac_x_fd.col(i) = (res-res0)/delta/iota;
        }
        
        // now compare the matrices
        BOOST_TEST_MESSAGE("** Checking Real Part of Jacobian with Imaginary Perturbation **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.real(),  jac_x.real(),  tol));
        
        BOOST_TEST_MESSAGE("** Checking Imaginary Part of Jacobian with Imaginary Perturbation **");
        BOOST_CHECK(MAST::compare_matrix(  jac_x_fd.imag(),  jac_x.imag(),  tol));
    }
}


BOOST_AUTO_TEST_CASE   (InternalForceJacobianZeroFreq) {
    
    check_internal_force_jacobian(*this, 0.);

}



BOOST_AUTO_TEST_CASE   (InternalForceJacobianNonZeroFreq) {
    
    check_internal_force_jacobian(*this, 200.);
    
}




BOOST_AUTO_TEST_CASE   (FarFieldJacobianZeroFreq) {
    
    check_far_field_jacobian(*this, 0.);
}



BOOST_AUTO_TEST_CASE   (FarFieldJacobianNonZeroFreq) {
    
    check_far_field_jacobian(*this, 200.);
}



BOOST_AUTO_TEST_CASE   (MovingWallJacobianZeroFreq) {
    
    check_moving_wall_jacobian(*this, 0.);
}


BOOST_AUTO_TEST_CASE   (MovingWallJacobianNonZeroFreq) {
    
    check_moving_wall_jacobian(*this, 200.);
}



BOOST_AUTO_TEST_CASE   (InternalForceFreqSens) {
    
    // set the frequency value to some arbitrary value
    *_omega = 100.;
    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(_mesh->n_elem(), 1);
    libMesh::Elem& e = **_mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(_discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*_fluid_sys, e, p));
    elem->freq                = _freq_function;
    elem->sensitivity_param   = _omega;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    const Complex
    iota(0., 1.);
    
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res_sens    = ComplexVectorX::Zero(ndofs),
    res_sens_fd = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = _base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    elem->set_complex_solution(x, true);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->internal_residual(true, res0, jac_x);
    
    // calculate the sensitivity of the residual
    elem->internal_residual_sensitivity(false, res_sens, dummy);
    
    // now perturb the frequency parameter
    *_omega = ((1. + delta) * (*_omega)());
    
    elem->internal_residual(false, res_sens_fd, dummy);
    res_sens_fd -= res0;
    res_sens_fd /= delta;

    // now compare the two sensitivity vectors
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
}

BOOST_AUTO_TEST_CASE   (FarFieldResidualFreqSens) {
    
    // set the frequency value to some arbitrary value
    *_omega = 100.;
    
    
    MAST::SideBCMapType bc_map;
    bc_map.insert(std::pair<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>
                  (0, _far_field));
    
    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(_mesh->n_elem(), 1);
    libMesh::Elem& e = **_mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(_discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*_fluid_sys, e, p));
    elem->freq                = _freq_function;
    elem->sensitivity_param   = _omega;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    const Complex
    iota(0., 1.);
    
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res_sens    = ComplexVectorX::Zero(ndofs),
    res_sens_fd = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = _base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    elem->set_complex_solution(x, true);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->side_external_residual(true, res0, jac_x, bc_map);
    
    // calculate the sensitivity of the residual
    elem->side_external_residual_sensitivity(false, res_sens, dummy, bc_map);
    
    // now perturb the frequency parameter
    *_omega = ((1. + delta) * (*_omega)());
    
    elem->side_external_residual(false, res_sens_fd, dummy, bc_map);
    res_sens_fd -= res0;
    res_sens_fd /= delta;

    // now compare the two sensitivity vectors
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
}


BOOST_AUTO_TEST_CASE   (MovingWallResidualFreqSens) {
    
    // set the frequency value to some arbitrary value
    *_omega = 100.;
    
    
    MAST::SideBCMapType bc_map;
    bc_map.insert(std::pair<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>
                  (0, _slip_wall));

    
    // make sure there is only one element in the mesh.
    libmesh_assert_equal_to(_mesh->n_elem(), 1);
    libMesh::Elem& e = **_mesh->local_elements_begin();
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(_discipline)->flight_condition();
    
    std::unique_ptr<MAST::FrequencyDomainLinearizedConservativeFluidElem>
    elem(new MAST::FrequencyDomainLinearizedConservativeFluidElem(*_fluid_sys, e, p));
    elem->freq                = _freq_function;
    elem->sensitivity_param   = _omega;
    
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;
    
    const Complex
    iota(0., 1.);
    
    
    // number of dofs in this element
    const unsigned int ndofs = 16;
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x_base      =RealVectorX::Zero(ndofs);
    
    ComplexVectorX
    x           = ComplexVectorX::Zero(ndofs),
    res0        = ComplexVectorX::Zero(ndofs),
    res_sens    = ComplexVectorX::Zero(ndofs),
    res_sens_fd = ComplexVectorX::Zero(ndofs);
    
    ComplexMatrixX
    jac_x       = ComplexMatrixX::Zero(ndofs, ndofs),
    jac_x_fd    = ComplexMatrixX::Zero(ndofs, ndofs),
    dummy;
    
    
    
    // set velocity to be zero
    elem->set_velocity(x_base);
    
    // initialize the base solution vector. The 2D elem has 4 variables
    for (unsigned int i=0; i<4; i++) {
        for (unsigned int j=0; j<4; j++) {
            x_base(i*4+j) = _base_sol(i);
        }
    }
    elem->set_solution(x_base);
    
    
    // both x and x0 are zero at this point.
    // tell the element about the solution and velocity
    elem->set_complex_solution(x);
    elem->set_complex_solution(x, true);
    
    // get the base residual vector and the Jacobians for numerical comparisons
    // later.
    elem->side_external_residual(true, res0, jac_x, bc_map);
    
    // calculate the sensitivity of the residual
    elem->side_external_residual_sensitivity(false, res_sens, dummy, bc_map);
    
    // now perturb the frequency parameter
    *_omega = ((1. + delta) * (*_omega)());
    
    elem->side_external_residual(false, res_sens_fd, dummy, bc_map);
    res_sens_fd -= res0;
    res_sens_fd /= delta;
    
    // now compare the two sensitivity vectors
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
    
    BOOST_TEST_MESSAGE("** Checking Real Part of Internal Residual Sensitivity **");
    BOOST_CHECK(MAST::compare_vector(  res_sens_fd.real(),  res_sens.real(),  tol));
}


BOOST_AUTO_TEST_SUITE_END()



