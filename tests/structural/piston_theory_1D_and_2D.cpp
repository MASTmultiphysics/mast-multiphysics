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



template <typename ValType>
void
check_piston_theory_jacobian (ValType& v) {

    const Real
    delta    = 1.e-5,
    tol      = 1.e-2;

    // tell the discipline about the section property and the piston theory
    // boundary condition
    v._discipline->add_volume_load(0, *v._p_theory);
    
    
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
}


BOOST_FIXTURE_TEST_SUITE  (Structural1DJacobianEvaluation, MAST::BuildStructural1DElem)


void
check_value_sensitivity(MAST::FieldFunction<Real>& func,
                        std::vector<MAST::Parameter*>& params) {
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-3;
    
    // calculate the base pressure
    Real
    p0    = 0.,
    dp    = 0.,
    dp_fd = 0.,
    v0    = 0.,
    dv    = 0.;
    
    libMesh::Point pt;
    
    func(pt, 0., p0);
    
    // now calculate sensitivity wrt the parameters
    for (unsigned int i=0; i<params.size(); i++) {
        
        MAST::Parameter& f = *params[i];
        
        // calculate the partial derivative
        func.derivative(*params[i],
                        pt,
                        0.,
                        dp);
        
        // now perturb the parameter and calculate the finite difference
        // sensitivity of pressure.
        v0           = f();
        (fabs(v0) > 0)?  dv=delta*v0 : dv=delta;
        f()         += dv;
        
        
        func(pt, 0., dp_fd);
        
        dp_fd -= p0;
        dp_fd /= dv;
        
        // reset the parameter value
        f()        = v0;
        
        BOOST_TEST_MESSAGE("  ** dfunc/dp (total) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value( dp_fd,  dp, tol));
    }
}


BOOST_AUTO_TEST_CASE   (PistonTheoryPressure) {
    
    const Real
    tol      = 1.e-2;

    this->init(false, false);

    // calculate the surface pressure and its sensitivity with respect to the
    // relevant parameters
    std::unique_ptr<MAST::FieldFunction<Real> >
    pressure       (_p_theory->get_pressure_function(*_dwdx_f, *_dwdt_f).release()),
    dpressure_dx   (_p_theory->get_dpdx_function    (*_dwdx_f, *_dwdt_f).release()),
    dpressure_dxdot(_p_theory->get_dpdxdot_function (*_dwdx_f, *_dwdt_f).release());

    
    // first make sure that the sensitivity from pressure wrt x and xdot
    // is the same as analysis from the other two respective functions
    Real
    dpdx1     = 0.,
    dpdx2     = 0.,
    dpdxdot1  = 0.,
    dpdxdot2  = 0.;
    
    libMesh::Point pt;
    
    // first the sensitivity results
    pressure->derivative(*_dwdx, pt, 0., dpdx1);
    pressure->derivative(*_dwdt, pt, 0., dpdxdot1);

    // next, the analysis from the derivative functions
    (*dpressure_dx)   (pt, 0., dpdx2);
    (*dpressure_dxdot)(pt, 0., dpdxdot2);
    
    BOOST_TEST_MESSAGE("  ***** dpress/dx (total)  **");
    BOOST_CHECK(MAST::compare_value( dpdx1,  dpdx2, tol));

    BOOST_TEST_MESSAGE("  ***** dpress/dxdot (total)  **");
    BOOST_CHECK(MAST::compare_value( dpdxdot1,  dpdxdot2, tol));

    
    // now get sensitivity wrt all parameters
    std::vector<MAST::Parameter*> params(4);
    params[0] = _velocity;
    params[1] = _mach;
    params[2] = _rho_air;
    params[3] = _gamma_air;

    // sensitivity of these two functions is only wrt 4 params
    BOOST_TEST_MESSAGE("  ***** dpress_dx/dp (total)  **");
    check_value_sensitivity(*dpressure_dx,    params);
    
    BOOST_TEST_MESSAGE("  ***** dpress_dxdot/dp (total)  **");
    check_value_sensitivity(*dpressure_dxdot, params);
    

    // now add the two parameters to check the sensitivity of the pressure
    params.push_back(_dwdx);
    params.push_back(_dwdt);

    BOOST_TEST_MESSAGE("  ***** dpress/dp (total)  **");
    check_value_sensitivity(*pressure,        params);
    
}



BOOST_AUTO_TEST_CASE   (PistonTheory1D) {
    
    this->init(false, false); // keeping nonlinear strain off, since it does not influence piston theory
    check_piston_theory_jacobian<MAST::BuildStructural1DElem>(*this);
    
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DJacobianEvaluation, MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (PistonTheory2DQUAD4) {
    
    this->init(false, false, libMesh::QUAD4);
    check_piston_theory_jacobian<MAST::BuildStructural2DElem>(*this);
}


BOOST_AUTO_TEST_CASE   (PistonTheory2DTRI3) {
    
    this->init(false, false, libMesh::TRI3);
    check_piston_theory_jacobian<MAST::BuildStructural2DElem>(*this);
}


BOOST_AUTO_TEST_SUITE_END()


