/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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


// libMesh includes
#include "libmesh/dof_map.h"

extern void
set_deformation(const unsigned int dim,
                RealVectorX& vec);

extern const Real
delta,
tol;


template <typename ValType>
void check_internal_force_and_jacobian_sensitivity (ValType& v) {
    
    // get reference to the element in this mesh
    const libMesh::Elem& elem = **(v._mesh->local_elements_begin());
    
    // now create the structural element
    std::auto_ptr<MAST::StructuralElementBase>
    e(MAST::build_structural_element(*v._structural_sys,
                                     elem,
                                     *v._p_card, false).release());
    
    
    // number of dofs in this element
    const libMesh::DofMap& dofmap = v._sys->get_dof_map();
    std::vector<unsigned int> dof_ids;
    dofmap.dof_indices(&elem, dof_ids);
    
    const unsigned int ndofs = (unsigned int)dof_ids.size();
    
    // now get the residual and Jacobian evaluations
    RealVectorX
    x           = RealVectorX::Zero(ndofs),
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
    set_deformation(elem.dim(), x);
    e->set_solution(x);
    e->internal_residual(true, res0, jac0, false);
    
    
    for (unsigned int i=0; i<v._params_for_sensitivity.size(); i++) {

        MAST::Parameter& f = *v._params_for_sensitivity[i];
        
        // set the sensitivity of solution to be zero
        e->sensitivity_param  = &f;

        // get the base residual vector and the Jacobians for numerical comparisons
        // later.
        dresdp.setZero();
        djacdp.setZero();
        e->internal_residual_sensitivity(true, dresdp, djacdp, false);
        
        // reset the sensitivity parameter
        e->sensitivity_param  = NULL;
        
        // now calculate the finite difference sensitivity
        
        // identify the perturbation in the parameter
        p0           = f();
        (p0 > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        dresdp_fd.setZero();
        djacdp_fd.setZero();
        e->internal_residual(true, dresdp_fd, djacdp_fd, false);
        
        // reset the parameter value
        f()        = p0;
        
        // calculate the finite-difference quantities
        dresdp_fd -= res0;
        dresdp_fd /= dp;
        djacdp_fd -= jac0;
        djacdp_fd /= dp;
        
        // now compare the matrices
        BOOST_TEST_MESSAGE("  ** dres/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(   dresdp,    dresdp_fd, tol));
        BOOST_TEST_MESSAGE("  ** djac/dp (partial) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_matrix(   djacdp,    djacdp_fd, tol));
    }
}



BOOST_FIXTURE_TEST_SUITE  (Structural1DInternalForceSensitivity,
                           MAST::BuildStructural1DElem)

BOOST_AUTO_TEST_CASE   (InternalForceSensitivity1D) {
    
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural1DElem>(*this);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE  (Structural2DInternalForceSensitivity,
                           MAST::BuildStructural2DElem)

BOOST_AUTO_TEST_CASE   (InternalForceSensitivity2D) {
    
    check_internal_force_and_jacobian_sensitivity<MAST::BuildStructural2DElem>(*this);
}

BOOST_AUTO_TEST_SUITE_END()

