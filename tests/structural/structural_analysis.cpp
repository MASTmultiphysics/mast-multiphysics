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
#include "examples/structural/bar_extension/bar_extension.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


extern const Real
delta,
tol;


BOOST_FIXTURE_TEST_SUITE  (Structural1DBarExtension,
                           MAST::BarExtension)

BOOST_AUTO_TEST_CASE   (BarExtensionSolution) {
    
    this->solve();
}


BOOST_AUTO_TEST_CASE   (BarExtensionSensitivity) {
    
    // verify the sensitivity solution of this system
    RealVectorX
    sol,
    dsol,
    dsol_fd;
    
    this->solve();
    const libMesh::NumericVector<Real>& sol_vec = this->solve();
    
    
    const unsigned int n = sol_vec.size();
    sol      =   RealVectorX::Zero(n);
    dsol     =   RealVectorX::Zero(n);
    
    // copy the solution for later use
    for (unsigned int i=0; i<n; i++)  sol(i)  = sol_vec(i);
    
    // now iterate over all the parameters and calculate the analytical
    // sensitivity and compare with the numerical sensitivity

    Real
    p0    = 0.,
    dp    = 0.;

    

    for (unsigned int i=0; i<_params_for_sensitivity.size(); i++ ) {

        MAST::Parameter& f = *this->_params_for_sensitivity[i];

        dsol     =   RealVectorX::Zero(n);
        dsol_fd  =   RealVectorX::Zero(n);

        // calculate the analytical sensitivity
        _discipline->add_parameter(f);
        const libMesh::NumericVector<Real>& dsol_vec = this->sensitivity_solve(f);
        for (unsigned int i=0; i<n; i++)  dsol(i)  = dsol_vec(i);
        _discipline->remove_parameter(f);
        
        // now calculate the finite difference sensitivity

        // identify the perturbation in the parameter
        p0           = f();
        (p0 > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        const libMesh::NumericVector<Real>& sol_vec1 = this->solve();
        for (unsigned int i=0; i<n; i++)  dsol_fd(i)  = sol_vec1(i);

        dsol_fd -= sol;
        dsol_fd /= dp;

        // reset the parameter value
        f()        = p0;

        // now compare the matrices
        BOOST_TEST_MESSAGE("  ** dX/dp (total) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector(   dsol,    dsol_fd, tol));
    }
    
}


BOOST_AUTO_TEST_SUITE_END()



