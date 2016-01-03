/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "examples/structural/beam_modal_analysis/beam_modal_analysis.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"

// MAST includes
#include "base/nonlinear_system.h"


BOOST_FIXTURE_TEST_SUITE  (Structural1DBeamModalAnalysis,
                           MAST::BeamModalAnalysis)

BOOST_AUTO_TEST_CASE   (BeamModalSolution) {
    

    const Real
    tol      = 1.e-2;
    
    this->solve();
    
    // check the solution
    // iterate over each node, and compare the nodal solution with the
    // expected anlaytical value
    Real
    th_y        = (*_thy)(),
    th_z        = (*_thz)(),
    Izz         = pow(th_z,1)*pow(th_y,3)/12.,
    A           = th_z*th_y,
    rho         = (*_rho)(),
    Eval        = (*_E)(),
    pi          = acos(-1.),
    analytical  = 0.,
    re          = 0.,
    im          = 0.;
    
    
    // analytical solution to the natural frequency of simply supported problem
    // is
    //   lambda = omega^2 = (n pi/L)^4 EI/(rho A)
    //
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());

    for (unsigned int i=0; i<nconv; i++) {
        
        analytical   = Eval*Izz/(rho*A) * pow((i+1)*pi/_length, 4);
        _sys->get_eigenpair(i, re, im, *_sys->solution);
        
        // compare the eigenvalues
        BOOST_CHECK(MAST::compare_value(analytical, re, tol));
    }
}


BOOST_AUTO_TEST_SUITE_END()



