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
#include "examples/structural/beam_piston_theory_flutter/beam_piston_theory_flutter.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "base/nonlinear_system.h"


BOOST_FIXTURE_TEST_SUITE  (Structural1DBeamPistonTheoryFlutterAnalysis,
                           MAST::BeamPistonTheoryFlutterAnalysis)

/*
BOOST_AUTO_TEST_CASE   (BeamPistonTheoryFlutterSolution) {
    
    
    const Real
    tol      = 1.e-2;
    
    std::vector<Real>
    eig;
    
    this->solve(false, &eig);
    
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
    analytical  = 0.;
    
    
    // analytical solution to the natural frequency of simply supported problem
    // is
    //   lambda = omega^2 = (n pi/L)^4 EI/(rho A)
    //
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    
    for (unsigned int i=0; i<nconv; i++) {
        
        analytical   = Eval*Izz/(rho*A) * pow((i+1)*pi/_length, 4);
        
        // compare the eigenvalues
        BOOST_CHECK(MAST::compare_value(analytical, eig[i], tol));
    }
}
*/


BOOST_AUTO_TEST_CASE    (BeamPistonTheoryFlutterSolutionSensitivity) {
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-3;
    
    Real
    V0       = 0.,
    dV       = 0.,
    dV_fd    = 0.;
    
    this->init(libMesh::EDGE2, false);

    std::vector<Real>
    eig_vec,
    deig_vec;
    
    V0       = this->solve(false);
    
    // now iterate over all the parameters and calculate the analytical
    // sensitivity and compare with the numerical sensitivity
    
    Real
    p0    = 0.,
    dp    = 0.;
    
    /////////////////////////////////////////////////////////
    // now evaluate the direct sensitivity
    /////////////////////////////////////////////////////////
    
    for (unsigned int i=0; i<this->_params_for_sensitivity.size(); i++ ) {
        
        MAST::Parameter& f = *this->_params_for_sensitivity[i];
        
        // calculate the analytical sensitivity
        // analysis is required at the baseline before sensitivity solution
        // and the solution has changed after the previous perturbed solution
        this->solve(false, 1.e-4, 100);
        dV = this->sensitivity_solve(f);
        
        // now calculate the finite difference sensitivity
        
        // identify the perturbation in the parameter
        p0           = f();
        (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        // solve at the perturbed parameter value
        dV_fd        = this->solve(false);

        dV_fd       -= V0;
        dV_fd       /= dp;
        
        // reset the parameter value
        f()        = p0;
        
        // now compare the eigenvalue sensitivity
        BOOST_TEST_MESSAGE("  ** dV_F/dp (total) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_value( dV_fd,  dV, tol));
    }
}

BOOST_AUTO_TEST_SUITE_END()



