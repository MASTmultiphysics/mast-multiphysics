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
#include "examples/structural/plate_modal_analysis/plate_modal_analysis.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "base/nonlinear_system.h"


BOOST_FIXTURE_TEST_SUITE  (StructuralPlateModalAnalysis,
                           MAST::PlateModalAnalysis)

BOOST_AUTO_TEST_CASE   (PlateModalSolution) {
    
    // initialize plate object
    this->init(libMesh::QUAD4, false);
    
    const Real
    tol      = 1.e-2;
    
    std::vector<Real>
    eig;
    
    this->solve(false, &eig);
    
    // check the solution
    // iterate over each node, and compare the nodal solution with the
    // expected anlaytical value
    Real
    th          = (*_th)(),
    rho         = (*_rho)(),
    Eval        = (*_E)(),
    nu          = (*_nu)(),
    D           = pow(th,3)/12.*Eval/(1-nu*nu),
    pi          = acos(-1.),
    analytical  = 0.;
    
    
    // analytical solution to the natural frequency of simply supported problem
    // is
    //   lambda = omega^2 = (n pi/L)^4 EI/(rho A)
    //
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    
    std::set<Real> factors_set;
    
    for (unsigned int n=1; n<4; n++) {
        for (unsigned int m=1; m<4; m++) {
            factors_set.insert((pow( n/_length,4) + pow( m/_width,4) +
                                2.* pow( n/_length,2) * pow( m/_width,2)));
        }
    }
    
    // convert this to a vector of sorted factors
    std::vector<Real> factors;
    std::set<Real>::const_iterator
    it  = factors_set.begin(),
    end = factors_set.end();
    
    for ( ; it != end; it++) {
        factors.push_back(*it);
    }
    
    // make sure that the size of factors is greater thn nconv
    libmesh_assert_greater(factors.size(), nconv);

    BOOST_TEST_MESSAGE("  ** Plate Bending Eigenvalues **");
    for (unsigned int i=0; i<nconv; i++) {
        
        analytical   = D*pow(pi,4)/rho/th * factors[i];
        
        // compare the eigenvalues
        
        BOOST_CHECK(MAST::compare_value(analytical, eig[i], tol));
    }
}



BOOST_AUTO_TEST_CASE    (PlateModalSolutionSensitivity) {
    
    // initialize plate object
    this->init(libMesh::QUAD4, false);

    const Real
    delta    = 1.e-5,
    tol      = 1.e-3;
    
    std::vector<Real>
    eig_vec,
    deig_vec;
    
    this->solve(false, &eig_vec);
    
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    
    
    // verify the sensitivity solution of this system
    RealVectorX
    eig     = RealVectorX::Zero(nconv),
    deig    = RealVectorX::Zero(nconv),
    deig_fd = RealVectorX::Zero(nconv);
    
    for (unsigned int i=0; i<nconv; i++) eig(i) = eig_vec[i];
    
    
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
        this->solve(false, &eig_vec);
        std::vector<Real> deig_vec(nconv);
        this->sensitivity_solve(f, deig_vec);
        for (unsigned int i=0; i<nconv; i++) deig(i) = deig_vec[i];
        
        // now calculate the finite difference sensitivity
        
        // identify the perturbation in the parameter
        p0           = f();
        (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
        f()         += dp;
        
        // solve at the perturbed parameter value
        this->solve(false, &eig_vec);
        for (unsigned int i=0; i<nconv; i++) deig_fd(i) = eig_vec[i];
        
        deig_fd -= eig;
        deig_fd /= dp;
        
        // reset the parameter value
        f()        = p0;
        
        // now compare the eigenvalue sensitivity
        BOOST_TEST_MESSAGE("  ** dlambda/dp (total) wrt : " << f.name() << " **");
        BOOST_CHECK(MAST::compare_vector( deig_fd,  deig, tol));
    }
}

BOOST_AUTO_TEST_SUITE_END()



