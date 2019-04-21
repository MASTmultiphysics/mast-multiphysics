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
#include "examples/structural/plate_thermally_stressed_modal_analysis/plate_thermally_stressed_modal_analysis.h"
#include "tests/base/test_comparisons.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_discipline.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "elasticity/structural_element_base.h"
#include "base/nonlinear_system.h"


BOOST_FIXTURE_TEST_SUITE  (StructuralPlateThermallyStressedModalAnalysis,
                           MAST::PlateThermallyStressedModalAnalysis)



BOOST_AUTO_TEST_CASE    (PlateThermallyStressedModalSolutionSensitivity) {
    
    // initialize plate object. Initialize with nonlinear strain
    this->init(libMesh::QUAD4, true);
    
    const Real
    delta    = 1.e-6,
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
        for (unsigned int i=0; i<deig.size(); i++)
            std::cout << deig_fd(i) << "   " << deig(i) << std::endl;
        BOOST_CHECK(MAST::compare_vector( deig_fd,  deig, tol));
    }
}

BOOST_AUTO_TEST_SUITE_END()



