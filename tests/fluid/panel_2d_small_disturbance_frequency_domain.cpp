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
#include "examples/fluid/panel_small_disturbance_frequency_domain_analysis_2D/panel_small_disturbance_frequency_domain_analysis_2d.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"
#include "tests/base/test_comparisons.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


BOOST_FIXTURE_TEST_SUITE  (PanelSmallDisturbanceFrequencyDomain2D,
                           MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis)

BOOST_AUTO_TEST_CASE   (FreqDomainSensitivityWrtOmega) {
    
    MAST::Parameter& f = *_omega;
    f = 100.;  // some arbitrary value
    
    const Real
    delta    = 1.e-5,
    tol      = 1.e-3;
    
    // verify the sensitivity solution of this system
    RealVectorX
    sol_re,
    sol_im,
    dsol_re,
    dsol_im,
    dsol_fd_re,
    dsol_fd_im;
    
    solve();
    
    // get the solutions
    std::string
    nm_re      = _sys->name() + "real_sol",
    nm_im      = _sys->name() + "imag_sol",
    nm_re_sens = _sys->name() + "real_sol_sens",
    nm_im_sens = _sys->name() + "imag_sol_sens";
    
    libMesh::NumericVector<Real>
    & sol_vec_re = _sys->get_vector(nm_re),
    & sol_vec_im = _sys->get_vector(nm_im);
    
    const unsigned int
    n_dofs     = sol_vec_re.size();
    
    sol_re       =   RealVectorX::Zero(n_dofs);
    sol_im       =   RealVectorX::Zero(n_dofs);
    dsol_re      =   RealVectorX::Zero(n_dofs);
    dsol_im      =   RealVectorX::Zero(n_dofs);
    dsol_fd_re   =   RealVectorX::Zero(n_dofs);
    dsol_fd_im   =   RealVectorX::Zero(n_dofs);
    
    // copy the solution for later use
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol_re(i)  = sol_vec_re(i);
        sol_im(i)  = sol_vec_im(i);
    }
    
    // now iterate over all the parameters and calculate the analytical
    // sensitivity and compare with the numerical sensitivity
    
    Real
    p0    = 0.,
    dp    = 0.;
    
    /////////////////////////////////////////////////////////
    // now evaluate the direct sensitivity
    /////////////////////////////////////////////////////////
    
    //for (unsigned int i=0; i<_params_for_sensitivity.size(); i++ ) {
    
    //MAST::Parameter& f = *_params_for_sensitivity[i];
    
    dsol_re      =   RealVectorX::Zero(n_dofs);
    dsol_im      =   RealVectorX::Zero(n_dofs);
    dsol_fd_re   =   RealVectorX::Zero(n_dofs);
    dsol_fd_im   =   RealVectorX::Zero(n_dofs);
    
    // calculate the analytical sensitivity
    // analysis is required at the baseline before sensitivity solution
    // and the solution has changed after the previous perturbed solution
    solve();
    sensitivity_solve(f);
    
    libMesh::NumericVector<Real>
    & dsol_vec_re = _sys->get_vector(nm_re_sens),
    & dsol_vec_im = _sys->get_vector(nm_im_sens);

    // copy the sensitivity solution
    for (unsigned int j=0; j<n_dofs; j++) {
        
        dsol_re(j)  = dsol_vec_re(j);
        dsol_im(j)  = dsol_vec_im(j);
    }

    // now calculate the finite difference sensitivity
    
    // identify the perturbation in the parameter
    p0           = f();
    (fabs(p0) > 0)?  dp=delta*p0 : dp=delta;
    f()         += dp;
    
    // solve at the perturbed parameter value
    solve();
    
    libMesh::NumericVector<Real>
    & dsol_vec_re_fd = _sys->get_vector(nm_re),
    & dsol_vec_im_fd = _sys->get_vector(nm_im);

    // copy the perturbed solution
    for (unsigned int j=0; j<n_dofs; j++) {
        
        dsol_fd_re(j)  = dsol_vec_re_fd(j);
        dsol_fd_im(j)  = dsol_vec_im_fd(j);
    }
    
    // calculate the finite difference sensitivity for solution
    dsol_fd_re -= sol_re;
    dsol_fd_im -= sol_im;
    dsol_fd_re /= dp;
    dsol_fd_im /= dp;
    
    // reset the parameter value
    f()        = p0;
    
    // now compare the solution sensitivity
    BOOST_TEST_MESSAGE("  ** dX_re/dp (total) wrt : " << f.name() << " **");
    BOOST_CHECK(MAST::compare_vector( dsol_fd_re,  dsol_re, tol));

    BOOST_TEST_MESSAGE("  ** dX_im/dp (total) wrt : " << f.name() << " **");
    BOOST_CHECK(MAST::compare_vector( dsol_fd_im,  dsol_im, tol));

    //}
    
}



BOOST_AUTO_TEST_SUITE_END()

