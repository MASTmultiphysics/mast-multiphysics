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


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MAST_TESTS
#include <boost/test/unit_test.hpp>


// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"


libMesh::LibMeshInit     *_libmesh_init         = nullptr;
const Real                _tol                  = 1.e-6;

struct GlobalTestFixture {
    
    GlobalTestFixture() {
        
        // create the libMeshInit function
        _libmesh_init =
        new libMesh::LibMeshInit(boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv);
    }
    
    ~GlobalTestFixture() {
        
        delete _libmesh_init;
    }
    
};


BOOST_TEST_GLOBAL_FIXTURE( GlobalTestFixture );



// Test includes
#include "fluid/base/fluid_elem_initialization.h"
#include "base/test_comparisons.h"


BOOST_FIXTURE_TEST_SUITE(ConservativeFluidElemJacobians, BuildFluidElem)


// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
BOOST_AUTO_TEST_CASE(SlipWallJacobian) {
    
    Real
    frac  = 1.e-7,
    delta = 0.;
    
    unsigned int
    n = _sys->n_dofs();
    
    RealMatrixX
    jac0    = RealMatrixX::Zero(n, n),
    jac_fd  = RealMatrixX::Zero(n, n);
    
    RealVectorX
    v0      = RealVectorX::Zero(n),
    x       = RealVectorX::Zero(n),
    x0      = RealVectorX::Zero(n),
    f0      = RealVectorX::Zero(n),
    f       = RealVectorX::Zero(n);
    
    init_solution_for_elem(x0);
    
    // velocity is set to assuming a linear variation of the state from zero
    // over dt = 1.e-2
    v0 = x0/1.e-2;
    
    _fluid_elem->set_solution(x0);
    _fluid_elem->set_velocity(v0);
    _fluid_elem->far_field_surface_residual(true,
                                            f0,
                                            jac0,
                                            0,
                                            *_far_field_bc);
    
    for (unsigned int i=0; i<n; i++) {
        
        x = x0;
        
        if (fabs(x0(i)) > 0.) {
            delta = x0(i)*frac;
        }
        else {
            delta = frac;
        }
        
        
        x(i) += delta;
        
        _fluid_elem->set_solution(x);
        _fluid_elem->set_velocity(v0);
        
        // get the new residual
        f.setZero();
        _fluid_elem->far_field_surface_residual(false,
                                                f,
                                                jac0,
                                                0,
                                                *_far_field_bc);
        
        // set the i^th column of the finite-differenced Jacobian
        jac_fd.col(i) = (f-f0)/delta;
    }
    
    
    BOOST_CHECK(MAST::compare_matrix(jac0, jac_fd, _tol));
}

BOOST_AUTO_TEST_SUITE_END()
