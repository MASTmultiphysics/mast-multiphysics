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
BOOST_AUTO_TEST_CASE(FarFieldJacobian) {

    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            tol;
        unsigned int    side;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->far_field_surface_residual(jac, f, j, side,
                                                       *e->_far_field_bc);
        }
    };

    Check val;
    val.jac_xdot = false;
    val.frac     = _tol*1.e-1;
    val.tol      = _tol;
    val.e        = this;

    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().n_sides(); val.side++)
        BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(SlipWallJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            tol;
        unsigned int    side;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->slip_wall_surface_residual(jac, f, j, side,
                                                       *e->_slip_wall_bc);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _tol*1.e-1;
    val.tol      = _tol;
    val.e        = this;

    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().n_sides(); val.side++)
        BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(InternalResidualJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->internal_residual(jac, f, j);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _tol*1.e-1;
    // a smaller tolerance is required for the internal resisudal since
    // an exact Jacobian is not computed for the stabilization terms
    val.tol  = 1.e-2;
    val.e        = this;
    
    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(VelocityResidualXJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            RealMatrixX jj = j; // dummy for jac wrt x_dot
            e->_fluid_elem->velocity_residual(jac, f, jj, j);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _tol*1.e-1;
    val.tol      = _tol;
    val.e        = this;

    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(VelocityResidualXdotJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            RealMatrixX jj = j; // dummy for jac wrt x
            e->_fluid_elem->velocity_residual(jac, f, j, jj);
        }
    };
    
    Check val;
    val.jac_xdot = true;
    val.frac     = _tol*1.e-1;
    val.tol      = _tol;
    val.e        = this;

    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_SUITE_END()
