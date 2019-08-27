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
#include <boost/version.hpp>

// C++ includes
#include <memory>

// MAST includes
#include "base/mast_data_types.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/libmesh.h"


libMesh::LibMeshInit     *_libmesh_init         = nullptr;
const Real                _frac                 = 1.e-4;
const Real                _delta                = 1.e-4;
const Real                _tol                  = 1.e-5;

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


#if BOOST_VERSION > 106100
BOOST_TEST_GLOBAL_FIXTURE( GlobalTestFixture );
#else
BOOST_GLOBAL_FIXTURE( GlobalTestFixture );
#endif


// Test includes
#include "fluid/base/fluid_elem_initialization.h"
#include "base/test_comparisons.h"


BOOST_FIXTURE_TEST_SUITE(ConservativeFluidElemJacobians, BuildFluidElem)


// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
BOOST_AUTO_TEST_CASE(FarFieldJacobian) {

    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
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
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->init(false);

    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().get_reference_elem().n_sides(); val.side++)
        BOOST_CHECK(check_jacobian(val));
    
}



BOOST_AUTO_TEST_CASE(FarFieldResSens) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
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
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->init(false);
    
    Real
    rho0          = this->_flight_cond->gas_property.rho;
    
    RealVectorX
    sol           = RealVectorX::Zero(this->_sys->n_dofs()),
    res_sens      = RealVectorX::Zero(this->_sys->n_dofs()),
    res_up        = RealVectorX::Zero(this->_sys->n_dofs()),
    res_lo        = RealVectorX::Zero(this->_sys->n_dofs());
    
    RealMatrixX
    jac;
    
    MAST::Parameter p("dummy", 0.);
    
    this->init_solution_for_elem(sol);
    this->_fluid_elem->set_solution(sol);
    
    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().get_reference_elem().n_sides(); val.side++) {
        
        res_sens.setZero();
        res_up.setZero();
        res_lo.setZero();
        
        this->_flight_cond->gas_property.rho = rho0;
        this->_flight_cond->init();
        this->_fluid_elem->far_field_surface_residual_sensitivity(p, false, res_sens, jac, val.side, *this->_far_field_bc);
        
        
        
        this->_flight_cond->gas_property.rho = rho0 * (1. + _frac);
        this->_flight_cond->init();
        this->_fluid_elem->far_field_surface_residual(false, res_up, jac, val.side, *this->_far_field_bc);

        this->_flight_cond->gas_property.rho = rho0 * (1. - _frac);
        this->_flight_cond->init();
        this->_fluid_elem->far_field_surface_residual(false, res_lo, jac, val.side, *this->_far_field_bc);
        
        res_up -= res_lo;
        res_up /= (2.*_frac*rho0);
        
        BOOST_CHECK(MAST::compare_vector(res_up, res_sens, _tol));
    }
}


BOOST_AUTO_TEST_CASE(SlipWallJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
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
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->init(false);

    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().get_reference_elem().n_sides(); val.side++)
        BOOST_CHECK(check_jacobian(val));
    
}



BOOST_AUTO_TEST_CASE(NoSlipWallJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
        Real            tol;
        unsigned int    side;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->noslip_wall_surface_residual(jac, f, j, side,
                                                         *e->_noslip_wall_bc);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->_delta = 0.e-4;
    this->init(true);
    
    // check for each side.
    for (val.side=0; val.side<_fluid_elem->elem().get_reference_elem().n_sides(); val.side++)
        BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(InternalResidualJacobianInviscid) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->internal_residual(jac, f, j);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _frac;
    val.delta    = _delta;
    // a smaller tolerance is required for the internal resisudal since
    // an exact Jacobian is not computed for the stabilization terms
    val.tol      = 1.e-2;
    val.e        = this;
    this->init(false);

    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(InternalResidualJacobianViscous) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            e->_fluid_elem->internal_residual(jac, f, j);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _frac;
    val.delta    = _delta;
    // a smaller tolerance is required for the internal resisudal since
    // an exact Jacobian is not computed for the stabilization terms
    val.tol      = 1.e-2;
    val.e        = this;
    this->_delta = 1.e-4;
    this->init(true);
    
    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(VelocityResidualXJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            RealMatrixX jj = j; // dummy for jac wrt x_dot
            e->_fluid_elem->velocity_residual(jac, f, jj, j);
        }
    };
    
    Check val;
    val.jac_xdot = false;
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->init(false);

    //BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_CASE(VelocityResidualXdotJacobian) {
    
    struct Check {
        bool            jac_xdot;
        Real            frac;
        Real            delta;
        Real            tol;
        BuildFluidElem* e;
        void compute(bool jac, RealVectorX& f, RealMatrixX& j) {
            RealMatrixX jj = j; // dummy for jac wrt x
            e->_fluid_elem->velocity_residual(jac, f, j, jj);
        }
    };
    
    Check val;
    val.jac_xdot = true;
    val.frac     = _frac;
    val.delta    = _delta;
    val.tol      = _tol;
    val.e        = this;
    this->init(false);

    BOOST_CHECK(check_jacobian(val));
    
}


BOOST_AUTO_TEST_SUITE_END()
