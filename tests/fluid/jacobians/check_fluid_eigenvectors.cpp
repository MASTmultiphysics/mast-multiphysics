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
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/libmesh.h"


libMesh::LibMeshInit     *_libmesh_init         = nullptr;
const Real                _frac                 = 1.e-6;
const Real                _delta                = 1.e-6;
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


#if BOOST_VERSION > 106100
BOOST_TEST_GLOBAL_FIXTURE( GlobalTestFixture );
#else
BOOST_GLOBAL_FIXTURE( GlobalTestFixture );
#endif



// Test includes
#include "fluid/base/fluid_elem_initialization.h"
#include "base/test_comparisons.h"


BOOST_FIXTURE_TEST_SUITE(ConservativeFluidElemEigenvectors, BuildFluidElem)


// sensitivity wrt rho is all zero. Hence, the first one tested is wrt u1
BOOST_AUTO_TEST_CASE(Eigenvectors) {
    
    this->init(false);
 
    MAST::PrimitiveSolution p_sol;

    libMesh::Point
    nvec;
    
    RealVectorX
    eig_vals_from_vec =  RealVectorX::Zero(this->_dim+2),
    eig_vals          =  RealVectorX::Zero(this->_dim+2);
    
    RealMatrixX
    tmp               =  RealMatrixX::Zero(this->_dim+2, this->_dim+2),
    jac               =  RealMatrixX::Zero(this->_dim+2, this->_dim+2),
    l_eig_mat         =  RealMatrixX::Zero(this->_dim+2, this->_dim+2),
    l_eig_mat_inv_tr  =  RealMatrixX::Zero(this->_dim+2, this->_dim+2);
    
    std::vector<RealMatrixX>
    jacx(this->_dim);

    this->init_primitive_sol(p_sol);
    
    for (unsigned int i=0; i<this->_dim; i++) {
        
        jacx[i] = RealMatrixX::Zero(this->_dim+2, this->_dim+2);
        this->_fluid_elem->calculate_advection_flux_jacobian(i, p_sol, jacx[i]);
    }
    
    for (unsigned int i=0; i<this->_dim; i++) {
        
        // zero all values
        nvec.zero();
        jac.setZero();
        eig_vals.setZero();
        l_eig_mat.setZero();
        l_eig_mat_inv_tr.setZero();
        
        nvec(i) = 1;
        
        for (unsigned int j=0 ; j<this->_dim; j++) {
            jac += nvec(j) * jacx[j];
        }
        
        this->_fluid_elem->calculate_advection_left_eigenvector_and_inverse_for_normal
        (p_sol, nvec, eig_vals, l_eig_mat, l_eig_mat_inv_tr);
        
        // compare eigenvalues
        tmp = l_eig_mat.transpose() *  jac * l_eig_mat_inv_tr;
        BOOST_CHECK(MAST::compare_vector(eig_vals, tmp.diagonal(), _tol));
        
        // compare the jacobian reconstruction
        for (unsigned int j=0; j<l_eig_mat_inv_tr.cols(); j++)
            l_eig_mat_inv_tr.col(j) *= eig_vals(j);
        tmp =   l_eig_mat_inv_tr * l_eig_mat.transpose();
        BOOST_CHECK(MAST::compare_matrix(jac, tmp, _tol));
    }
}

BOOST_AUTO_TEST_SUITE_END()
