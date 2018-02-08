/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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

// MAST includes
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "base/elem_base.h"



MAST::NonlinearImplicitAssemblyElemOperations::NonlinearImplicitAssemblyElemOperations():
MAST::AssemblyElemOperations() {
    
}


MAST::NonlinearImplicitAssemblyElemOperations::
~NonlinearImplicitAssemblyElemOperations() {
    
}



namespace MAST {
    
    bool
    is_numerical_zero(const Real v, const Real eps) {
        
        return fabs(v) <= eps;
    }
    
    
    bool
    compare(const Real v1, const Real v2, const Real tol) {
        
        const Real
        eps      = 1.0e-7;
        
        bool rval = false;
        
        // check to see if the values are both small enough
        // to be zero
        if (MAST::is_numerical_zero(v1, eps) &&
            MAST::is_numerical_zero(v2, eps))
            rval = true;
        // check to see if the absolute difference is small enough
        else if (MAST::is_numerical_zero(v1-v2, eps))
            rval = true;
        // check to see if the relative difference is small enough
        else if (fabs(v1) > 0)
            rval = fabs((v1-v2)/v1) <= tol;
        
        return rval;
    }
    
    bool
    compare_matrix(const RealMatrixX& m0, const RealMatrixX& m, const Real tol) {
        
        unsigned int
        m0_rows = (unsigned int) m0.rows(),
        m0_cols = (unsigned int) m0.cols();
        libmesh_assert_equal_to(m0_rows,  m.rows());
        libmesh_assert_equal_to(m0_cols,  m.cols());
        
        
        bool pass = true;
        for (unsigned int i=0; i<m0_rows; i++) {
            for (unsigned int j=0; j<m0_cols; j++)
                if (!MAST::compare(m0(i,j), m(i,j), tol)) {
                    libMesh::out << "Failed comparison at (i,j) = ("
                    << i << ", " << j << ") : "
                    << "expected: " << m0(i,j) << "  , "
                    << "computed: " << m(i,j) << " : "
                    << "diff: " << m0(i,j) - m(i,j) << " , "
                    << "tol: " << tol << std::endl;
                    pass = false;
                }
        }
        
        return pass;
    }
}



void
MAST::NonlinearImplicitAssemblyElemOperations::
check_element_numerical_jacobian(RealVectorX& sol) {
    RealVectorX
    dsol,
    res0,
    dres;
    
    RealMatrixX
    jac0,
    jac,
    dummy;
    
    unsigned int ndofs = (unsigned int)sol.size();
    res0.setZero(ndofs);
    dres.setZero(ndofs);
    jac0.setZero(ndofs, ndofs);
    jac.setZero(ndofs, ndofs);
    
    this->set_elem_solution(sol);
    this->elem_calculations(true, res0, jac0);
    Real delta = 1.0e-8;
    
    for (unsigned int i=0; i<sol.size(); i++) {
        dsol = sol;
        dsol(i) += delta;
        
        this->set_elem_solution(dsol);
        this->elem_calculations(false, dres, dummy);
        jac.col(i) = (dres-res0)/delta;
    }
    
    // write the numerical and analytical jacobians
    libMesh::out
    << "Analytical Jacobian: " << std::endl
    << jac0
    << std::endl << std::endl
    << "Numerical Jacobian: " << std::endl
    << jac
    << std::endl << std::endl;
    
    MAST::compare_matrix(jac, jac0, 1.0e-5);
    // set the original solution vector for the element
    this->set_elem_solution(sol);
}


