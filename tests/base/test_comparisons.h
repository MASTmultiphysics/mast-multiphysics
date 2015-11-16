/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast_test_comparisons_h__
#define __mast_test_comparisons_h__

// MAST includes
#include "base/mast_data_types.h"


namespace MAST {
    
    inline bool
    compare_value(const Real v0, const Real v, const Real tol) {
        
        bool pass = true;
        
        if (!boost::test_tools::check_is_close(v0,
                                               v,
                                               boost::test_tools::percent_tolerance<Real>(tol))) {
            BOOST_TEST_MESSAGE ("Failed comparison: "
                                << "expected: " << v0<< "  , "
                                << "computed: " << v << " : "
                                << "diff: " << v0 - v << " , "
                                << "tol: " << tol);
            pass = false;
        }
        
        return pass;
    }
    
    
    
    inline bool
    compare_vector(const RealVectorX& v0, const RealVectorX& v, const Real tol) {
        
        unsigned int
        v0_size = (unsigned int) v0.rows();
        libmesh_assert_equal_to(v0_size,  v.size());
        
        
        bool pass = true;
        for (unsigned int i=0; i<v0_size; i++) {
            if (!boost::test_tools::check_is_close(v0(i),
                                                   v(i),
                                                   boost::test_tools::percent_tolerance(tol))) {
                BOOST_TEST_MESSAGE("Failed comparison at i = ("
                                   << i << ") : "
                                   << "expected: " << v0(i) << "  , "
                                   << "computed: " << v(i) << " : "
                                   << "diff: " << v0(i) - v(i) << " , "
                                   << "tol: " << tol);
                pass = false;
            }
        }
        
        return pass;
    }
    
    
    
    
    inline bool
    compare_matrix(const RealMatrixX& m0, const RealMatrixX& m, const Real tol) {
        
        unsigned int
        m0_rows = (unsigned int) m0.rows(),
        m0_cols = (unsigned int) m0.cols();
        libmesh_assert_equal_to(m0_rows,  m.rows());
        libmesh_assert_equal_to(m0_cols,  m.cols());
        
        
        bool pass = true;
        for (unsigned int i=0; i<m0_rows; i++) {
            for (unsigned int j=0; j<m0_cols; j++)
                if (!boost::test_tools::check_is_close(m0(i,j),
                                                       m(i,j),
                                                       boost::test_tools::percent_tolerance(tol))) {
                    BOOST_TEST_MESSAGE("Failed comparison at (i,j) = ("
                                       << i << ", " << j << ") : "
                                       << "expected: " << m0(i,j) << "  , "
                                       << "computed: " << m(i,j) << " : "
                                       << "diff: " << m0(i,j) - m(i,j) << " , "
                                       << "tol: " << tol);
                    pass = false;
                }
        }
        
        return pass;
    }
    
}

#endif //__mast_test_compasisons_h__
