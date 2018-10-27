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

#ifndef __mast_utility__
#define __mast_utility__


// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/parallel.h"


namespace MAST {
    
    template <typename ValType>
    void transform_to_elem_vector(libMesh::DenseVector<ValType>& v,
                                  const DenseRealVector& v_real);
    
    
    template <>
    inline void transform_to_elem_vector(DenseRealVector& v,
                                         const DenseRealVector& v_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int n = v.size();
        libmesh_assert_equal_to(v_real.size(), n);
        v = v_real;
    }
    
    
    template <>
    inline void transform_to_elem_vector(DenseComplexVector& v,
                                         const DenseRealVector& v_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int n = v.size();
        libmesh_assert_equal_to(v_real.size(), 2*n);
        
        for (unsigned int i=0; i<n; i++)
            v(i) = Complex(v_real(i), v_real(i+n));
    }

    
    
    template <typename ValType>
    void transform_to_elem_matrix(libMesh::DenseMatrix<ValType>& m,
                                  const DenseRealMatrix& m_real);
    
    
    template <>
    inline void transform_to_elem_matrix(DenseRealMatrix& m,
                                         const DenseRealMatrix& m_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int mm = m.m(), nn = m.n();
        libmesh_assert_equal_to(m_real.m(), mm);
        libmesh_assert_equal_to(m_real.n(), nn);
        m = m_real;
    }
    
    
    template <>
    inline void transform_to_elem_matrix(DenseComplexMatrix& m,
                                         const DenseRealMatrix& m_real) {
        // make sure that the real vector is twice the size of the dense vector
        const unsigned int mm = m.m(), nn = m.n();
        libmesh_assert_equal_to(m_real.m(), 2*mm);
        libmesh_assert_equal_to(m_real.n(), 2*nn);

        for (unsigned int i=0; i<mm; i++)
            for (unsigned int j=0; j<nn; j++)
                m(i,j) = Complex(m_real(i,j), -m_real(i,j+nn));
    }

    
    /*!
     *    All calculations in MAST are done using Real numbers. The complex
     *    variables are divided into two unknowns, one each for the real and
     *    imaginary variables. This provides a template method to add a real
     *    or complex vector to the assembled vector.
     */
    template <typename ValType>
    void add_to_assembled_vector(RealVectorX& assembled_vec,
                                 const ValType& elem_vec);
    
    
    /*!
     *    All calculations in MAST are done using Real numbers. The complex
     *    variables are divided into two unknowns, one each for the real and
     *    imaginary variables. This provides a template method to add a real
     *    or complex matrix to the assembled matrix.
     */
    template <typename ValType>
    inline void add_to_assembled_matrix(RealMatrixX& assembled_mat,
                                        const ValType& elem_mat);
    
    
    template <>
    inline void
    add_to_assembled_matrix(RealMatrixX& assembled_mat,
                            const RealMatrixX& elem_mat) {
        assembled_mat += elem_mat;
    }
    
    
    template <>
    inline void
    add_to_assembled_vector(RealVectorX& assembled_vec,
                            const RealVectorX& elem_vec) {
        assembled_vec += elem_vec;
    }
    
    
    template <>
    inline void
    add_to_assembled_matrix(RealMatrixX& assembled_mat,
                            const ComplexMatrixX& elem_mat) {
        
        // make sure the the assembled mat is twice the size of the elem mat
        const unsigned int
        m = (unsigned int)elem_mat.rows(),
        n = (unsigned int)elem_mat.cols();
        libmesh_assert_equal_to(assembled_mat.rows(), m*2);
        libmesh_assert_equal_to(assembled_mat.cols(), n*2);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++) {
                assembled_mat(i,j)     +=  std::real(elem_mat(i,j));
                assembled_mat(i+m,j+n) +=  std::real(elem_mat(i,j));
                assembled_mat(i,j+n)   += -std::imag(elem_mat(i,j));
                assembled_mat(i+m,j)   +=  std::imag(elem_mat(i,j));
            }
    }
    
    
    template <>
    inline void
    add_to_assembled_vector(RealVectorX& assembled_vec,
                            const ComplexVectorX& elem_vec) {
        // make sure the the assembled mat is twice the size of the elem mat
        const unsigned int n = (unsigned int)elem_vec.size();
        libmesh_assert_equal_to(assembled_vec.size(), n*2);
        
        for (unsigned int i=0; i<n; i++) {
            assembled_vec(i)     +=   std::real(elem_vec(i));
            assembled_vec(i+n)   +=  std::imag(elem_vec(i));
        }
    }
    
    
    
    inline void
    copy (DenseRealMatrix& m1, const RealMatrixX& m2) {
        
        const unsigned int m=(unsigned int)m2.rows(), n=(unsigned int)m2.cols();
        m1.resize(m, n);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                m1(i,j) = m2(i,j);
    }

    

    inline void
    copy (RealMatrixX& m2, const DenseRealMatrix& m1) {
        
        const unsigned int m=(unsigned int)m1.m(), n=(unsigned int)m1.n();
        m2.setZero(m, n);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                m2(i,j) = m1(i,j);
    }

    
    
    inline void
    copy (DenseRealVector& v1, const RealVectorX& v2) {
        
        const unsigned int m=(unsigned int)v2.rows();
        v1.resize(m);
        for (unsigned int i=0; i<m; i++)
            v1(i) = v2(i);
    }

    
    inline void
    copy (RealVectorX& v1, const DenseRealVector& v2) {
        
        const unsigned int m=(unsigned int)v2.size();
        
        v1 = RealVectorX::Zero(m);
        
        for (unsigned int i=0; i<m; i++)
            v1(i) = v2(i);
    }


    
    inline void
    parallel_sum (const libMesh::Parallel::Communicator& c,
                  RealMatrixX& mat) {
        
        const unsigned int
        m =  (unsigned int)mat.rows(),
        n =  (unsigned int)mat.cols();
    
        std::vector<Real> vals(m*n);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                vals[m*j+i] = mat(i,j);
        
        c.sum(vals);
        
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                mat(i,j) = vals[m*j+i];
    }

    
    
    inline void
    parallel_sum (const libMesh::Parallel::Communicator& c,
                  ComplexMatrixX& mat) {
        
        const unsigned int
        m =  (unsigned int)mat.rows(),
        n =  (unsigned int)mat.cols();
        
        std::vector<Complex> vals(m*n);
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                vals[m*j+i] = mat(i,j);
        
        c.sum(vals);
        
        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                mat(i,j) = vals[m*j+i];
        
    }
}




#endif // __mast_utility__

