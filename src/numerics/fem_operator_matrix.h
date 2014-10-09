/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
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

#ifndef __mast__fem_operator_matrix__
#define __mast__fem_operator_matrix__

// C++ includes
#include <vector>


// MAST includes
#include "base/mast_data_types.h"




namespace MAST {
    
    class FEMOperatorMatrix
    {
    public:
        FEMOperatorMatrix();
        
        
        virtual ~FEMOperatorMatrix();
        
        
        /*!
         *   clears the data structures
         */
        void clear();
        
        
        unsigned int m() const {return _n_interpolated_vars;}
        
        unsigned int n() const {return _n_discrete_vars*_n_dofs_per_var;}
        
        /*!
         *   this initializes the operator for number of rows and variables, assuming
         *   that all variables has the same number of dofs. This is typically the case
         *   for structural strain operator matrices. Note that when this method is used
         *   the user must set the matrix entries by calling set_shape_functions
         */
        void reinit(unsigned int n_interpolated_vars,
                    unsigned int n_discrete_vars,
                    unsigned int n_discrete_dofs_per_var);
        
        /*!
         *   sets the shape function values for the block corresponding to
         *   \par interpolated_var and \par discrete_var. This means that the row
         *   \par interpolated_var, the value in columns
         *   \par discrete_vars*n_discrete_dofs_per_var - (discrete_vars+1)*n_discrete_dofs_per_var-1)
         *    will be set equal to \par shape_func .
         */
        void set_shape_function(unsigned int interpolated_var,
                                unsigned int discrete_var,
                                const RealVectorX& shape_func);
        
        /*!
         *   this initializes all variables to use the same interpolation function.
         *   It is assumed that the number of discrete vars is same as the number of
         *   interpolated vars. This is typically the case for fluid elements and
         *   for structural element inertial matrix calculations
         */
        void reinit(unsigned int n_interpolated_vars,
                    const RealVectorX& shape_func);
        
        /*!
         *   res = [this] * v
         */
        template <typename T>
        void vector_mult(T& res, const T& v) const;
        
        
        /*!
         *   res = v^T * [this]
         */
        template <typename T>
        void vector_mult_transpose(T& res, const T& v) const;
        
        
        /*!
         *   [R] = [this] * [M]
         */
        template <typename T>
        void right_multiply(T& r, const T& m) const;
        
        
        /*!
         *   [R] = [this]^T * [M]
         */
        template <typename T>
        void right_multiply_transpose(T& r, const T& m) const;
        
        
        /*!
         *   [R] = [this]^T * [M]
         */
        template <typename T>
        void right_multiply_transpose(T& r, const MAST::FEMOperatorMatrix& m) const;
        
        
        /*!
         *   [R] = [M] * [this]
         */
        template <typename T>
        void left_multiply(T& r, const T& m) const;
        
        
        /*!
         *   [R] = [M] * [this]^T
         */
        template <typename T>
        void left_multiply_transpose(T& r, const T& m) const;
        
        
    protected:
        
        /*!
         *    number of rows of the operator
         */
        unsigned int _n_interpolated_vars;
        
        /*!
         *    number of discrete variables in the system
         */
        unsigned int _n_discrete_vars;
        
        /*!
         *    number of dofs for each variable
         */
        unsigned int _n_dofs_per_var;
        
        /*!
         *    stores the shape function values that defines the coupling
         *    of i_th interpolated var and j_th discrete var. Stored in
         *    column major format. NULL, if values are zero, otherwise the
         *    value is set in the vector.
         */
        std::vector<RealVectorX*>  _var_shape_functions;
    };
    
}



inline
void
MAST::FEMOperatorMatrix::clear() {
    
    _n_interpolated_vars = 0;
    _n_discrete_vars     = 0;
    _n_dofs_per_var      = 0;
    
    // iterate over the shape function entries and delete the non-NULL values
    std::vector<RealVectorX*>::iterator it = _var_shape_functions.begin(),
    end = _var_shape_functions.end();
    
    for ( ; it!=end; it++)
        if ( *it != NULL)
            delete *it;
    
    _var_shape_functions.clear();
}




inline
void
MAST::FEMOperatorMatrix::
reinit(unsigned int n_interpolated_vars,
       unsigned int n_discrete_vars,
       unsigned int n_discrete_dofs_per_var) {
    
    this->clear();
    _n_interpolated_vars = n_interpolated_vars;
    _n_discrete_vars = n_discrete_vars;
    _n_dofs_per_var = n_discrete_dofs_per_var;
    _var_shape_functions.resize(_n_interpolated_vars*_n_discrete_vars);
    for (unsigned int i=0; i<_var_shape_functions.size(); i++)
        _var_shape_functions[i] = NULL;
}



inline
void
MAST::FEMOperatorMatrix::
set_shape_function(unsigned int interpolated_var,
                   unsigned int discrete_var,
                   const RealVectorX& shape_func) {
    
    // make sure that reinit has been called.
    libmesh_assert(_var_shape_functions.size());
    
    // also make sure that the specified indices are within bounds
    libmesh_assert(interpolated_var < _n_interpolated_vars);
    libmesh_assert(discrete_var < _n_discrete_vars);
    
    RealVectorX* vec =
    _var_shape_functions[discrete_var*_n_interpolated_vars+interpolated_var];
    
    if (!vec)
    {
        vec = new RealVectorX;
        _var_shape_functions[discrete_var*_n_interpolated_vars+interpolated_var] = vec;
    }
    
    *vec = shape_func;
}



inline
void
MAST::FEMOperatorMatrix::
reinit(unsigned int n_vars,
       const RealVectorX& shape_func) {
    
    this->clear();
    
    _n_interpolated_vars = n_vars;
    _n_discrete_vars = n_vars;
    _n_dofs_per_var = (unsigned int)shape_func.size();
    _var_shape_functions.resize(n_vars*n_vars);
    
    for (unsigned int i=0; i<n_vars; i++)
    {
        RealVectorX*  vec = new RealVectorX;
        *vec = shape_func;
        _var_shape_functions[i*n_vars+i] = vec;
    }
}



template <typename T>
inline
void
MAST::FEMOperatorMatrix::
vector_mult(T& res, const T& v) const {
    
    libmesh_assert_equal_to(res.size(), _n_interpolated_vars);
    libmesh_assert_equal_to(v.size(), n());
    
    res.zero();
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) // check if this is non-NULL
                for (unsigned int k=0; k<_n_dofs_per_var; k++)
                    res(i) +=
                    (*_var_shape_functions[index])(k) * v(j*_n_dofs_per_var+k);
        }
}


template <typename T>
inline
void
MAST::FEMOperatorMatrix::
vector_mult_transpose(T& res, const T& v) const {
    
    libmesh_assert_equal_to(res.size(), n());
    libmesh_assert_equal_to(v.size(), _n_interpolated_vars);
    
    res.setZero(res.size());
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) // check if this is non-NULL
                for (unsigned int k=0; k<_n_dofs_per_var; k++)
                    res(j*_n_dofs_per_var+k) +=
                    (*_var_shape_functions[index])(k) * v(i);
        }
}



template <typename T>
inline
void
MAST::FEMOperatorMatrix::
right_multiply(T& r, const T& m) const {
    
    libmesh_assert_equal_to(r.rows(), _n_interpolated_vars);
    libmesh_assert_equal_to(r.cols(), m.cols());
    libmesh_assert_equal_to(m.rows(), n());
    
    r.setZero();
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.cols(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(i,l) +=
                        (*_var_shape_functions[index])(k) * m(j*_n_dofs_per_var+k,l);
            }
        }
}




template <typename T>
inline
void
MAST::FEMOperatorMatrix::
right_multiply_transpose(T& r, const T& m) const {
    
    libmesh_assert_equal_to(r.rows(), n());
    libmesh_assert_equal_to(r.cols(), m.cols());
    libmesh_assert_equal_to(m.rows(), _n_interpolated_vars);
    
    r.setZero(r.rows(), r.cols());
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.cols(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(j*_n_dofs_per_var+k,l) +=
                        (*_var_shape_functions[index])(k) * m(i,l);
            }
        }
}



template <typename T>
inline
void
MAST::FEMOperatorMatrix::
right_multiply_transpose(T& r, const MAST::FEMOperatorMatrix& m) const {
    
    libmesh_assert_equal_to(r.rows(), n());
    libmesh_assert_equal_to(r.cols(), m.n());
    libmesh_assert_equal_to(_n_interpolated_vars, m._n_interpolated_vars);
    
    r.setZero();
    unsigned int index_i, index_j = 0;
    
    for (unsigned int i=0; i<_n_discrete_vars; i++) // row of result
        for (unsigned int j=0; j<m._n_discrete_vars; j++) // column of result
            for (unsigned int k=0; k<_n_interpolated_vars; k++) {
                index_i = i*_n_interpolated_vars+k;
                index_j = j*m._n_interpolated_vars+k;
                if (_var_shape_functions[index_i] &&
                    m._var_shape_functions[index_j]) { // if shape function exists for both
                    const RealVectorX &n1 = *_var_shape_functions[index_i],
                    &n2 = *m._var_shape_functions[index_j];
                    for (unsigned int i_n1=0; i_n1<n1.size(); i_n1++)
                        for (unsigned int i_n2=0; i_n2<n2.size(); i_n2++)
                            r (i*_n_dofs_per_var+i_n1,
                               j*m._n_dofs_per_var+i_n2) += n1(i_n1) * n2(i_n2);
                }
            }
}




template <typename T>
inline
void
MAST::FEMOperatorMatrix::
left_multiply(T& r, const T& m) const {
    
    libmesh_assert_equal_to(r.rows(), m.rows());
    libmesh_assert_equal_to(r.cols(), n());
    libmesh_assert_equal_to(m.cols(), _n_interpolated_vars);
    
    r.setZero(r.rows(), r.cols());
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.rows(); l++) // rows of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(l,j*_n_dofs_per_var+k) +=
                        (*_var_shape_functions[index])(k) * m(l,i);
            }
        }
}



template <typename T>
inline
void
MAST::FEMOperatorMatrix::
left_multiply_transpose(T& r, const T& m) const {
    
    libmesh_assert_equal_to(r.rows(), m.rows());
    libmesh_assert_equal_to(r.cols(), _n_interpolated_vars);
    libmesh_assert_equal_to(m.cols(), n());
    
    r.zero();
    unsigned int index = 0;
    
    for (unsigned int i=0; i<_n_interpolated_vars; i++) // row
        for (unsigned int j=0; j<_n_discrete_vars; j++) { // column of operator
            index = j*_n_interpolated_vars+i;
            if (_var_shape_functions[index]) { // check if this is non-NULL
                for (unsigned int l=0; l<m.rows(); l++) // column of matrix
                    for (unsigned int k=0; k<_n_dofs_per_var; k++)
                        r(l,i) +=
                        (*_var_shape_functions[index])(k) * m(l,j*_n_dofs_per_var+k);
            }
        }
}




#endif // __mast__fem_operator_matrix__

