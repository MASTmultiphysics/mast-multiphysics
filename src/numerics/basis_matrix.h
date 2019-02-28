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


#ifndef __mast__structural_basis_matrix_h__
#define __mast__structural_basis_matrix_h__

// C++ includes
#include <vector>

// libmesh includes
#include "libmesh/shell_matrix.h"
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    template <typename T>
    class BasisMatrix:
    public libMesh::ShellMatrix<T>
    {
    public:
        BasisMatrix(const libMesh::Parallel::Communicator &comm_in);
        
        virtual ~BasisMatrix();
        
        
        /**
         * @returns \p m, the row-dimension of the matrix where the marix is
         * \f$ M \times N \f$.
         */
        virtual libMesh::numeric_index_type m () const {
            
            libmesh_assert(modes.size() > 0);
            return modes[0]->size();
        }
        
        /**
         * @returns \p n, the column-dimension of the matrix where the marix
         * is \f$ M \times N \f$.
         */
        virtual libMesh::numeric_index_type n () const {
            
            libmesh_assert(modes.size() > 0);
            return (libMesh::numeric_index_type) modes.size();
        }
        
        
        /**
         * Multiplies the matrix with \p arg and stores the result in \p
         * dest.
         */
        template <typename VecType>
        void vector_mult (libMesh::NumericVector<T>& dest,
                          const VecType& arg) const {
            
            libmesh_assert(modes.size() > 0);
            libmesh_assert_equal_to(m(), dest.size());
            libmesh_assert_equal_to(n(), arg.size());
            
            dest.zero();
            for (unsigned int i=0; i<n(); i++)
                dest.add(arg(i), *(modes[i]));
        }
        
        /**
         * Multiplies the matrix with \p arg and stores the result in \p
         * dest.
         */
        virtual void
        vector_mult (libMesh::NumericVector<T>& dest,
                     const libMesh::NumericVector<T>& arg) const {
            
            // not defined for multiplcation with libMesh::NumericVector
            libmesh_assert(false);
        }
        
        
        /**
         * Multiplies the matrix with \p arg and stores the result in \p
         * dest.
         */
        template <typename VecType>
        void
        vector_mult_transpose (VecType& dest,
                               const libMesh::NumericVector<T>& arg) const {
            
            libmesh_assert(modes.size() > 0);
            libmesh_assert_equal_to(m(), arg.size());
            libmesh_assert_equal_to(n(), dest.size());
            
            dest.setZero(n());
            for (unsigned int i=0; i<n(); i++)
                dest(i) = arg.dot(*(modes[i]));
        }
        
        
        /**
         * Multiplies the transpose of matrix with \p arg and stores the
         * result in \p dest.
         */
        virtual void
        vector_mult_transpose (libMesh::NumericVector<T>& dest,
                               const libMesh::NumericVector<T>& arg) const {
            
            // not defined for multiplcation with libMesh::NumericVector
            libmesh_assert(false);
        }
        
        /**
         * Multiplies the matrix with \p arg and adds the result to \p dest.
         */
        virtual void
        vector_mult_add (libMesh::NumericVector<T>& dest,
                         const libMesh::NumericVector<T>& arg) const {
            
            // not defined for multiplcation with libMesh::NumericVector
            libmesh_assert(false);
        }
        
        
        /**
         * Copies the diagonal part of the matrix into \p dest.
         */
        virtual void
        get_diagonal (libMesh::NumericVector<T>& dest) const {
            
            // not defined for multiplcation with libMesh::NumericVector
            libmesh_assert(false);
        }
        
        
        /*!
         *  Returns the vector that defines the \p i^th basis vector
         */
        
        virtual libMesh::NumericVector<T>&
        basis(unsigned int i) {
            
            // not defined for multiplcation with libMesh::NumericVector
            libmesh_assert(modes.size() > 0);
            libmesh_assert_less(i, modes.size());
            
            return *(modes[i]);
        }
        
        /*!
         *   vector of modes
         */
        std::vector<libMesh::NumericVector<T>*> modes;
    };
    
}

#endif // __mast__structural_basis_matrix_h__
