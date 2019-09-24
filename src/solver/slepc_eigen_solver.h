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

#ifndef __mast__slepc_eigen_solver__
#define __mast__slepc_eigen_solver__

// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/slepc_eigen_solver.h"


namespace MAST {
    
    /*!
     *  This class inherits from libMesh::SlepcEigenSolver<Real> and implements a
     *  method for retriving the real and imaginary components of the eigenvector, 
     *  which the libMesh interface does not provide.
     */
    
    class  SlepcEigenSolver:
    public libMesh::SlepcEigenSolver<Real> {
      
    public:
        SlepcEigenSolver(const libMesh::Parallel::Communicator & comm_in);
        
        /**
         * This function returns the real and imaginary part of the
         * ith eigenvalues.
         */
        virtual std::pair<Real, Real>
        get_eigenvalue (unsigned int i);

        
        /**
         * This function returns the real and imaginary part of the
         * ith eigenvalue and copies the respective eigenvector to the
         * solution vector. Note that also in case of purely real matrix
         * entries the eigenpair may be complex values.
         */
        virtual std::pair<Real, Real>
        get_eigenpair (unsigned int i,
                       libMesh::NumericVector<Real> &eig_vec,
                       libMesh::NumericVector<Real> *eig_vec_im = libmesh_nullptr);


    };
}


#endif // __mast__slepc_eigen_solver__

