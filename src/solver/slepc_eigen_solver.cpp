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

// MAST includes
#include "solver/slepc_eigen_solver.h"

// libMesh includes
#include "libmesh/petsc_vector.h"
#include "libmesh/enum_eigen_solver_type.h"


MAST::SlepcEigenSolver::SlepcEigenSolver(const libMesh::Parallel::Communicator & comm_in):
libMesh::SlepcEigenSolver<Real>(comm_in) {
    
}




std::pair<Real, Real>
MAST::SlepcEigenSolver::get_eigenvalue(unsigned int i) {
    
    PetscErrorCode ierr=0;
    
    PetscReal re, im;
    
    ierr = EPSGetEigenvalue(eps(), i, &re, &im);
    
    CHKERRABORT(this->comm().get(), ierr);
    
    return std::make_pair(re, im);
}




std::pair<Real, Real>
MAST::SlepcEigenSolver::get_eigenpair(unsigned int i,
                                      libMesh::NumericVector<Real> &eig_vec,
                                      libMesh::NumericVector<Real> *eig_vec_im) {
    
    // make sure that for non-Hermitian problems with real matrices
    // a vector is provided for imaginary part.
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(eig_vec_im == NULL);
#else
    if (this->eigen_problem_type() == libMesh::HEP ||
        this->eigen_problem_type() == libMesh::GHEP)
        libmesh_assert(eig_vec_im == NULL);
#endif
    
    PetscErrorCode ierr=0;
    
    PetscReal re, im;
    
    // Make sure the NumericVector passed in is really a PetscVector
    libMesh::PetscVector<Real>
    *v_re = libMesh::cast_ptr<libMesh::PetscVector<Real>*>(&eig_vec),
    *v_im = NULL;
    
    // real and imaginary part of the ith eigenvalue.
    PetscScalar kr, ki;
    
    eig_vec.close();
    if (eig_vec_im) {
        
        eig_vec_im = libMesh::libmesh_cast_ptr<libMesh::PetscVector<Real>*>(eig_vec_im);
        eig_vec_im->close();
        
        // now get the eigenvector
        ierr = EPSGetEigenpair(eps(), i, &kr, &ki, v_re->vec(), v_im->vec());
    }
    else
        ierr = EPSGetEigenpair(eps(), i, &kr, &ki, v_re->vec(), PETSC_NULL);
    
    CHKERRABORT(this->comm().get(), ierr);
    
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif
    
    return std::make_pair(re, im);
}

