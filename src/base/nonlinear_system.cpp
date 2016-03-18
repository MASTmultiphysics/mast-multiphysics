/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

// C++ includes
#include <vector>

// MAST includes
#include "base/nonlinear_system.h"
#include "base/physics_discipline_base.h"
#include "base/eigensystem_assembly.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/slepc_eigen_solver.h"


MAST::NonlinearSystem::NonlinearSystem(libMesh::EquationSystems& es,
                                       const std::string& name,
                                       const unsigned int number):
libMesh::NonlinearImplicitSystem(es, name, number),
_initialize_B_matrix                  (false),
matrix_A                              (NULL),
matrix_B                              (NULL),
condensed_matrix_A                    (NULL),
condensed_matrix_B                    (NULL),
eigen_solver                          (NULL),
_condensed_dofs_initialized           (false),
_exchange_A_and_B                     (false),
_n_requested_eigenpairs               (0),
_n_converged_eigenpairs               (0),
_n_iterations                         (0),
_is_generalized_eigenproblem          (false),
_eigen_problem_type                   (libMesh::NHEP),
_eigenproblem_assemble_system_object  (NULL) {
    
}



MAST::NonlinearSystem::~NonlinearSystem() {
    
    this->clear();
}





void
MAST::NonlinearSystem::clear() {
    
    
    // delete the matricies
    delete matrix_A;
    if (matrix_B) delete matrix_B;
    
    // NULL-out the matricies.
    matrix_A = NULL;
    matrix_B = NULL;
    
    // clear the solver
    eigen_solver->clear();
    
    libMesh::NonlinearImplicitSystem::clear();
}


void
MAST::NonlinearSystem::set_eigenproblem_type (libMesh::EigenProblemType ept) {
    
    _eigen_problem_type = ept;
    
}




void
MAST::NonlinearSystem::init_data () {
    
    // initialize parent data
    libMesh::NonlinearImplicitSystem::init_data();
    
    // define the type of eigenproblem
    if (_eigen_problem_type == libMesh::GNHEP ||
        _eigen_problem_type == libMesh::GHEP  ||
        _eigen_problem_type == libMesh::GHIEP)
        _is_generalized_eigenproblem = true;
    
    
    libMesh::DofMap& dof_map = this->get_dof_map();
    dof_map.compute_sparsity(this->get_mesh());
    
    // build the system matrix
    matrix_A = libMesh::SparseMatrix<Real>::build(this->comm()).release();
    dof_map.attach_matrix(*matrix_A);
    matrix_A->init();
    matrix_A->zero();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix) {
        
        matrix_B = libMesh::SparseMatrix<Real>::build(this->comm()).release();
        dof_map.attach_matrix(*matrix_B);
        matrix_B->init();
        matrix_B->zero();
    }
    
    this->condensed_matrix_A.reset(libMesh::SparseMatrix<Real>::build(this->comm()).release());
    this->condensed_matrix_B.reset(libMesh::SparseMatrix<Real>::build(this->comm()).release());
    
    eigen_solver.reset(libMesh::EigenSolver<Real>::build(this->comm()).release());
    if (libMesh::on_command_line("--solver_system_names")) {
        
        EPS eps =  dynamic_cast<libMesh::SlepcEigenSolver<Real>*>(eigen_solver.get())->eps();
        std::string nm = this->name() + "_";
        EPSSetOptionsPrefix(eps, nm.c_str());
    }
    eigen_solver->set_eigenproblem_type(_eigen_problem_type);
    
}




void MAST::NonlinearSystem::reinit () {
    // initialize parent data
    libMesh::NonlinearImplicitSystem::reinit();
    
    // Clear the matrices
    matrix_A->clear();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix)
        matrix_B->clear();
    
    libMesh::DofMap& dof_map = this->get_dof_map();
    
    // Clear the sparsity pattern
    dof_map.clear_sparsity();
    
    // Compute the sparsity pattern for the current
    // mesh and DOF distribution.  This also updates
    // both matrices, \p DofMap now knows them
    dof_map.compute_sparsity(this->get_mesh());
    
    matrix_A->init();
    matrix_A->zero();
    
    if (_is_generalized_eigenproblem || _initialize_B_matrix)
    {
        matrix_B->init();
        matrix_B->zero();
    }
}




void
MAST::NonlinearSystem::eigenproblem_solve() {
    
    
    START_LOG("eigensolve()", "NonlinearSystem");
    
    // A reference to the EquationSystems
    libMesh::EquationSystems& es = this->get_equation_systems();
    
    // check that necessary parameters have been set
    libmesh_assert(_n_requested_eigenpairs);
    
    es.parameters.set<unsigned int>("eigenpairs")    =   _n_requested_eigenpairs;
    es.parameters.set<unsigned int>("basis vectors") = 5*_n_requested_eigenpairs;
    
    // Get the tolerance for the solver and the maximum
    // number of iterations. Here, we simply adopt the linear solver
    // specific parameters.
    const Real
    tol    = es.parameters.get<Real>("linear solver tolerance");
    
    const unsigned int
    maxits = es.parameters.get<unsigned int>("linear solver maximum iterations"),
    nev    = es.parameters.get<unsigned int>("eigenpairs"),
    ncv    = es.parameters.get<unsigned int>("basis vectors");
    
    std::pair<unsigned int, unsigned int> solve_data;
    
    libMesh::SparseMatrix<Real>
    *eig_A  = NULL,
    *eig_B  = NULL;
    
    // assemble the matrices
    this->assemble_eigensystem();
    
    // If we haven't initialized any condensed dofs,
    // just use the default eigen_system
    if (!_condensed_dofs_initialized) {
        
        // call the solver depending on the type of eigenproblem
        if (generalized())  {
            
            //in case of a generalized eigenproblem
            
            // exchange the matrices if requested by the user
            if (!_exchange_A_and_B) {
                eig_A  =  matrix_A;
                eig_B  =  matrix_B;
            }
            else {
                eig_B  =  matrix_A;
                eig_A  =  matrix_B;
            }
            
            solve_data = eigen_solver->solve_generalized (*eig_A,
                                                          *eig_B,
                                                          nev,
                                                          ncv,
                                                          tol,
                                                          maxits);
        }
        else {
            
            libmesh_assert (!matrix_B);
            
            //in case of a standard eigenproblem
            solve_data = eigen_solver->solve_standard (*matrix_A,
                                                       nev,
                                                       ncv,
                                                       tol,
                                                       maxits);
        }
    }
    else {
        
        // If we reach here, then there should be some non-condensed dofs
        libmesh_assert(!_local_non_condensed_dofs_vector.empty());
        
        // Now condense the matrices
        matrix_A->create_submatrix(*condensed_matrix_A,
                                   _local_non_condensed_dofs_vector,
                                   _local_non_condensed_dofs_vector);
        
        if (generalized()) {
            
            matrix_B->create_submatrix(*condensed_matrix_B,
                                       _local_non_condensed_dofs_vector,
                                       _local_non_condensed_dofs_vector);
        }
        
        // call the solver depending on the type of eigenproblem
        if ( generalized() ) {
            
            //in case of a generalized eigenproblem
            
            // exchange the matrices if requested by the user
            if (!_exchange_A_and_B) {
                eig_A  =  condensed_matrix_A.get();
                eig_B  =  condensed_matrix_B.get();
            }
            else {
                eig_B  =  condensed_matrix_A.get();
                eig_A  =  condensed_matrix_B.get();
            }
            
            solve_data = eigen_solver->solve_generalized(*eig_A,
                                                         *eig_B,
                                                         nev,
                                                         ncv,
                                                         tol,
                                                         maxits);
        }
        else {
            
            libmesh_assert (!matrix_B);
            
            //in case of a standard eigenproblem
            solve_data = eigen_solver->solve_standard (*condensed_matrix_A,
                                                       nev,
                                                       ncv,
                                                       tol,
                                                       maxits);
        }
    }
    
    _n_converged_eigenpairs = solve_data.first;
    _n_iterations           = solve_data.second;
    
    STOP_LOG("eigensolve()", "NonlinearSystem");
    
    return;
}




void
MAST::NonlinearSystem::get_eigenpair(unsigned int i,
                                     Real&  re,
                                     Real&  im,
                                     libMesh::NumericVector<Real>& vec_re,
                                     libMesh::NumericVector<Real>* vec_im) {
    
    START_LOG("get_eigenpair()", "NonlinearSystem");
    std::pair<Real, Real>
    val;
    
    // imaginary component is needed only when the problem is non-Hermitian
    if (_eigen_problem_type == libMesh::HEP ||
        _eigen_problem_type == libMesh::GHEP)
        libmesh_assert (vec_im == NULL);
    
    
    // If we haven't initialized any condensed dofs,
    // just use the default eigen_system
    if (!_condensed_dofs_initialized) {
        
        // call the eigen_solver get_eigenpair method
        val   = this->eigen_solver->get_eigenpair (i, vec_re, vec_im);
        
        if (!_exchange_A_and_B) {
            re   = val.first;
            im   = val.second;
        }
        else {
            Complex complex_val (val.first, val.second);
            complex_val = 1./complex_val;
            re = complex_val.real();
            im = complex_val.imag();
        }
    }
    else {
        
        // If we reach here, then there should be some non-condensed dofs
        libmesh_assert(!_local_non_condensed_dofs_vector.empty());
        
        std::auto_ptr< libMesh::NumericVector<Real> >
        temp_re(libMesh::NumericVector<Real>::build(this->comm()).release()),
        temp_im;
        
        unsigned int
        n_local   = (unsigned int)_local_non_condensed_dofs_vector.size(),
        n         = n_local;
        this->comm().sum(n);
        
        // initialize the vectors
        temp_re->init (n, n_local, false, libMesh::PARALLEL);
        
        // imaginary only if the problem is non-Hermitian
        if (vec_im) {
            
            temp_im.reset(libMesh::NumericVector<Real>::build(this->comm()).release());
            temp_im->init (n, n_local, false, libMesh::PARALLEL);
        }
        

        // call the eigen_solver get_eigenpair method
        val   = this->eigen_solver->get_eigenpair (i, *temp_re, temp_im.get());
        
        if (!_exchange_A_and_B) {
            re   = val.first;
            im   = val.second;
        }
        else {
            Complex complex_val (val.first, val.second);
            complex_val = 1./complex_val;
            re = complex_val.real();
            im = complex_val.imag();
        }

        
        // Now map temp to solution. Loop over local entries of local_non_condensed_dofs_vector
        // the real part
        vec_re.zero();
        for (unsigned int j=0; j<_local_non_condensed_dofs_vector.size(); j++) {
            
            unsigned int index = _local_non_condensed_dofs_vector[j];
            vec_re.set(index,(*temp_re)(temp_re->first_local_index()+j));
        }
        vec_re.close();
        
        // now the imaginary part if it was provided
        if (vec_im) {
            
            vec_im->zero();
            
            for (unsigned int j=0; j<_local_non_condensed_dofs_vector.size(); j++) {
                
                unsigned int index = _local_non_condensed_dofs_vector[j];
                vec_im->set(index,(*temp_im)(temp_im->first_local_index()+j));
            }
            
            vec_im->close();
        }
    }
    
    // scale the eigenvector so that it has a unit inner product
    // with the B matrix.
    switch (_eigen_problem_type) {
        case libMesh::HEP: {
            Real
            v = vec_re.dot(vec_re);
            
            // make sure that v is a nonzero value
            libmesh_assert_greater(v, 0.);
            vec_re.scale(1./std::sqrt(v));
        }
            break;
            
        case libMesh::GHEP: {
            
            std::auto_ptr<libMesh::NumericVector<Real> >
            tmp(vec_re.zero_clone().release());
            
            // inner product with respect to B matrix
            matrix_B->vector_mult(*tmp, vec_re);

            Real
            v = tmp->dot(vec_re);
            
            // make sure that v is a nonzero value
            libmesh_assert_greater(v, 0.);
            vec_re.scale(1./std::sqrt(v));
        }
            break;
            
            
        case libMesh::NHEP:   // to be implemented
        case libMesh::GNHEP:  // to be implemented
        case libMesh::GHIEP:  // to be implemented
        default:
            libmesh_error(); // should not get here
    }

    this->update();
    
    STOP_LOG("get_eigenpair()", "NonlinearSystem");
}





void
MAST::NonlinearSystem::
eigenproblem_sensitivity_solve (const libMesh::ParameterVector& parameters,
                                std::vector<Real>& sens) {
    
    // make sure that eigensolution is already available
    libmesh_assert(_n_converged_eigenpairs);
    
    // the sensitivity is calculated based on the inner product of the left and
    // right eigen vectors.
    //
    //        y^T [A] x - lambda y^T [B] x = 0
    //    where y and x are the left and right eigenvectors, respectively.
    //    Therefore,
    //        d lambda/dp = (y^T (d[A]/dp - lambda d[B]/dp) x) / (y^T [B] x)
    //
    //    the denominator remain constant for all sensitivity calculations.
    //
    std::vector<Real> denom(_n_converged_eigenpairs, 0.);
    sens.resize(_n_converged_eigenpairs*parameters.size(), 0.);
    
    std::vector<libMesh::NumericVector<Real>*>
    x_right (_n_converged_eigenpairs),
    x_left  (_n_converged_eigenpairs);
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    tmp     (this->solution->zero_clone().release());

    std::vector<Real>
    eig (_n_converged_eigenpairs);
    
    Real
    re  = 0.,
    im  = 0.;
    
    
    for (unsigned int i=0; i<_n_converged_eigenpairs; i++) {

        x_right[i] = (this->solution->zero_clone().release());
        
        
        switch (_eigen_problem_type) {
                
            case libMesh::HEP: {
                // right and left eigenvectors are same
                // imaginary part of eigenvector for real matrices is zero
                this->get_eigenpair(i, re, im, *x_right[i], NULL);
                denom[i] = x_right[i]->dot(*x_right[i]);               // x^H x
                eig[i]   = re;
            }
                break;
                
            case libMesh::GHEP: {
                // imaginary part of eigenvector for real matrices is zero
                this->get_eigenpair(i, re, im, *x_right[i], NULL);
                matrix_B->vector_mult(*tmp, *x_right[i]);
                denom[i] = x_right[i]->dot(*tmp);                  // x^H B x
                eig[i]   = re;
            }
                break;
                
            default:
                // to be implemented for the non-Hermitian problems
                libmesh_error();
                break;
        }
    }
    
    unsigned int
    num = 0;
    
    for (unsigned int p=0; p<parameters.size(); p++) {
        
        // calculate sensitivity of matrix quantities
        this->solution->zero();
        this->assemble_eigensystem_sensitivity(parameters, p);
        
        // now calculate sensitivity of each eigenvalue for the parameter
        for (unsigned int i=0; i<_n_converged_eigenpairs; i++) {
            
            num = p*_n_converged_eigenpairs+i;

            switch (_eigen_problem_type) {
                    
                case libMesh::HEP: {
                    
                    matrix_A->vector_mult(*tmp, *x_right[i]);
                    sens[num] = x_right[i]->dot(*tmp);                  // x^H A' x
                    sens[num]-= eig[i] * x_right[i]->dot(*x_right[i]);  // - lambda x^H x
                    sens[num] /= denom[i];                              // x^H x
                }
                    break;
                    
                case libMesh::GHEP: {
                    
                    matrix_A->vector_mult(*tmp, *x_right[i]);
                    sens[num] = x_right[i]->dot(*tmp);              // x^H A' x
                    matrix_B->vector_mult(*tmp, *x_right[i]);
                    sens[num]-= eig[i] * x_right[i]->dot(*tmp);     // - lambda x^H B' x
                    sens[num] /= denom[i];                          // x^H B x
                }
                    break;
                    
                default:
                    // to be implemented for the non-Hermitian problems
                    libmesh_error();
                    break;
            }
        }
    }
}




void
MAST::NonlinearSystem::
assemble_eigensystem () {
    
    libmesh_assert(_eigenproblem_assemble_system_object);
    
    _eigenproblem_assemble_system_object->eigenproblem_assemble(matrix_A,
                                                                matrix_B);
    
}



void
MAST::NonlinearSystem::
assemble_eigensystem_sensitivity(const libMesh::ParameterVector& parameters,
                                 const unsigned int p) {
    
    libmesh_assert(_eigenproblem_assemble_system_object);
    
    _eigenproblem_assemble_system_object->eigenproblem_sensitivity_assemble
    (parameters, p, matrix_A, matrix_B);
    
}




void
MAST::NonlinearSystem::
attach_eigenproblem_assemble_object(MAST::EigenSystemAssembly& assemble) {
    
    libmesh_assert(!_eigenproblem_assemble_system_object);
    
    _eigenproblem_assemble_system_object = &assemble;
}




void
MAST::NonlinearSystem::reset_eigenproblem_assemble_object () {
    
    _eigenproblem_assemble_system_object   = NULL;
}




void
MAST::NonlinearSystem::
initialize_condensed_dofs(MAST::PhysicsDisciplineBase& physics) {
    
    
    std::set<unsigned int> global_dirichlet_dofs_set;
    
    physics.get_system_dirichlet_bc_dofs(*this, global_dirichlet_dofs_set);
    
    
    // First, put all local dofs into non_dirichlet_dofs_set and
    std::set<unsigned int> local_non_condensed_dofs_set;
    for(unsigned int i=this->get_dof_map().first_dof();
        i<this->get_dof_map().end_dof();
        i++)
        local_non_condensed_dofs_set.insert(i);
    
    // Now erase the condensed dofs
    std::set<unsigned int>::iterator
    iter     = global_dirichlet_dofs_set.begin(),
    iter_end = global_dirichlet_dofs_set.end();
    
    for ( ; iter != iter_end ; ++iter) {
        
        unsigned int condensed_dof_index = *iter;
        
        if ( (this->get_dof_map().first_dof() <= condensed_dof_index) &&
            (condensed_dof_index < this->get_dof_map().end_dof()) )
            local_non_condensed_dofs_set.erase(condensed_dof_index);
    }
    
    // Finally, move local_non_condensed_dofs_set over to a vector for
    // convenience in solve()
    iter     = local_non_condensed_dofs_set.begin();
    iter_end = local_non_condensed_dofs_set.end();
    
    _local_non_condensed_dofs_vector.clear();
    
    for ( ; iter != iter_end; ++iter)
        _local_non_condensed_dofs_vector.push_back(*iter);
    
    _condensed_dofs_initialized = true;
}



