/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "solver/continuation_solver_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/assembly_elem_operation.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/linear_solver.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

// PETSc includes
#include <petscmat.h>


MAST::ContinuationSolverBase::ContinuationSolverBase():
max_it                    (10),
abs_tol                   (1.e-8),
rel_tol                   (1.e-8),
arc_length                (0.),
min_step                  (2.),
max_step                  (20.),
step_size_change_exponent (0.5),
step_desired_iters        (5),
schur_factorization       (true),
_initialized              (false),
_elem_ops                 (nullptr),
_assembly                 (nullptr),
_p                        (nullptr),
_p0                       (0.),
_X_scale                  (0.),
_p_scale                  (0.) {
    
}



MAST::ContinuationSolverBase::~ContinuationSolverBase() {
    
}



void
MAST::ContinuationSolverBase::
set_assembly_and_load_parameter(MAST::AssemblyElemOperations&   elem_ops,
                                MAST::AssemblyBase&             assembly,
                                MAST::Parameter&                p) {

    libmesh_assert(!_elem_ops);
    
    _elem_ops  =  &elem_ops;
    _assembly  =  &assembly;
    _p         =  &p;
    
    
}


void
MAST::ContinuationSolverBase::clear_assembly_and_load_parameters() {

    _elem_ops  =  nullptr;
    _assembly  =  nullptr;
    _p         =  nullptr;

}



void
MAST::ContinuationSolverBase::solve()  {
    
    libmesh_assert(_initialized);
    
    unsigned int
    iter    = 0;
    
    libMesh::NumericVector<Real>
    &X = *_assembly->system().solution;

    _p0     = (*_p)();
    _X0.reset(X.clone().release());
    

    // save data for possible reuse if the iterations are restarted.
    _save_iteration_data();

    Real
    norm0   = _res_norm(X, *_p),
    norm    = norm0;
    
    bool
    cont    = true;
    
    while (cont) {
        
        libMesh::out
        << std::setw(10) << "iter: "
        << std::setw(5)  << iter
        << std::setw(10) << "res-l2: "
        << std::setw(15) << norm
        << std::setw(20) << "relative res-l2: "
        << std::setw(15) << norm/norm0 << std::endl;
        
        _solve_NR_iterate(X, *_p);
        norm = _res_norm(X, *_p);
        iter++;
        
        if (norm < abs_tol)       cont = false;
        if (norm/norm0 < rel_tol) cont = false;
        if (iter >= max_it) {
            if (arc_length > min_step)   {
                
                // reduce step-size if possible, otherwise terminate
                libMesh::out
                << "   Retrying with smaller step-size" << std::endl;
                
                // reanalyze with a smaller step-size
                iter = 0;
                arc_length =  std::max(min_step, .5*arc_length);
                X.zero();
                X.add(1., *_X0);
                X.close();
                *_p = _p0;
                _reset_iterations();
                cont = true;
            }
            else
                cont = false;
        }
    }
        
    libMesh::out
    << std::setw(10) << "iter: "
    << std::setw(5)  << iter
    << std::setw(10) << "res-l2: "
    << std::setw(15) << norm
    << std::setw(20) << "relative res-l2: "
    << std::setw(15) << norm/norm0
    << std::setw(20) << "Terminated"  << std::endl;
    
    if (iter) {
        Real
        factor   = std::pow((1.*step_desired_iters)/(1.*iter+1.), step_size_change_exponent);
        if (factor > 1.) {
            arc_length = std::min(max_step, factor*arc_length);
            libMesh::out
            << std::setw(30) << "increased step-size: "
            << std::setw(15) << arc_length << std::endl;
        }
        else if (factor < 1.) {
            arc_length = std::max(min_step, factor*arc_length);
            libMesh::out
            << std::setw(30) << "reduced step-size: "
            << std::setw(15) << arc_length << std::endl;
        }
    }
    
}


void
MAST::ContinuationSolverBase::_solve(const libMesh::NumericVector<Real>  &X,
                                     const MAST::Parameter               &p,
                                     libMesh::NumericVector<Real>        &f,
                                     bool                                update_f,
                                     libMesh::NumericVector<Real>        &dfdp,
                                     bool                                update_dfdp,
                                     const libMesh::NumericVector<Real>  &dgdX,
                                     const Real                          dgdp,
                                     const Real                          g,
                                     libMesh::NumericVector<Real>        &dX,
                                     Real                                &dp) {
    
    libmesh_assert(_elem_ops);
    
    MAST::NonlinearSystem
    &system    = _assembly->system();
    
    libmesh_assert(system.operation() == MAST::NonlinearSystem::NONE);
    
    //
    //    {  f(X, p) }   =   { 0 }
    //    {  g(X, p) }   =   { 0 }
    //
    //     [df/dX    df/dp]  { dX } = { -f}
    //     [dg/dX    dg/dp]  { dp } = { -g}
    //
    
    // create matrix, vector and solver objects for the enlarged system
    // first, the sparsity pattern.
    // We will append one row and column to the sparsity pattern of the
    // system. The highest rank cpu will own this new unknown so that
    // the dof ids do not change

    /////////////////////////////////////////////////////////////////
    // create the new matrix and vector quantities for the
    // combined system
    /////////////////////////////////////////////////////////////////
    libMesh::DofMap& dof_map = system.get_dof_map();
    const libMesh::Parallel::Communicator
    &comm = system.comm();
    
    PetscInt
    rank      = comm.rank(),
    my_m      = dof_map.n_dofs(),
    my_n      = my_m,
    total_m   = my_m + 1,  // size of the system with the constraint equation
    total_n   = total_m,
    n_l       = dof_map.n_dofs_on_processor(system.processor_id()),
    m_l       = n_l,
    total_n_l = n_l,
    total_m_l = n_l;
    
    std::vector<libMesh::dof_id_type>
    n_nz       = dof_map.get_n_nz(),  // number of non-zeros per row in the diagonal block on this rank
    n_oz       = dof_map.get_n_oz();  // number of non-zeros per row in the off-diagonal blocks on other ranks

    
    if (rank == comm.size()-1) {

        // the last rank will own the new dof, so add one more non-zero to
        // the rows on this rank
        total_n_l++;
        total_m_l++;

        // add one more non-zero to the block matrix
        for (unsigned int i=0; i<n_nz.size(); i++)
            n_nz[i]++;
        // the final row on this rank will be a full row, that is all
        // entries of this row will be non-zero
        n_nz.push_back(total_n_l);
        // likewise, all entries in the off-diagonal block will be nonzero
        n_oz.push_back(total_n - total_n_l);
    }
    else {
        
        // add an extra non-zero in the off-diagonal block
        for (unsigned int i=0; i<n_oz.size(); i++)
            n_oz[i]++;
    }
    
    // create the matrix
    PetscErrorCode   ierr;
    Mat              mat;
    
    ierr = MatCreate(comm.get(), &mat);                      CHKERRABORT(comm.get(), ierr);
    ierr = MatSetSizes(mat,
                       total_m_l, total_n_l,
                       total_m, total_n);                             CHKERRABORT(comm.get(), ierr);
    
    if (libMesh::on_command_line("--solver_system_names")) {
        
        std::string nm = _assembly->system().name() + "_continuation_";
        MatSetOptionsPrefix(mat, nm.c_str());
    }
    ierr = MatSetFromOptions(mat);                                 CHKERRABORT(comm.get(), ierr);
    
    ierr = MatSeqAIJSetPreallocation(mat,
                                     my_m,
                                     (PetscInt*)&n_nz[0]);         CHKERRABORT(comm.get(), ierr);
    ierr = MatMPIAIJSetPreallocation(mat,
                                     0,
                                     (PetscInt*)&n_nz[0],
                                     0,
                                     (PetscInt*)&n_oz[0]);         CHKERRABORT(comm.get(), ierr);
    ierr = MatSeqBAIJSetPreallocation (mat, 1,
                                       0, (PetscInt*)&n_nz[0]);    CHKERRABORT(comm.get(), ierr);
    ierr = MatMPIBAIJSetPreallocation (mat, 1,
                                       0, (PetscInt*)&n_nz[0],
                                       0, (PetscInt*)&n_oz[0]);    CHKERRABORT(comm.get(), ierr);
    
    // now add the block entries
    ierr = MatSetOption(mat,
                        MAT_NEW_NONZERO_ALLOCATION_ERR,
                        PETSC_TRUE);                               CHKERRABORT(comm.get(), ierr);

    // now create the vectors
    Vec              res_vec, sol_vec;
    
    ierr = MatCreateVecs(mat, &res_vec, PETSC_NULL);               CHKERRABORT(comm.get(), ierr);
    ierr = MatCreateVecs(mat, &sol_vec, PETSC_NULL);               CHKERRABORT(comm.get(), ierr);
    
    
    std::unique_ptr<libMesh::SparseMatrix<Real> >
    jac_mat(new libMesh::PetscMatrix<Real>(mat, comm));
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    res(new libMesh::PetscVector<Real>(res_vec, comm)),
    sol(new libMesh::PetscVector<Real>(sol_vec, comm));

    /////////////////////////////////////////////////////////////////
    // now, assemble the information into the new matrix/vector
    /////////////////////////////////////////////////////////////////
    // compute the dfdp term if asked for, since that goes into the matrix.
    _assembly->set_elem_operation_object(*_elem_ops);
    system.set_operation(MAST::NonlinearSystem::FORWARD_SENSITIVITY_SOLVE);
    if (update_dfdp)
        _assembly->sensitivity_assemble(*system.solution, true, p, dfdp);
    
    // now compute the jacobian
    system.set_operation(MAST::NonlinearSystem::NONLINEAR_SOLVE);
    _assembly->close_matrix = false;
    _assembly->residual_and_jacobian(X,
                                     update_f?res.get():nullptr,
                                     jac_mat.get(),
                                     system);
    _assembly->close_matrix = true;
    
    _assembly->clear_elem_operation_object();
    system.set_operation(MAST::NonlinearSystem::NONE);
    
    // first set the data for the last column
    if (!update_f) {
        for (unsigned int i=dof_map.first_dof(rank); i<dof_map.end_dof(rank); i++)
            res->set(i, f.el(i));
    }
    
    // diagonal entry on the last rank
    if (rank == comm.size()-1)
        res->set(total_m-1, g);
    
    // finish assembling the residual vector if it was not updated by
    // the residual and jacobian routine
    res->close();
    res->scale(-1.);
    res->close();


    // first set the data for the last column
    for (unsigned int i=dof_map.first_dof(rank);
         i<dof_map.end_dof(rank); i++) {
        
        jac_mat->set(        i, total_m-1,   dfdp.el(i));
        jac_mat->set(total_m-1,         i,   dgdX.el(i));
    }
    
    // diagonal entry on the last rank
    if (rank == comm.size()-1)
        jac_mat->set(total_m-1, total_m-1, dgdp);

    jac_mat->close();
    
    
    /////////////////////////////////////////////////////////////////
    // now, copy the information to the new matrix/vector
    /////////////////////////////////////////////////////////////////
    KSP        ksp;
    PC         pc;
    
    // setup the KSP
    ierr = KSPCreate(comm.get(), &ksp);       CHKERRABORT(comm.get(), ierr);
    
    if (libMesh::on_command_line("--solver_system_names")) {
        
        std::string nm = _assembly->system().name() + "_continuation_";
        KSPSetOptionsPrefix(ksp, nm.c_str());
    }

    ierr = KSPSetOperators(ksp, mat, mat);    CHKERRABORT(comm.get(), ierr);
    ierr = KSPSetFromOptions(ksp);            CHKERRABORT(comm.get(), ierr);
    
    // setup the PC
    ierr = KSPGetPC(ksp, &pc);                CHKERRABORT(comm.get(), ierr);
    ierr = PCSetFromOptions(pc);              CHKERRABORT(comm.get(), ierr);
    
    // now solve
    ierr = KSPSolve(ksp, res_vec, sol_vec);

    // copy the solution back to the system
    dX.zero();
    for (unsigned int i=dof_map.first_dof(rank);
         i<dof_map.end_dof(rank); i++)
        dX.set(i, sol->el(i));
    dX.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
    system.get_dof_map().enforce_constraints_exactly (system,
                                                      &dX,
                                                      /* homogeneous = */ true);
#endif
    
    std::vector<Real> val(1, 0.);
    std::vector<libMesh::numeric_index_type> dofs(1, total_m-1);
    sol->localize(val, dofs);
    dp = val[0];
    
    // destroy the objects
    ierr = KSPDestroy(&ksp);
    ierr = MatDestroy(&mat);
    ierr = VecDestroy(&res_vec);
    ierr = VecDestroy(&sol_vec);
}




void
MAST::ContinuationSolverBase::
_solve_schur_factorization(const libMesh::NumericVector<Real>  &X,
                           const MAST::Parameter               &p,
                           libMesh::SparseMatrix<Real>         &jac,
                           bool                                update_jac,
                           libMesh::NumericVector<Real>        &f,
                           bool                                update_f,
                           libMesh::NumericVector<Real>        &dfdp,
                           bool                                update_dfdp,
                           libMesh::NumericVector<Real>        &dXdp,
                           bool                                update_dXdp,
                           const libMesh::NumericVector<Real>  &dgdX,
                           const Real                          dgdp,
                           const Real                          g,
                           libMesh::NumericVector<Real>        &dX,
                           Real                                &dp) {
    
    libmesh_assert(_elem_ops);
    
    MAST::NonlinearSystem
    &system    = _assembly->system();
    
    libmesh_assert(system.operation() == MAST::NonlinearSystem::NONE);
    
    //
    //    {  f(X, p) }   =   { 0 }
    //    {  g(X, p) }   =   { 0 }
    //
    //     [df/dX    df/dp]  { dX } = { -f}
    //     [dg/dX    dg/dp]  { dp } = { -g}
    //
    //   Schur factorization:
    //     df/dX  dX =    -f - df/dp dp
    //
    //   Substitute in second equation
    //     dg/dp  dp - dg/dX inv(df/dX) ( f + df/dp dp) = -g
    // =>  [dg/dp - dg/dX inv(df/dX) df/dp ] dp = -g + dg/dX inv(df/dX) f
    // =>  [dg/dp + dg/dX dX/dp ] dp = -g + dg/dX r1
    //
    //   1.  solve   r1        =  inv(df/dX) f
    //   2.  solve   dXdp      = -inv(df/dX) df/dp
    //   3.  solve   a         = dg/dp + dg/dX dXdp
    //   4.  solve   a  dp     = -g + dg/dX r1
    //   5.  solve   df/dX  dX =    -f - df/dp dp
    //
    
    
    std::pair<unsigned int, Real>
    solver_params = system.get_linear_solve_parameters(),
    rval;
    
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    r1(X.zero_clone().release());
    
    libMesh::SparseMatrix<Real>
    *pc  = system.request_matrix("Preconditioner");
    
    Real
    a    = 0.;
    
    //////////////////////////////////////////////////////////
    //         STEP 1:  r1  = inv(df/dX) f
    //////////////////////////////////////////////////////////
    system.set_operation(MAST::NonlinearSystem::NONLINEAR_SOLVE);
    _assembly->set_elem_operation_object(*_elem_ops);
    
    if (update_f || update_jac)
        _assembly->residual_and_jacobian(X,
                                         update_f?     &f:nullptr,
                                         update_jac? &jac:nullptr,
                                         system);
    rval = system.linear_solver->solve (jac, pc,
                                        *r1,
                                        f,
                                        solver_params.second,
                                        solver_params.first);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
    system.get_dof_map().enforce_constraints_exactly (system,
                                                      r1.get(),
                                                      /* homogeneous = */ true);
#endif


    //////////////////////////////////////////////////////////
    //         STEP 2:  dXdp  = - inv(df/dX) df/dp
    //////////////////////////////////////////////////////////
    system.set_operation(MAST::NonlinearSystem::FORWARD_SENSITIVITY_SOLVE);
    if (update_dfdp)
        _assembly->sensitivity_assemble(*system.solution, true, p, dfdp);

    if (update_dfdp || update_dXdp) {
        
        rval = system.linear_solver->solve (jac, pc,
                                            dXdp,
                                            dfdp,
                                            solver_params.second,
                                            solver_params.first);
        
        dXdp.scale(-1.);
        dXdp.close();
        
        // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        system.get_dof_map().enforce_constraints_exactly (system,
                                                          &dXdp,
                                                          /* homogeneous = */ true);
#endif
    }

    //////////////////////////////////////////////////////////
    //         STEP 3:  a   = dg/dp + dg/dX dXdp
    //////////////////////////////////////////////////////////
    a   = dgdp + dgdX.dot(dXdp);
    
    //////////////////////////////////////////////////////////
    //         STEP 4:  a  dp     = -g + dg/dX r1
    //////////////////////////////////////////////////////////
    dp  = 1./a * (- g + dgdX.dot(*r1));
    
    //////////////////////////////////////////////////////////
    //         STEP 5:  df/dX  dX =    -f - df/dp dp
    //////////////////////////////////////////////////////////
    r1->zero();
    r1->add(1., f);
    r1->add(dp, dfdp);
    r1->scale(-1.);
    r1->close();
    rval = system.linear_solver->solve (jac, pc,
                                        dX,
                                        *r1,
                                        solver_params.second,
                                        solver_params.first);
    
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    system.get_dof_map().enforce_constraints_exactly (system,
                                                      &dX,
                                                      /* homogeneous = */ true);
#endif

    _assembly->clear_elem_operation_object();
    system.set_operation(MAST::NonlinearSystem::NONE);
}


Real
MAST::ContinuationSolverBase::
_res_norm(const libMesh::NumericVector<Real>                  &X,
          const MAST::Parameter                               &p) {
    
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem
    &system    = _assembly->system();

    libmesh_assert(system.operation() == MAST::NonlinearSystem::NONE);

    Real
    val = 0.;

    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release());

    system.set_operation(MAST::NonlinearSystem::NONLINEAR_SOLVE);
    _assembly->set_elem_operation_object(*_elem_ops);
    
    _assembly->residual_and_jacobian(X, f.get(), nullptr, system);

    _assembly->clear_elem_operation_object();
    system.set_operation(MAST::NonlinearSystem::NONE);

    val = std::pow(std::pow(_g(X, p), 2) + std::pow(f->l2_norm(), 2), 0.5);
    
    return val;
}

