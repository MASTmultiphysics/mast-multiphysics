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
#include "solver/multiphysics_nonlinear_solver.h"
#include "base/transient_assembly.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"


//---------------------------------------------------------------
// method for matrix vector multiplicaiton of the off-diagonal component y=Ax
struct
__mast_multiphysics_petsc_shell_context {
    unsigned int                                i;
    unsigned int                                j;
    SNES                                     snes;
    MAST::MultiphysicsNonlinearSolverBase* solver;
};


PetscErrorCode
__mast_multiphysics_petsc_mat_mult(Mat mat,Vec dx,Vec y) {

    LOG_SCOPE("mat_mult()", "PetscNonlinearSolver");
    
    PetscErrorCode ierr=0;
    
    libmesh_assert(mat);
    libmesh_assert(dx);
    libmesh_assert(y);
    
    
    void * ctx = PETSC_NULL;
    
    // get the matrix context. This is for the i^th row-block and j^th column
    // block, which uses the linear perturbation in the j^th block.
    ierr = MatShellGetContext(mat, &ctx);
    
    __mast_multiphysics_petsc_shell_context
    *mat_ctx = static_cast<__mast_multiphysics_petsc_shell_context*> (ctx);


    MAST::MultiphysicsNonlinearSolverBase
    *solver = mat_ctx->solver;
    
    // get the current nonlinear solution
    Vec x;
    
    ierr = SNESGetSolution(mat_ctx->snes, &x);              CHKERRABORT(solver->comm().get(), ierr);

    
    const unsigned int
    nd = solver->n_disciplines();
    
    std::vector<Vec>
    sol (nd);
    
    std::vector<libMesh::NumericVector<Real>*>
    sys_sols (nd, nullptr),
    sys_dsols(nd, nullptr);
    
    
    //////////////////////////////////////////////////////////////////
    // get the subvectors for each discipline
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        
        // extract the subvector for this system
        ierr = VecGetSubVector( x, sys_is,  &sol[i]);      CHKERRABORT(solver->comm().get(), ierr);

        // use the nonlinear sol of all disciplines, by the perturbed sol of only
        // one should be nonzero.
        sys_sols[i]   = new libMesh::PetscVector<Real>( sol[i], sys.comm());
        if (i == mat_ctx->j) {
            sys_dsols[i]  = new libMesh::PetscVector<Real>(dx, sys.comm());
            
            // Enforce constraints (if any) exactly on the
            // current_local_solution.  This is the solution vector that is
            // actually used in the computation of the residual below, and is
            // not locked by debug-enabled PETSc the way that "x" is.
            sys.get_dof_map().enforce_constraints_exactly(sys, sys_dsols[i],
                                                          true /* homogeneous = true */);
        }
        else
            sys_dsols[i]  = sys_sols[i]->zero_clone().release();
    }
    
    //////////////////////////////////////////////////////////////////
    // initialize the data structures before calculation of residuals
    //////////////////////////////////////////////////////////////////
    if (solver->get_pre_residual_update_object())
        solver->get_pre_residual_update_object()->update_at_perturbed_solution(sys_sols,
                                                                               sys_dsols);
    
    //////////////////////////////////////////////////////////////////
    // calculate the matrix-vector product
    //////////////////////////////////////////////////////////////////
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    res(new libMesh::PetscVector<Real>(y, solver->comm()));
    
    // system for this discipline
    MAST::NonlinearSystem& sys = solver->get_system_assembly(mat_ctx->i).system();

    solver->get_system_assembly(mat_ctx->i).linearized_jacobian_solution_product
    (*sys_sols [mat_ctx->i],
     *sys_dsols[mat_ctx->i],
     *res,
     sys);
    
    res->close();
    
    //////////////////////////////////////////////////////////////////
    // resotre the subvectors
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        
        // delete the NumericVector wrappers
        delete sys_sols[i];
        delete sys_dsols[i];
        
        // now restore the subvectors
        ierr = VecRestoreSubVector(x, sys_is, &sol[i]);  CHKERRABORT(solver->comm().get(), ierr);
    }
    
    return ierr;
}


//---------------------------------------------------------------
// this function is called by PETSc to evaluate the residual at X
PetscErrorCode
__mast_multiphysics_petsc_snes_residual (SNES snes, Vec x, Vec r, void * ctx) {
    
    LOG_SCOPE("residual()", "PetscMultiphysicsNonlinearSolver");
    
    PetscErrorCode ierr=0;
    
    libmesh_assert(x);
    libmesh_assert(r);
    libmesh_assert(ctx);
    
    MAST::MultiphysicsNonlinearSolverBase * solver =
    static_cast<MAST::MultiphysicsNonlinearSolverBase*> (ctx);

    const unsigned int
    nd = solver->n_disciplines();
    
    std::vector<Vec>
    sol(nd),
    res(nd);

    std::vector<libMesh::NumericVector<Real>*>
    sys_sols(nd, nullptr),
    sys_res (nd, nullptr);
    

    //////////////////////////////////////////////////////////////////
    // get the subvectors for each discipline
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        
        // extract the subvector for this system
        ierr = VecGetSubVector(x, sys_is, &sol[i]);      CHKERRABORT(solver->comm().get(), ierr);
        ierr = VecGetSubVector(r, sys_is, &res[i]);      CHKERRABORT(solver->comm().get(), ierr);
        
        sys_sols[i] = new libMesh::PetscVector<Real>(sol[i], sys.comm()),
        sys_res[i]  = new libMesh::PetscVector<Real>(res[i], sys.comm());

        // Enforce constraints (if any) exactly on the
        // current_local_solution.  This is the solution vector that is
        // actually used in the computation of the residual below, and is
        // not locked by debug-enabled PETSc the way that "x" is.
        sys.get_dof_map().enforce_constraints_exactly(sys, sys_sols[i]);
    }

    //////////////////////////////////////////////////////////////////
    // initialize the data structures before calculation of residuals
    //////////////////////////////////////////////////////////////////
    if (solver->get_pre_residual_update_object())
        solver->get_pre_residual_update_object()->update_at_solution(sys_sols);
    
    //////////////////////////////////////////////////////////////////
    // calculate the residuals
    //////////////////////////////////////////////////////////////////
    std::vector<Real>
    l2_norms(nd, 0.);
    
    Real
    global_l2 = 0.;

    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        solver->get_system_assembly(i).residual_and_jacobian (*sys_sols[i],
                                                              sys_res[i],
                                                              nullptr,
                                                              sys);
        
        sys_res[i]->close();
        
        // now calculate the norms
        l2_norms[i]  = sys_res[i]->l2_norm();
        global_l2   += pow(l2_norms[i], 2);
    }
    
    global_l2 = pow(global_l2, 0.5);
    
    // write the global and component-wise l2 norms of the residual
    libMesh::out
    << "|| R ||_2 = " << global_l2 << "  : || R_i ||_2 = ( ";
    for (unsigned int i=0; i<nd; i++) {
        libMesh::out << l2_norms[i];
        if (i < nd-1)
            libMesh::out << " , ";
    }
    libMesh::out << " )" << std::endl;
    
    //////////////////////////////////////////////////////////////////
    // resotre the subvectors
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
                
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];

        // delete the NumericVector wrappers
        delete sys_sols[i];
        delete sys_res[i];
        
        // now restore the subvectors
        ierr = VecRestoreSubVector(x, sys_is, &sol[i]);  CHKERRABORT(solver->comm().get(), ierr);
        ierr = VecRestoreSubVector(r, sys_is, &res[i]);  CHKERRABORT(solver->comm().get(), ierr);
    }
    
    return ierr;
}



//---------------------------------------------------------------
// this function is called by PETSc to evaluate the Jacobian at X
PetscErrorCode
__mast_multiphysics_petsc_snes_jacobian(SNES snes, Vec x, Mat jac, Mat pc, void * ctx)
{
    LOG_SCOPE("jacobian()", "PetscMultiphysicsNonlinearSolver");
    
    PetscErrorCode ierr=0;
    
    libmesh_assert(x);
    libmesh_assert(jac);
    libmesh_assert(pc);
    libmesh_assert(ctx);

    
    MAST::MultiphysicsNonlinearSolverBase * solver =
    static_cast<MAST::MultiphysicsNonlinearSolverBase*> (ctx);
    
    const unsigned int
    nd = solver->n_disciplines();
    
    std::vector<Vec>
    sol(nd);
    
    std::vector<libMesh::NumericVector<Real>*>
    sys_sols(nd, nullptr);
    
    
    //////////////////////////////////////////////////////////////////
    // get the subvectors for each discipline
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        
        // extract the subvector for this system
        ierr = VecGetSubVector(x, sys_is, &sol[i]);      CHKERRABORT(solver->comm().get(), ierr);
        
        sys_sols[i] = new libMesh::PetscVector<Real>(sol[i], sys.comm());

        
        // Enforce constraints (if any) exactly on the
        // current_local_solution.  This is the solution vector that is
        // actually used in the computation of the residual below, and is
        // not locked by debug-enabled PETSc the way that "x" is.
        sys.get_dof_map().enforce_constraints_exactly(sys, sys_sols[i]);
    }
    
    
    //////////////////////////////////////////////////////////////////
    // initialize the data structures before calculation of residuals
    //////////////////////////////////////////////////////////////////
    if (solver->get_pre_residual_update_object())
        solver->get_pre_residual_update_object()->update_at_solution(sys_sols);
    
    //////////////////////////////////////////////////////////////////
    // calculate the residuals
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys = solver->get_system_assembly(i).system();
        
        //PetscMatrix<Number> PC(pc, sys.comm());
        //PC.attach_dof_map(sys.get_dof_map());
        //PC.close();

        solver->get_system_assembly(i).residual_and_jacobian (*sys_sols[i],
                                                              nullptr,
                                                              sys.matrix,
                                                              sys);
        
        sys.matrix->close();
    }

    //////////////////////////////////////////////////////////////////
    // resotre the subvectors
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i< nd; i++) {
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        
        // delete the NumericVector wrappers
        delete sys_sols[i];
        
        // now restore the subvectors
        ierr = VecRestoreSubVector(x, sys_is, &sol[i]);  CHKERRABORT(solver->comm().get(), ierr);
    }

    
    // call assembly of the global matrix
    ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);  CHKERRABORT(solver->comm().get(), ierr);
    ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);    CHKERRABORT(solver->comm().get(), ierr);
    
    
    return ierr;
}





MAST::MultiphysicsNonlinearSolverBase::
MultiphysicsNonlinearSolverBase(const libMesh::Parallel::Communicator& comm_in,
                                const std::string& nm,
                                unsigned int n):
libMesh::ParallelObject       (comm_in),
_name                         (nm),
_n_disciplines                (n),
_update                       (nullptr),
_discipline_assembly          (n, nullptr),
_is                           (_n_disciplines, PETSC_NULL),
_sub_mats                     (_n_disciplines*_n_disciplines, PETSC_NULL),
_n_dofs                       (0) {
    
}





MAST::MultiphysicsNonlinearSolverBase::~MultiphysicsNonlinearSolverBase() {
    
}




void
MAST::MultiphysicsNonlinearSolverBase::
set_system_assembly(unsigned int i,
                    MAST::TransientAssembly& assembly) {
    
    // make sure that the index is within bounds
    libmesh_assert_less(i, _n_disciplines);
    
    // also make sure that the specific discipline has not been set already
    libmesh_assert(!_discipline_assembly[i]);
    
    _discipline_assembly[i] = &assembly;
}



MAST::TransientAssembly&
MAST::MultiphysicsNonlinearSolverBase::
get_system_assembly(unsigned int i) {
    
    // make sure that the index is within bounds
    libmesh_assert_less(i, _n_disciplines);
    
    // also make sure that the specific discipline has not been set already
    libmesh_assert(_discipline_assembly[i]);
    
    return *_discipline_assembly[i];
}



void
MAST::MultiphysicsNonlinearSolverBase::solve() {
    
    
    // make sure that all systems have been specified
    bool p = true;
    for (unsigned int i=0; i<_n_disciplines; i++)
        p = ((_discipline_assembly[i] != nullptr) && p);
    libmesh_assert(p);
    
    
    //////////////////////////////////////////////////////////////////////
    // create the solver context. This will be partially initialized now
    // since it is needed in the shell matrix context
    //////////////////////////////////////////////////////////////////////
    PetscErrorCode   ierr;
    SNES             snes;

    // setup the solver context
    ierr = SNESCreate(this->comm().get(), &snes);      CHKERRABORT(this->comm().get(), ierr);

    //////////////////////////////////////////////////////////////////////
    // create a petsc nested matrix of size NxN
    //////////////////////////////////////////////////////////////////////
    const bool       sys_name = libMesh::on_command_line("--solver_system_names");
    std::string      nm;
    unsigned int     n_local_dofs = 0;
    std::vector<__mast_multiphysics_petsc_shell_context>
    mat_ctx(_n_disciplines*_n_disciplines);
    
    
    // all diagonal blocks use the system matrcices, while shell matrices
    // are created for the off-diagonal terms.
    _n_dofs = 0.;
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        MAST::NonlinearSystem& sys = _discipline_assembly[i]->system();
        
        // add the number of dofs in this system to the global count
        _n_dofs      +=  sys.n_dofs();
        n_local_dofs +=  sys.n_local_dofs();
        
        for (unsigned int j=0; j<_n_disciplines; j++) {
            
            if (i==j) {
                
                // the diagonal matrix
                _sub_mats[i*_n_disciplines+j] =
                dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix)->mat();
            }
            else {
                
                MAST::NonlinearSystem& sys_j = _discipline_assembly[j]->system();
                
                PetscInt
                mat_i_m    = sys.get_dof_map().n_dofs(),
                mat_i_n_l  = sys.get_dof_map().n_dofs_on_processor(sys.processor_id()),
                mat_j_m    = sys_j.get_dof_map().n_dofs(),
                mat_j_n_l  = sys_j.get_dof_map().n_dofs_on_processor(sys.processor_id());
                
                // the off-diagonal matrix
                ierr = MatCreateShell(this->comm().get(),
                                      mat_i_n_l,
                                      mat_j_n_l,
                                      mat_i_m,
                                      mat_j_m,
                                      PETSC_NULL,
                                      &_sub_mats[i*_n_disciplines+j]);
                CHKERRABORT(this->comm().get(), ierr);
                
                // initialize the context and tell the matrix about it
                mat_ctx[i*_n_disciplines+j].i      = i;
                mat_ctx[i*_n_disciplines+j].j      = j;
                mat_ctx[i*_n_disciplines+j].snes   = snes;
                mat_ctx[i*_n_disciplines+j].solver = this;
                
                ierr = MatShellSetContext(_sub_mats[i*_n_disciplines+j],
                                          &mat_ctx[i*_n_disciplines+j]);
                CHKERRABORT(this->comm().get(), ierr);
                
                // set the mat-vec multiplication operation for this
                // shell matrix
                ierr = MatShellSetOperation(_sub_mats[i*_n_disciplines+j],
                                            MATOP_MULT,
                                            (void(*)(void))__mast_multiphysics_petsc_mat_mult);
                CHKERRABORT(this->comm().get(), ierr);
            }
        }
    }
    
    
    ierr = MatCreateNest(this->comm().get(),
                         _n_disciplines, PETSC_NULL,
                         _n_disciplines, PETSC_NULL,
                         &_sub_mats[0],
                         &_mat);
    CHKERRABORT(this->comm().get(), ierr);

    // we need to turn off MULT_TRANSPOSE operator for the global matrix
    // since that is not implemented for the shell matrices pushes PETSc 3.7.4
    // to produce an error about it.
    ierr = MatShellSetOperation(_mat,
                                MATOP_MULT_TRANSPOSE,
                                PETSC_NULL);
    CHKERRABORT(this->comm().get(), ierr);

    if (sys_name) {
        
        nm = this->name() + "_";
        MatSetOptionsPrefix(_mat, nm.c_str());
    }
    ierr = MatSetFromOptions(_mat);
    CHKERRABORT(this->comm().get(), ierr);
    
    // get the IS belonging to each block
    ierr  =    MatNestGetISs(_mat, &_is[0], PETSC_NULL);
    CHKERRABORT(this->comm().get(), ierr);
    
    //////////////////////////////////////////////////////////////////////
    // setup the vector for solution
    //////////////////////////////////////////////////////////////////////
    ierr = VecCreate(this->comm().get(), &_sol);        CHKERRABORT(this->comm().get(), ierr);
    ierr = VecSetSizes(_sol, n_local_dofs, _n_dofs);    CHKERRABORT(this->comm().get(), ierr);
    ierr = VecSetType(_sol, VECMPI);                    CHKERRABORT(this->comm().get(), ierr);
    ierr = VecDuplicate(_sol, &_res);                   CHKERRABORT(this->comm().get(), ierr);
    
    ierr = MatShellSetOperation(_mat,
                                MATOP_MULT_TRANSPOSE_ADD,
                                PETSC_NULL);
    CHKERRABORT(this->comm().get(), ierr);
    ierr = MatShellSetOperation(_mat,
                                MATOP_TRANSPOSE,
                                PETSC_NULL);
    CHKERRABORT(this->comm().get(), ierr);
    
    //////////////////////////////////////////////////////////////////////
    // now initialize the vector from the system solutions, which will serve
    // as the initial solution for the coupled system
    //////////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i<_n_disciplines; i++) {

        // solution from the provided system
        libMesh::NumericVector<Real>
        &sys_sol = *_discipline_assembly[i]->system().solution;
        
        // limiting indices for the system, and multiphysics assembly. Should
        // be the same
        int
        multiphysics_first = 0,
        multiphysics_last  = 0,
        first              = sys_sol.first_local_index(),
        last               = sys_sol.last_local_index();
        
        // get the subvector corresponding to the the discipline
        Vec sub_vec;

        ierr = VecGetSubVector(_sol, _is[i], &sub_vec);  CHKERRABORT(this->comm().get(), ierr);
        ierr = VecGetOwnershipRange(sub_vec,
                                    &multiphysics_first,
                                    &multiphysics_last); CHKERRABORT(this->comm().get(), ierr);
        
        // the first and last indices must match between the two representations
        libmesh_assert_equal_to(multiphysics_first, first);
        libmesh_assert_equal_to( multiphysics_last,  last);
        
        std::unique_ptr<libMesh::NumericVector<Real> >
        multiphysics_sol(new libMesh::PetscVector<Real>(sub_vec, this->comm()));
        
        for (unsigned int i=first; i<last; i++)
            multiphysics_sol->set(i, sys_sol(i));

        multiphysics_sol->close();
        
        ierr = VecRestoreSubVector(_sol, _is[i], &sub_vec);  CHKERRABORT(this->comm().get(), ierr);
    }
    
    
    
    
    //////////////////////////////////////////////////////////////////////
    // initialize the solver context
    //////////////////////////////////////////////////////////////////////
    KSP            ksp;
    PC             pc;
    
    
    
    // tell the solver where to store the Jacobian and how to calculate it
    ierr = SNESSetFunction (snes,
                            _res,
                            __mast_multiphysics_petsc_snes_residual,
                            this);
    ierr = SNESSetJacobian(snes,
                           _mat,
                           _mat,
                           __mast_multiphysics_petsc_snes_jacobian,
                           this);
    
    
    if (sys_name) {
        
        nm = this->name() + "_";
        SNESSetOptionsPrefix(snes, nm.c_str());
    }
    
    
    
    // setup the ksp
    ierr = SNESGetKSP (snes, &ksp);                   CHKERRABORT(this->comm().get(), ierr);
    ierr = SNESSetFromOptions(snes);                  CHKERRABORT(this->comm().get(), ierr);
    
    
    
    // setup the pc
    ierr = KSPGetPC(ksp, &pc);                        CHKERRABORT(this->comm().get(), ierr);
    
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        if (sys_name) {
            
            nm = _discipline_assembly[i]->system().name();
            ierr = PCFieldSplitSetIS(pc, nm.c_str(), _is[i]); CHKERRABORT(this->comm().get(), ierr);
        }
        else
            ierr = PCFieldSplitSetIS(pc, nullptr, _is[i]);CHKERRABORT(this->comm().get(), ierr);
    }
    

    //ierr = SNESSetSolution(snes, _sol);
    //this->verify_gateaux_derivatives(snes);
    
    //////////////////////////////////////////////////////////////////////
    // now, solve
    //////////////////////////////////////////////////////////////////////
    START_LOG("SNESSolve", this->name()+"_MultiphysicsSolve");
    
    // now solve
    ierr = SNESSolve(snes, PETSC_NULL, _sol);
    
    STOP_LOG("SNESSolve", this->name()+"_MultiphysicsSolve");
    
    
    //////////////////////////////////////////////////////////////////////
    // now copy the solution back to the system solution vector
    //////////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        // solution from the provided system
        libMesh::NumericVector<Real>
        &sys_sol = *_discipline_assembly[i]->system().solution;
        
        // limiting indices for the system, and multiphysics assembly. Should
        // be the same
        int
        multiphysics_first = 0,
        multiphysics_last  = 0,
        first              = sys_sol.first_local_index(),
        last               = sys_sol.last_local_index();
        
        // get the subvector corresponding to the the discipline
        Vec sub_vec;
        
        ierr = VecGetSubVector(_sol, _is[i], &sub_vec);  CHKERRABORT(this->comm().get(), ierr);
        ierr = VecGetOwnershipRange(sub_vec,
                                    &multiphysics_first,
                                    &multiphysics_last); CHKERRABORT(this->comm().get(), ierr);
        
        // the first and last indices must match between the two representations
        libmesh_assert_equal_to(multiphysics_first, first);
        libmesh_assert_equal_to( multiphysics_last,  last);
        
        std::unique_ptr<libMesh::NumericVector<Real> >
        multiphysics_sol(new libMesh::PetscVector<Real>(sub_vec, this->comm()));
        
        for (unsigned int i=first; i<last; i++)
            sys_sol.set(i, (*multiphysics_sol)(i));

        ierr = VecRestoreSubVector(_sol, _is[i], &sub_vec);  CHKERRABORT(this->comm().get(), ierr);
    }

    
    
    // destroy the Petsc contexts
    ierr = SNESDestroy(&snes);                        CHKERRABORT(this->comm().get(), ierr);
    for (unsigned int i=0; i<_n_disciplines; i++)
        for (unsigned int j=0; j<_n_disciplines; j++)
            if (i != j) {
                ierr = MatDestroy(&_sub_mats[i*_n_disciplines+j]);
                CHKERRABORT(this->comm().get(), ierr);
            }
    
    ierr = MatDestroy(&_mat);                          CHKERRABORT(this->comm().get(), ierr);
    ierr = VecDestroy(&_sol);                          CHKERRABORT(this->comm().get(), ierr);
    ierr = VecDestroy(&_res);                          CHKERRABORT(this->comm().get(), ierr);
    
}




void
MAST::MultiphysicsNonlinearSolverBase::verify_gateaux_derivatives(SNES snes) {
    
    PetscErrorCode ierr=0;

    // create vectors for the sol, dsol, and res of the global and
    // disciplinary systems
    std::unique_ptr<libMesh::NumericVector<Real> >
    global_sol  (new libMesh::PetscVector<Real>(_sol, this->comm())),
    global_res  (new libMesh::PetscVector<Real>(_res, this->comm())),
    global_res0 (global_sol->zero_clone().release());
    
    // first calculate the baseline residual for the global solution vector
    __mast_multiphysics_petsc_snes_residual(snes,
                                            _sol,
                                            dynamic_cast<libMesh::PetscVector<Real>*>(global_res0.get())->vec(),
                                            this);
    
    // value of perturbation used
    const Real
    delta = 1.0e-4;
    
    // now, iterate over each dof on the disciplines, and calculate the
    // Gateaux derivative
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        // system for this discipline
        MAST::NonlinearSystem& sys_i = this->get_system_assembly(i).system();
        
        // get the IS for this system
        IS sys_is_i = this->index_sets()[i];
        
        // verify the derivative wrt dofs from other disciplines
        for (unsigned int j=0; j<_n_disciplines; j++) {

            if (i != j) {
                
                IS sys_is_j = this->index_sets()[j];
                
                // system for this discipline
                MAST::NonlinearSystem& sys_j = this->get_system_assembly(j).system();
                
                std::unique_ptr<libMesh::NumericVector<Real> >
                dsol_j     (sys_j.solution->zero_clone().release()),
                dJac_ij_dXj(sys_i.solution->zero_clone().release());
                
                for (unsigned int k=0; k<sys_j.n_dofs(); k++) {
                    
                    dsol_j->zero();
                    dsol_j->add(k, delta);
                    dsol_j->close();
                    
                    Vec
                    vecx = dynamic_cast<libMesh::PetscVector<Real>*>(dsol_j.get())->vec(),
                    vecy = dynamic_cast<libMesh::PetscVector<Real>*>(dJac_ij_dXj.get())->vec();
                    
                    __mast_multiphysics_petsc_mat_mult(_sub_mats[i*_n_disciplines+j], vecx, vecy);

                    // now perform the same calculation with the finite differencing
                    // using residual calculation
                    // extract the subvector for this system
                    Vec
                    res_i,
                    sol_j;
                    ierr = VecGetSubVector(_sol, sys_is_j, &sol_j);      CHKERRABORT(this->comm().get(), ierr);
                    
                    Real
                    old_val = 0.;
                    
                    // create a vector for modification of the ith value of this solution
                    {
                        libMesh::PetscVector<Real> v(sol_j, this->comm());
                        
                        // perturb the solution vector
                        old_val = v.el(k);
                        v.add(k, delta);
                        v.close();
                    }
                    
                    ierr = VecRestoreSubVector( _sol, sys_is_j, &sol_j);      CHKERRABORT(this->comm().get(), ierr);
                    
                    // now calculate the residual
                    __mast_multiphysics_petsc_snes_residual(snes,
                                                            _sol,
                                                            _res,
                                                            this);

                    // reset the global solution vector
                    ierr = VecGetSubVector(_sol, sys_is_j, &sol_j);      CHKERRABORT(this->comm().get(), ierr);
                    
                    // create a vector for modification of the ith value of this solution
                    {
                        libMesh::PetscVector<Real> v(sol_j, this->comm());
                        
                        // perturb the solution vector
                        v.set(k, old_val);
                        v.close();
                    }
                    
                    ierr = VecRestoreSubVector( _sol, sys_is_j, &sol_j);      CHKERRABORT(this->comm().get(), ierr);
                    
                    global_res->close();
                    global_res->add(-1., *global_res0);
                    global_res->close();
                    

                    // get the finite differenced residual
                    ierr = VecGetSubVector(_res, sys_is_i, &res_i);      CHKERRABORT(this->comm().get(), ierr);
                    
                    // create a vector for comparison
                    {
                        libMesh::PetscVector<Real> v(res_i, this->comm());
                        for (unsigned int l=0; l<sys_i.n_dofs(); l++) {
                            
                            if (dJac_ij_dXj->el(l) != v.el(l))
                                libMesh::out
                                << k << "   "
                                << l << "   "
                                << dJac_ij_dXj->el(l) << "   "
                                << v.el(l) << std::endl;
                        }
                        
                        libMesh::out << std::endl;
                    }
                    
                    ierr = VecRestoreSubVector( _res, sys_is_i, &res_i);      CHKERRABORT(this->comm().get(), ierr);

                }
            }
        }
    }
}



