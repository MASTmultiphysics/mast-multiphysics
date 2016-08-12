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


// MAST includes
#include "solver/multiphysics_nonlinear_solver.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"


//---------------------------------------------------------------
// method for matrix vector multiplicaiton of the off-diagonal component y=Ax
struct
__mast_multiphysics_petsc_shell_context {
    unsigned int i;
    unsigned int j;
    MAST::MultiphysicsNonlinearSolverBase* solver;
};


PetscErrorCode
__mast_multiphysics_petsc_mat_mult(Mat mat,Vec x,Vec y) {

    LOG_SCOPE("mat_mult()", "PetscNonlinearSolver");
    
    PetscErrorCode ierr=0;
    
    libmesh_assert(mat);
    libmesh_assert(x);
    libmesh_assert(y);
    
    
    void * ctx = PETSC_NULL;
    
    // get the context
    ierr = MatShellGetContext(mat, ctx);
    
    __mast_multiphysics_petsc_shell_context
    *mat_ctx = static_cast<__mast_multiphysics_petsc_shell_context*> (ctx);
    
    // now complete the matrix vector multiplication operator
    
    
    return ierr;
}


//---------------------------------------------------------------
// this function is called by PETSc to evaluate the residual at X
PetscErrorCode
__mast_multiphysics_petsc_snes_residual (SNES snes, Vec x, Vec r, void * ctx)
{
    LOG_SCOPE("residual()", "PetscMultiphysicsNonlinearSolver");
    
    PetscErrorCode ierr=0;
    
    libmesh_assert(x);
    libmesh_assert(r);
    libmesh_assert(ctx);
    
    // No way to safety-check this cast, since we got a void *...
    MAST::MultiphysicsNonlinearSolverBase * solver =
    static_cast<MAST::MultiphysicsNonlinearSolverBase*> (ctx);
    
    // global communicator context
    MPI_Comm g_comm = solver->comm();
    libMesh::Parallel::Communicator comm(g_comm);
    
    // iterate over each system and evaluate the residual
    for (unsigned int i=0; i< solver->n_disciplines(); i++) {
        
        // system for this discipline
        libMesh::NonlinearImplicitSystem&
        sys = dynamic_cast<libMesh::NonlinearImplicitSystem&>
        (solver->get_system_assembly(i).system());
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        Vec sol, res;
        
        // extract the subvector for this system
        ierr = VecGetSubVector(x, sys_is, &sol);      CHKERRABORT(g_comm, ierr);
        ierr = VecGetSubVector(r, sys_is, &res);      CHKERRABORT(g_comm, ierr);
        
        std::auto_ptr<libMesh::NumericVector<Real> >
        sol_vec(new libMesh::PetscVector<Real>(sol, sys.comm())),
        res_vec(new libMesh::PetscVector<Real>(res, sys.comm()));
        
        // Use the system's update() to get a good local version of the
        // parallel solution.  This operation does not modify the incoming
        // "x" vector, it only localizes information from "x" into
        // sys.current_local_solution.
        //X_global.swap(X_sys);
        //sys.update();
        //X_global.swap(X_sys);
        
        // Enforce constraints (if any) exactly on the
        // current_local_solution.  This is the solution vector that is
        // actually used in the computation of the residual below, and is
        // not locked by debug-enabled PETSc the way that "x" is.
        //sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());
        
        solver->get_system_assembly(i).residual_and_jacobian (*sol_vec,
                                                              res_vec.get(),
                                                              nullptr,
                                                              sys);
        
        res_vec->close();
        
        // now restore the subvectors
        ierr = VecRestoreSubVector(x, sys_is, &sol);  CHKERRABORT(g_comm, ierr);
        ierr = VecRestoreSubVector(r, sys_is, &res);  CHKERRABORT(g_comm, ierr);
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
    
    // No way to safety-check this cast, since we got a void *...
    MAST::MultiphysicsNonlinearSolverBase * solver =
    static_cast<MAST::MultiphysicsNonlinearSolverBase*> (ctx);
    
    // global communicator context
    MPI_Comm g_comm = solver->comm();
    libMesh::Parallel::Communicator comm(g_comm);
    
    // iterate over each system and evaluate the residual
    for (unsigned int i=0; i< solver->n_disciplines(); i++) {
        
        // system for this discipline
        libMesh::NonlinearImplicitSystem&
        sys = dynamic_cast<libMesh::NonlinearImplicitSystem&>
        (solver->get_system_assembly(i).system());
        
        // get the IS for this system
        IS sys_is = solver->index_sets()[i];
        Vec sol;
        
        
        // extract the subvector for this system
        ierr = VecGetSubVector(x, sys_is, &sol);      CHKERRABORT(g_comm, ierr);
        
        std::auto_ptr<libMesh::NumericVector<Real> >
        sol_vec(new libMesh::PetscVector<Real>(sol, sys.comm()));
        
        // Use the system's update() to get a good local version of the
        // parallel solution.  This operation does not modify the incoming
        // "x" vector, it only localizes information from "x" into
        // sys.current_local_solution.
        //X_global.swap(X_sys);
        //sys.update();
        //X_global.swap(X_sys);
        
        // Enforce constraints (if any) exactly on the
        // current_local_solution.  This is the solution vector that is
        // actually used in the computation of the residual below, and is
        // not locked by debug-enabled PETSc the way that "x" is.
        //sys.get_dof_map().enforce_constraints_exactly(sys, sys.current_local_solution.get());

        //PetscMatrix<Number> PC(pc, sys.comm());
        //PC.attach_dof_map(sys.get_dof_map());
        //PC.close();
        
        solver->get_system_assembly(i).residual_and_jacobian (*sol_vec,
                                                              nullptr,
                                                              sys.matrix,
                                                              sys);
        
        sys.matrix->close();
        
        // now restore the subvectors
        ierr = VecRestoreSubVector(x, sys_is, &sol);  CHKERRABORT(g_comm, ierr);
    }
    
    return ierr;
}





MAST::MultiphysicsNonlinearSolverBase::
MultiphysicsNonlinearSolverBase(const std::string& nm,
                                unsigned int n):
_name(nm),
_n_disciplines(n),
_discipline_assembly(n, nullptr),
_is(_n_disciplines, PETSC_NULL),
_sub_mats(_n_disciplines*_n_disciplines, PETSC_NULL),
_n_dofs(0) {
    
}





MAST::MultiphysicsNonlinearSolverBase::~MultiphysicsNonlinearSolverBase() {
    
}




void
MAST::MultiphysicsNonlinearSolverBase::
set_system_assembly(unsigned int i,
                    MAST::NonlinearImplicitAssembly& assembly) {
    
    // make sure that the index is within bounds
    libmesh_assert_less(i, _n_disciplines);
    
    // also make sure that the specific discipline has not been set already
    libmesh_assert(!_discipline_assembly[i]);
    
    _discipline_assembly[i] = &assembly;
}



MAST::NonlinearImplicitAssembly&
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
    bool p = false;
    for (unsigned int i=0; i<_n_disciplines; i++)
        p = _discipline_assembly[i] != nullptr;
    libmesh_assert(p);
    
    //////////////////////////////////////////////////////////////////////
    // create a global communicator spanning all disciplines, which will be
    // used for the snes context
    //////////////////////////////////////////////////////////////////////
    
    MPI_Comm_group(_discipline_assembly[0]->system().comm().get(), &_g_union);
    for (unsigned int i=1; i<_n_disciplines; i++) {
        
        MPI_Group  g1, g2;
        
        MPI_Comm_group(_discipline_assembly[i]->system().comm().get(), &g1);
        MPI_Group_union(g1, _g_union, &g2);
        MPI_Group_free(&g1);
        MPI_Group_free(&_g_union);
        _g_union = g2;
    }
    
    // g_union should now have the union of all groups over all disciplines.
    // Use it to create a new communicator for the snes context
    MPI_Comm_create_group(MPI_COMM_WORLD, _g_union, 0, &_g_comm);

    
    
    
    //////////////////////////////////////////////////////////////////////
    // create a petsc nested matrix of size NxN
    //////////////////////////////////////////////////////////////////////
    PetscErrorCode   ierr;
    const bool       sys_name = libMesh::on_command_line("--solver_system_names");
    std::string      nm;
    std::vector<__mast_multiphysics_petsc_shell_context>
    mat_ctx(_n_disciplines*_n_disciplines);
    
    
    // all diagonal blocks use the system matrcices, while shell matrices
    // are created for the off-diagonal terms.
    for (unsigned int i=0; i<_n_disciplines; i++)
        for (unsigned int j=0; j<_n_disciplines; j++) {
    
        libMesh::NonlinearImplicitSystem& sys =
        dynamic_cast<libMesh::NonlinearImplicitSystem&>(_discipline_assembly[i]->system());
        
        // add the number of dofs in this system to the global count
        _n_dofs   +=  sys.n_dofs();
        
            if (i==j) {
                
                // the diagonal matrix
                _sub_mats[i*_n_disciplines+i] =
                dynamic_cast<libMesh::PetscMatrix<Real>*>(sys.matrix)->mat();
            }
            else {
                
                libMesh::System&
                sys_j = _discipline_assembly[j]->system();
                
                PetscInt
                mat_i_m    = sys.get_dof_map().n_dofs(),
                mat_i_n_l  = sys.get_dof_map().n_dofs_on_processor(sys.processor_id()),
                mat_j_m    = sys_j.get_dof_map().n_dofs(),
                mat_j_n_l  = sys_j.get_dof_map().n_dofs_on_processor(sys.processor_id());
                
                // the off-diagonal matrix
                ierr = MatCreateShell(_g_comm,
                                      mat_i_n_l,
                                      mat_j_n_l,
                                      mat_i_m,
                                      mat_j_n_l,
                                      PETSC_NULL,
                                      &_sub_mats[i*_n_disciplines+j]);
                CHKERRABORT(_g_comm, ierr);
                
                // initialize the context and tell the matrix about it
                mat_ctx[i*_n_disciplines+j].i      = i;
                mat_ctx[i*_n_disciplines+j].j      = j;
                mat_ctx[i*_n_disciplines+j].solver = this;
                
                ierr = MatShellSetContext(_sub_mats[i*_n_disciplines+j],
                                          &mat_ctx[i*_n_disciplines+j]);
                CHKERRABORT(_g_comm, ierr);
                
                // set the mat-vec multiplication operation for this
                // shell matrix
                ierr = MatShellSetOperation(_sub_mats[i*_n_disciplines+j],
                                            MATOP_MULT,
                                            (void(*)(void))__mast_multiphysics_petsc_mat_mult);
                CHKERRABORT(_g_comm, ierr);
            }
    }

    
    ierr = MatCreateNest(_g_comm,
                         _n_disciplines, PETSC_NULL,
                         _n_disciplines, PETSC_NULL,
                         &_sub_mats[0],
                         &_mat);
    CHKERRABORT(_g_comm, ierr);
    
    if (sys_name) {
        
        nm = this->name() + "_";
        MatSetOptionsPrefix(_mat, nm.c_str());
    }
    ierr = MatSetFromOptions(_mat);                    CHKERRABORT(_g_comm, ierr);
    
    
    // get the IS belonging to each block
    ierr  =    MatNestGetISs(_mat, &_is[0], PETSC_NULL);CHKERRABORT(_g_comm, ierr);
    
    
    
    
    //////////////////////////////////////////////////////////////////////
    // setup the vector for solution
    //////////////////////////////////////////////////////////////////////
    ierr = VecCreate(_g_comm, &_sol);                   CHKERRABORT(_g_comm, ierr);
    ierr = VecSetSizes(_sol, PETSC_DECIDE, _n_dofs);    CHKERRABORT(_g_comm, ierr);
    ierr = VecSetType(_sol, VECMPI);                   CHKERRABORT(_g_comm, ierr);
    ierr = VecDuplicate(_sol, &_res);                   CHKERRABORT(_g_comm, ierr);
    
    
    //////////////////////////////////////////////////////////////////////
    // now initialize the vector from the system solutions, which will serve
    // as the initial solution for the coupled system
    //////////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        
    }
    
    
    ierr = MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);  CHKERRABORT(_g_comm, ierr);
    ierr = MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);    CHKERRABORT(_g_comm, ierr);
    ierr = VecAssemblyBegin(_res);                      CHKERRABORT(_g_comm, ierr);
    ierr = VecAssemblyEnd(_res);                        CHKERRABORT(_g_comm, ierr);
    
    
    
    //////////////////////////////////////////////////////////////////////
    // create the solver context
    //////////////////////////////////////////////////////////////////////
    SNES           snes;
    KSP            ksp;
    PC             pc;

    
    // setup the solver context
    ierr = SNESCreate(_g_comm, &snes);                 CHKERRABORT(_g_comm, ierr);
    
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
    ierr = SNESSetFromOptions(snes);                  CHKERRABORT(_g_comm, ierr);

    
    
    // setup the ksp
    ierr = SNESGetKSP (snes, &ksp);                   CHKERRABORT(_g_comm, ierr);
    ierr = KSPSetFromOptions(ksp);                    CHKERRABORT(_g_comm, ierr);


    
    // setup the pc
    ierr = KSPGetPC(ksp, &pc);                        CHKERRABORT(_g_comm, ierr);
    
    for (unsigned int i=0; i<_n_disciplines; i++) {
        
        if (sys_name) {
            
            nm = _discipline_assembly[i]->system().name();
            ierr = PCFieldSplitSetIS(pc, nm.c_str(), _is[i]); CHKERRABORT(_g_comm, ierr);
        }
        else
            ierr = PCFieldSplitSetIS(pc, NULL, _is[i]);CHKERRABORT(_g_comm, ierr);
    }
    ierr = PCSetFromOptions(pc);                      CHKERRABORT(_g_comm, ierr);
    
    //////////////////////////////////////////////////////////////////////
    // now, solve
    //////////////////////////////////////////////////////////////////////
    START_LOG("SNESSolve", this->name()+"_MultiphysicsSolve");
    
    // now solve
    ierr = SNESSolve(snes, PETSC_NULL, _sol);
    
    STOP_LOG("SNESSolve", this->name()+"_MultiphysicsSolve");
    
    
    // destroy the Petsc contexts
    ierr = SNESDestroy(&snes);                        CHKERRABORT(_g_comm, ierr);
    for (unsigned int i=0; i<_n_disciplines; i++)
        for (unsigned int j=0; j<_n_disciplines; j++)
            if (i != j) {
                ierr = MatDestroy(&_sub_mats[i*_n_disciplines+j]);
                CHKERRABORT(_g_comm, ierr);
            }
    
    ierr = MatDestroy(&_mat);                          CHKERRABORT(_g_comm, ierr);
    ierr = VecDestroy(&_sol);                          CHKERRABORT(_g_comm, ierr);
    ierr = VecDestroy(&_res);                          CHKERRABORT(_g_comm, ierr);
    
    
    // destroy the MPI contexts
    MPI_Group_free(&_g_union);
    MPI_Comm_free(&_g_comm);
}

