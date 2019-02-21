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
#include "solver/continuation_solver_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/assembly_elem_operation.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/linear_solver.h"


MAST::ContinuationSolverBase::ContinuationSolverBase():
_initialized (false),
max_it       (20),
abs_tol      (1.e-8),
rel_tol      (1.e-8),
arc_length   (0.),
_elem_ops    (nullptr),
_assembly    (nullptr),
_p           (nullptr),
_p0          (0.) {
    
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
        if (iter >= max_it)       cont = false;
    }
        
    libMesh::out
    << std::setw(10) << "iter: "
    << std::setw(5)  << iter
    << std::setw(10) << "res-l2: "
    << std::setw(15) << norm
    << std::setw(20) << "relative res-l2: "
    << std::setw(15) << norm/norm0
    << std::setw(20) << "Terminated"  << std::endl;
    
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
    r1(X.clone().release());
    
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
        _assembly->sensitivity_assemble(p, dfdp);

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
                                                          dXdp,
                                                          /* homogeneous = */ true);
#endif
    }

    //////////////////////////////////////////////////////////
    //         STEP 3:  a   = dg/dp + dg/dX dXdp
    //////////////////////////////////////////////////////////
    a   = dgdp + dgdX.dot(dXdp);
    libmesh_assert_greater(a, 0.);
    
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
                                                      dX.get(),
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

