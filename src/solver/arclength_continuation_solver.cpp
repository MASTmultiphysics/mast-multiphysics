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
#include "solver/arclength_continuation_solver.h"
#include "base/assembly_base.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/linear_solver.h"
#include "libmesh/dof_map.h"

MAST::ArclengthContinuationSolver::ArclengthContinuationSolver():
MAST::ContinuationSolverBase(),
_dpds_sign  (1.) {
    
}


MAST::ArclengthContinuationSolver::~ArclengthContinuationSolver() {
    
}


void
MAST::ArclengthContinuationSolver::initialize(Real dp) {
    
    libmesh_assert(!_initialized);
    libmesh_assert(_elem_ops);
    
    (*_p)() += dp;
    if (dp < 0.) _dpds_sign = -1.;
    
    MAST::NonlinearSystem&
    system = _assembly->system();
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    dX(system.solution->clone().release());
    
    _assembly->system().solve(*_elem_ops, *_assembly);
    
    dX->add(-1., *system.solution);
    dX->scale(-1.);
    dX->close();

    // initialize scaling factors
    _X_scale   = 1./dX->l2_norm();
    _p_scale   = 1./std::fabs(dp);
    
    arc_length = std::sqrt(std::pow(_p_scale, 2) * dp*dp +
                           std::pow(_X_scale, 2) * std::pow(dX->l2_norm(), 2));
    
    libmesh_assert_greater(arc_length, 0.);
    
    _initialized = true;
}


void
MAST::ArclengthContinuationSolver::
_solve_NR_iterate(libMesh::NumericVector<Real>       &X,
                  MAST::Parameter                    &p) {
    
    libmesh_assert(_initialized);

    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release()),
    dgdX(X.zero_clone().release()),
    dfdp(X.zero_clone().release()),
    dXdp(X.zero_clone().release()),
    dX(X.zero_clone().release());
    
    libMesh::SparseMatrix<Real>
    &jac = *_assembly->system().matrix;
    
    Real
    g    = 0.,
    dgdp = 0.,
    dp   = 0.;
    
    _g(X, p, *dfdp, *dXdp, g, dgdp, dgdX.get());
    
    if (!schur_factorization)
        _solve(X, p,
               *f,                           true,  // update f
               *dfdp,                        false, // update dfdp
               *dgdX, dgdp, g,
               *dX, dp);
    else
        _solve_schur_factorization(X, p,
                                   jac,                          true,  // update jac
                                   *f,                           true,  // update f
                                   *dfdp,                        false, // update dfdp
                                   *dXdp,                        false, // update dXdp
                                   *dgdX, dgdp, g,
                                   *dX, dp);

    // update the solution and load parameter
    p() += dp;
    X.add(1., *dX);
    X.close();
}


void
MAST::ArclengthContinuationSolver::
_dXdp(const libMesh::NumericVector<Real> &X,
      const MAST::Parameter              &p,
      libMesh::NumericVector<Real>       &dfdp,
      libMesh::NumericVector<Real>       &dXdp) {

    libmesh_assert(_elem_ops);
    
    MAST::NonlinearSystem
    &system    = _assembly->system();

    libmesh_assert(system.operation() == MAST::NonlinearSystem::NONE);

    _assembly->set_elem_operation_object(*_elem_ops);

    system.set_operation(MAST::NonlinearSystem::FORWARD_SENSITIVITY_SOLVE);

    _assembly->sensitivity_assemble(*system.solution, true, p, dfdp);

    libMesh::SparseMatrix<Real>
    *pc  = system.request_matrix("Preconditioner");

    std::pair<unsigned int, Real>
    solver_params = system.get_linear_solve_parameters(),
    rval;

    dXdp.zero();
    
    rval = system.linear_solver->solve (*system.matrix, pc,
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
    
    _assembly->clear_elem_operation_object();
    system.set_operation(MAST::NonlinearSystem::NONE);
}




Real
MAST::ArclengthContinuationSolver::_g(const libMesh::NumericVector<Real> &X,
                                      const MAST::Parameter              &p) {
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    dfdp(X.zero_clone().release()),
    dXdp(X.zero_clone().release());
    
    Real
    g    = 0.,
    dgdp = 0.;
    
    _g(X, p, *dfdp, *dXdp, g, dgdp, nullptr);
    
    return g;
}


void
MAST::ArclengthContinuationSolver::_g(const libMesh::NumericVector<Real> &X,
                                      const MAST::Parameter              &p,
                                      libMesh::NumericVector<Real>       &dfdp,
                                      libMesh::NumericVector<Real>       &dXdp,
                                      Real                               &g,
                                      Real                               &dgdp,
                                      libMesh::NumericVector<Real>       *dgdX) {

    libmesh_assert(_initialized);

    // update the constraint data
    _dXdp(X, p, dfdp, dXdp);
    
    // this includes scaling of X and p
    Real
    dpds = _dpds_sign * std::sqrt(1./ ( std::pow(_X_scale/_p_scale,2) * dXdp.dot(dXdp) + 1.));
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    dX(X.clone().release());
    
    dX->add(-1., *_X0);
    dX->close();
    
    // (dX/ds)_scaled = (dX/dp)_scaled * (dp/ds)_scaled
    //                = (dX_scaled/dX) * dX/dp * (dp/dp_scaled) * (dp/ds)_scaled
    //                = X_scale/p_scale * dX/dp * (dp/ds)_scaled
    g    = (_X_scale * (dX->dot(dXdp) * dpds * _X_scale/_p_scale) +
            _p_scale * (p() - _p0)    * dpds) - arc_length;
    dgdp = _p_scale * dpds;
    
    if (dgdX) {
        
        dgdX->zero();
        dgdX->add(_X_scale * (dpds * _X_scale/_p_scale), dXdp);
        dgdX->close();
    }
}

