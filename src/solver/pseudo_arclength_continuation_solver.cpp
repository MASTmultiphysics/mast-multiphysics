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
#include "solver/pseudo_arclength_continuation_solver.h"
#include "base/assembly_base.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/linear_solver.h"


MAST::PseudoArclengthContinuationSolver::PseudoArclengthContinuationSolver():
MAST::ContinuationSolverBase(),
_t0_p       (0.),
_t0_p_orig  (0.) {
    
}


MAST::PseudoArclengthContinuationSolver::~PseudoArclengthContinuationSolver() {
    
}


void
MAST::PseudoArclengthContinuationSolver::initialize(Real dp) {

    libmesh_assert(!_initialized);
    libmesh_assert(_elem_ops);
    
    (*_p)() += dp;

    MAST::NonlinearSystem&
    system = _assembly->system();
    
    _t0_X.reset(system.solution->clone().release());
    _t0_X_orig.reset(system.solution->zero_clone().release());
    
    _assembly->system().solve(*_elem_ops, *_assembly);
    
    _t0_X->add(-1., *system.solution);
    _t0_X->scale(-1.);
    _t0_X->close();
    
    // initialize scaling factors
    _X_scale   = 1./_t0_X->l2_norm();
    _p_scale   = 1./std::fabs(dp);

    arc_length = std::sqrt(std::pow(_p_scale, 2) * dp*dp +
                           std::pow(_X_scale, 2) * std::pow(_t0_X->l2_norm(), 2));

    libmesh_assert_greater(arc_length, 0.);
    
    // scale _t0_X such that {t0_X, t0_p} is a unit vector
    _t0_X->scale(_X_scale/arc_length);
    _t0_X->close();
    _t0_p = _p_scale*dp/arc_length;

    _save_iteration_data();
    
    _initialized = true;
}


void
MAST::PseudoArclengthContinuationSolver::
_solve_NR_iterate(libMesh::NumericVector<Real>       &X,
                  MAST::Parameter                    &p) {

    libmesh_assert(_initialized);

    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release()),
    dfdp(X.zero_clone().release()),
    dXdp(X.zero_clone().release()),
    t1_X(X.zero_clone().release()),
    dX(X.zero_clone().release());
    
    libMesh::SparseMatrix<Real>
    &jac = *_assembly->system().matrix;
    
    Real
    t1_p = 0.,
    g    = 0.,
    dp   = 0.;
    
    _update_search_direction(X, p,
                             jac,
                             *t1_X,
                             t1_p);
    
    
    _g(X, p, *t1_X, t1_p, g);
    
    // scale the values so that they are in the scaled coordinates
    t1_X->scale(_X_scale);
    t1_X->close();
    t1_p *= _p_scale;
    
    if (!schur_factorization)
        _solve(X, p,
               *f,                           true,  // update f
               *dfdp,                        true,  // update dfdp
               *t1_X, t1_p, g,           // dgdX = X_scale * t1^X, dgdp = p_scale *t1^p
               *dX, dp);
    else
        _solve_schur_factorization(X, p,
                                   jac,                          false, // do not update jac
                                   *f,                           true,  // update f
                                   *dfdp,                        true,  // update dfdp
                                   *dXdp,                        true,  // update dXdp
                                   *t1_X, t1_p, g,           // dgdX = X_scale * t1^X, dgdp = p_scale *t1^p
                                   *dX, dp);
    
    // update the solution and load parameter
    p() += dp;
    X.add(1., *dX);
    X.close();
    
}


void
MAST::PseudoArclengthContinuationSolver::
_update_search_direction(const libMesh::NumericVector<Real> &X,
                         const MAST::Parameter              &p,
                         libMesh::SparseMatrix<Real>        &jac,
                         libMesh::NumericVector<Real>       &t1_X,
                         Real                               &t1_p) {
    
    libmesh_assert(_initialized);

    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release()),
    dfdp(X.zero_clone().release()),
    dXdp(X.zero_clone().release());

    MAST::NonlinearSystem
    &system    = _assembly->system();

    system.set_operation(MAST::NonlinearSystem::FORWARD_SENSITIVITY_SOLVE);
    _assembly->set_elem_operation_object(*_elem_ops);

    _assembly->sensitivity_assemble(*system.solution, true, p, *dfdp);
    dfdp->scale(_X_scale/_p_scale);
    dfdp->close();

    _assembly->clear_elem_operation_object();
    system.set_operation(MAST::NonlinearSystem::NONE);


    // first update the search direction
    if (!schur_factorization)
        _solve(X, p,
               *f,                false, // do not update f
               *dfdp,             false, // do not update dfdp
               *_t0_X, _t0_p,  -1.,  // dgdX = t0^X, dgdp = t0^p, g = -1
               t1_X, t1_p);
    else
        _solve_schur_factorization(X, p,
                                   jac,               true,  // update jac
                                   *f,                false, // do not update f
                                   *dfdp,             false, // do not update dfdp
                                   *dXdp,             true,  // update dXdp
                                   *_t0_X, _t0_p,  -1.,  // dgdX = t0^X, dgdp = t0^p, g = -1
                                   t1_X, t1_p);

    // now scale the vector for unit magnitude
    Real
    val = std::sqrt(t1_p*t1_p + std::pow(t1_X.l2_norm(), 2));
    
    t1_X.scale(1./val);
    t1_X.close();
    t1_p /= val;

    // store values for next iterate
    _t0_X->zero();
    _t0_X->add(1., t1_X);
    _t0_X->close();
    
    _t0_p = t1_p;
}


Real
MAST::PseudoArclengthContinuationSolver::_g(const libMesh::NumericVector<Real> &X,
                                            const MAST::Parameter              &p) {
    
    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release()),
    t1_X(X.zero_clone().release());
    
    Real
    g       = 0.,
    t1_p    = 0.;

    _update_search_direction(X, p,
                             *_assembly->system().matrix,
                             *t1_X,
                             t1_p);
    
    _g(X, p, *t1_X, t1_p, g);
    
    return g;
}


void
MAST::PseudoArclengthContinuationSolver::_g(const libMesh::NumericVector<Real> &X,
                                            const MAST::Parameter              &p,
                                            libMesh::NumericVector<Real>       &t1_X,
                                            Real                               &t1_p,
                                            Real                               &g) {

    libmesh_assert(_initialized);

    std::unique_ptr<libMesh::NumericVector<Real>>
    dX(X.clone().release());
    
    dX->add(-1., *_X0);
    dX->close();
    
    g    = (_X_scale *  dX->dot(t1_X) +
            _p_scale * (p() - _p0) * t1_p) - arc_length;
}


void
MAST::PseudoArclengthContinuationSolver::_save_iteration_data() {
    
    // copy for possible reuse in resetting step
    _t0_X_orig->zero();
    _t0_X_orig->add(1., *_t0_X);
    _t0_X_orig->close();
    _t0_p_orig = _t0_p;
}


void
MAST::PseudoArclengthContinuationSolver::_reset_iterations() {
    
    libmesh_assert(_initialized);
    
    _t0_X->zero();
    _t0_X->add(1., *_t0_X_orig);
    _t0_X->close();
    _t0_p = _t0_p_orig;
}
