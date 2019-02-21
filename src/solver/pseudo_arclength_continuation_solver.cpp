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
#include "solver/pseudo_arclength_continuation_solver.h"
#include "base/assembly_base.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/linear_solver.h"


MAST::PseudoArclengthContinuationSolver::PseudoArclengthContinuationSolver():
MAST::ContinuationSolverBase(),
_t0_p   (0.) {
    
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
    
    _assembly->system().solve(*_elem_ops, *_assembly);
    
    _t0_X->add(-1., *system.solution);
    _t0_X->scale(-1.);
    _t0_X->close();
    
    arc_length = std::sqrt(dp*dp + std::pow(_t0_X->l2_norm(), 2));

    libmesh_assert_greater(arc_length, 0.);
    
    // scale _t0_X such that {t0_X, t0_p} is a unit vector
    _t0_X->scale(1./arc_length);
    _t0_X->close();
    _t0_p = dp/arc_length;
    
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
                             *dfdp,
                             *dXdp,
                             *t1_X,
                             t1_p);
    
    
    _g(X, p, *t1_X, t1_p, g);
    
    _solve_schur_factorization(X, p,
                               jac,                          false, // do not update jac
                               *f,                           true,  // update f
                               *dfdp,                        false, // do not update dfdp
                               *dXdp,                        false, // do not update dXdp
                               *t1_X, t1_p, g,           // dgdX = t1^X, dgdp = t1^p
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
                         libMesh::NumericVector<Real>       &dfdp,
                         libMesh::NumericVector<Real>       &dXdp,
                         libMesh::NumericVector<Real>       &t1_X,
                         Real                               &t1_p) {
    
    libmesh_assert(_initialized);

    std::unique_ptr<libMesh::NumericVector<Real>>
    f(X.zero_clone().release());

    // first update the search direction
    _solve_schur_factorization(X, p,
                               jac,               true,  // update jac
                               *f,                false, // do not update f
                               dfdp,              true,  // update dfdp
                               dXdp,              true,  // update dXdp
                               *_t0_X, _t0_p, 1.,  // dgdX = t0^X, dgdp = t0^p, g = 1
                               t1_X, t1_p);

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
    dfdp(X.zero_clone().release()),
    dXdp(X.zero_clone().release()),
    t1_X(X.zero_clone().release());
    
    Real
    g       = 0.,
    t1_p    = 0.;

    _update_search_direction(X, p,
                             *_assembly->system().matrix,
                             *dfdp,
                             *dXdp,
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
    
    g    = dX->dot(t1_X) + (p() - _p0) * t1_p - arc_length;
}

