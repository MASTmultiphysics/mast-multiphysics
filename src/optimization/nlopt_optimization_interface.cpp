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
#include "optimization/nlopt_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "base/mast_config.h"

double
_mast_nlopt_objective_function(unsigned n,
                               const double* x,
                               double* grad,
                               void* f_data) {
    
    MAST::NLOptOptimizationInterface*
    opt_interface = static_cast<MAST::NLOptOptimizationInterface*>(f_data);

    return
    opt_interface->objective_evaluation(n, x, grad);
}


void
_mast_nlopt_ineq_constr_mfunc (unsigned m,
                               double *result,
                               unsigned n,
                               const double *x,
                               double *gradient, /* NULL if not needed */
                               void *func_data) {
    //
    //  \partial c_i/\partial x_j is stored in grad[i*n + j]
    //
    MAST::NLOptOptimizationInterface*
    opt_interface = static_cast<MAST::NLOptOptimizationInterface*>(func_data);
    
    opt_interface->inequality_constraint_evaluation(m, result, n, x, gradient);
}



MAST::NLOptOptimizationInterface::NLOptOptimizationInterface(nlopt_algorithm alg):
MAST::OptimizationInterface(),
_iter (0),
_alg  (alg) {
    
}


void
MAST::NLOptOptimizationInterface::optimize() {
    
    // make sure that all processes have the same problem setup
    _feval->sanitize_parallel();
    
    int
    n                  = _feval->n_vars(),
    m_eq               = _feval->n_eq(),
    m_ineq             = _feval->n_ineq(),
    n_rel_change_iters = _feval->n_iters_relative_change();
    _iter              = 0;
    
    // equality constraints are not currently handled in this API
    libmesh_assert_equal_to(m_eq, 0);
    
    Real
    obj  = 0.,
    tol  = 1.e-3;
    
    nlopt_result
    res;
    
    libmesh_assert_greater(n, 0);
    
    std::vector<Real>
    xval  (n, 0.),
    xmin  (n, 0.),
    xmax  (n, 0.);
    
    
    //////////////////////////////////////////////////////////
    // initialize the optimizer object
    //////////////////////////////////////////////////////////
    nlopt_opt
    opt  =  nlopt_create(_alg, n);
    
    //////////////////////////////////////////////////////////
    // set variable bounds
    //////////////////////////////////////////////////////////
    _feval->_init_dvar_wrapper(xval, xmin, xmax);

    res = nlopt_set_lower_bounds(opt, &xmin[0]);
    libmesh_assert_equal_to(res, NLOPT_SUCCESS);
    
    res = nlopt_set_upper_bounds(opt, &xmax[0]);
    libmesh_assert_equal_to(res, NLOPT_SUCCESS);
    
    //////////////////////////////////////////////////////////
    // attach the objective and constraint functions
    //////////////////////////////////////////////////////////
    res = nlopt_set_min_objective(opt, _mast_nlopt_objective_function, (void*)this);
    libmesh_assert_equal_to(res, NLOPT_SUCCESS);
    
    // nothing done about equality constraints so far.
    res = nlopt_add_inequality_mconstraint(opt,
                                           m_ineq,
                                           _mast_nlopt_ineq_constr_mfunc,
                                           (void*)this,
                                           &tol);
    libmesh_assert_equal_to(res, NLOPT_SUCCESS);
    
    
    //////////////////////////////////////////////////////////
    // optimize
    //////////////////////////////////////////////////////////
    res = nlopt_optimize(opt, &xval[0], &obj);
    libmesh_assert_equal_to(res, NLOPT_SUCCESS);
    
    //////////////////////////////////////////////////////////
    // clean up the optimizer object
    //////////////////////////////////////////////////////////
    nlopt_destroy(opt);
}



Real
MAST::NLOptOptimizationInterface::objective_evaluation(unsigned n,
                                                       const double* x,
                                                       double* grad) {
    
    libmesh_assert_equal_to(n, _feval->n_vars());
    
    unsigned int
    n_constr = _feval->n_eq() + _feval->n_ineq();
    
    std::vector<Real>
    xvals(x, x+n),
    fvals(n_constr, 0.),
    df0dx(n, 0.),
    dfdx (n*n_constr, 0.);
  
    std::vector<bool>
    eval_grads(n_constr, false);
    
    Real
    f0val = 0.;

    // output the data
    if (!grad) {
        _feval->_output_wrapper(_iter, xvals, f0val, fvals, true);
        _iter++;
    }
    
    _feval->_evaluate_wrapper(xvals,
                              f0val, grad!=nullptr, df0dx,
                              fvals, eval_grads, dfdx);

    if (grad)
        for (unsigned int i=0; i<n; i++) grad[i] = df0dx[i];
    
    return f0val;
}


void
MAST::NLOptOptimizationInterface::
inequality_constraint_evaluation(unsigned m,
                                 double *result,
                                 unsigned n,
                                 const double *x,
                                 double *gradient) {
    
    
    libmesh_assert_equal_to(n, _feval->n_vars());
    
    unsigned int
    n_constr = _feval->n_eq() + _feval->n_ineq();
    
    libmesh_assert_equal_to(m, n_constr);
    
    std::vector<Real>
    xvals(x, x+n),
    fvals(n_constr, 0.),
    df0dx(n, 0.),
    dfdx (n*n_constr, 0.);
    
    std::vector<bool>
    eval_grads(n_constr, gradient!=nullptr);
    
    Real
    f0val;
    
    _feval->_evaluate_wrapper(xvals,
                              f0val, false, df0dx,
                              fvals, eval_grads, dfdx);
    
    if (gradient) {

        //
        //  NLOpt requires the derivatives to be in this form
        //  \partial c_i/\partial x_j is stored in grad[i*n + j]
        //
        //   However, the function evaluation returns the following:
        //   grads(k): Derivative of f_i(x) with respect
        //   to x_j, where k = (j-1)*M + i.
        //
        //   Therefore, we translate the values into the proper order

        for (unsigned int i=0; i<m; i++)
            for (unsigned int j=0; j<n; j++)
                gradient[i*n+j] = dfdx[j*m+i];
    }
}

