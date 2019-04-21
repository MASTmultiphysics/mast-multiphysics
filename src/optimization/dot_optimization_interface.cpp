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
#include "optimization/dot_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "base/mast_config.h"


MAST::DOTOptimizationInterface::DOTOptimizationInterface():
MAST::OptimizationInterface() {
    
#if MAST_ENABLE_DOT == 0
    libmesh_error_msg("MAST configured without DOT support.");
#endif
}



void
MAST::DOTOptimizationInterface::optimize() {
    
#if MAST_ENABLE_DOT == 1
    int
    INFO    = 0,
    METHOD  = 0,      //  METHOD == 0 or 1 means MMFD for constrained
                      //  METHOD == 2      means SLP  for constrained
                      //  METHOD == 3      means SQP
                      //  METHOD == 0 or 1 means BFGS for unconstrained
                      //  METHOD == 2      means Fletcher-Reeves for unconstrained
    IPRINT  = 4,      //  IPRINT = 0, no output
                      //  IPRINT = 1, internal params, initial information and results
                      //  IPRINT = 2, same plus obj, X vec at each iter
                      //  IPRINT = 3, same plus g, and critical constraint numbers
                      //  IPRINT = 4, same plus gradients
                      //  IPRINT = 5, same plus search direction
                      //  IPRINT = 6, same plus set IPRM(11) = 1 and IPRM(12) = 1
                      //  IPRINT = 7, same except set IPRM(12) = 2
    NDV     = _feval->n_vars(),
    NCON    = _feval->n_eq() + _feval->n_ineq(),
    MINMAX  = 0,      //  MINMAX = 0,-1 for minimization, = 1 for maximization
    NRWK    = std::max(NDV*NCON*10, 100000), // add a factor of 10 to be safe.
    NRIWK   = std::max(NDV*NCON*10, 300); // May need to be changed for individual problems

    libmesh_assert_greater(NDV,  0);
    libmesh_assert_greater(NCON, 0);

    Real
    OBJ     = 0.;
    
    std::vector<int>
    IPRM   (20,   0),
    IWK    (NRIWK, 0);

    std::vector<double>
    X     (NDV,      0.),
    XL    (NDV,      0.),
    XU    (NDV,      0.),
    G     (NCON,     0.),
    RPRM  (20,       0.),
    WK    (NRWK,     0.),
    DF0DX (NDV,      0.),
    DFDX  (NDV*NCON, 0.);
    
    std::vector<bool>
    eval_grads(NCON, false);

    
    _feval->init_dvar(X, XL, XU);

    unsigned int
    ITER = 0;
    
    bool
    obj_grad = false,
    if_cont  = true;
    
    
    // user provided gradients
    IPRM[0]  = 1;
    
    while (if_cont) {

        dot_(&INFO,
             &METHOD,
             &IPRINT,
             &NDV,
             &NCON,
             &X[0],
             &XL[0],
             &XU[0],
             &OBJ,
             &MINMAX,
             &G[0],
             &RPRM[0],
             &IPRM[0],
             &WK[0],
             &NRWK,
             &IWK[0],
             &NRIWK);

        
        if (INFO == 0) // INFO == 0 means optimization is complete
                       // INFO == 0 and IPRM(18) > 0 means an error has occurred
            if_cont = false;
        else if (INFO == 1) { // evaluate objective and constraints
            std::fill(eval_grads.begin(), eval_grads.end(), false);
            obj_grad = false;
        }
        else if (INFO == 2) {
            std::fill(eval_grads.begin(), eval_grads.end(), true);
            obj_grad = true;
        }
        
        _feval->evaluate(X,
                         OBJ, obj_grad, DF0DX,
                         G, eval_grads, DFDX);
        
        // if gradients were requested, copy the data back to WK
        if (INFO == 2) {
            // first the objective function gradients
            for (unsigned int i=0; i<NDV; i++)   WK[i] = DF0DX[i];
            
            // next, the constraint gradients
            // DOT requires gradients for only the active constraints
            unsigned int
            n_active_constrs = IPRM[19],
            active_constr_id = 0;
            
            // the active constraints are listed in the first elements of
            // IWK
            for (unsigned int i=0; i<n_active_constrs; i++) {
                active_constr_id  = IWK[i];
                // now copy the gradient values of this constraint
                // wrt all DVs, which will be the row of DFDX
                for (unsigned int j=0; j<NDV; j++)
                    WK[(i+1)*NDV+j] = DFDX[j*NCON+active_constr_id-1];
            }
        }
        
        _feval->_output_wrapper(ITER, X, OBJ, G, true);
        ITER = ITER + 1;
    }
    
    // write the final iteration to the output
    _feval->_output_wrapper(ITER, X, OBJ, G, true);
    
#endif  // MAST_ENABLE_DOT 1
}
