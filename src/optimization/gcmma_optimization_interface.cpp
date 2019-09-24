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
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "base/mast_config.h"


MAST::GCMMAOptimizationInterface::GCMMAOptimizationInterface():
MAST::OptimizationInterface(),
_constr_penalty       (5.e1),
_initial_rel_step     (5.e-1),
_asymptote_reduction  (0.7),
_asymptote_expansion  (1.2),
_max_inner_iters      (15)   {
    
#if MAST_ENABLE_GCMMA == 0
    libmesh_error_msg("MAST configured without GCMMA support.");
#endif
}


void
MAST::GCMMAOptimizationInterface::set_real_parameter(const std::string &nm, Real val) {

    if (nm == "constraint_penalty") {
        
        libmesh_assert_greater(val, 0.);
        
        _constr_penalty = val;
    }
    else if (nm == "initial_rel_step") {
        
        libmesh_assert_greater(val, 0.);
        
        _initial_rel_step = val;
    }
    else if (nm == "asymptote_reduction") {
        
        libmesh_assert_greater(val, 0.);
        
        _asymptote_reduction = val;
    }
    else if (nm == "asymptote_expansion") {
        
        libmesh_assert_greater(val, 0.);
        
        _asymptote_expansion = val;
    }
    else
        libMesh::out
        << "Unrecognized real parameter: " << nm << std::endl;
}


void
MAST::GCMMAOptimizationInterface::set_integer_parameter(const std::string &nm, int val) {
    
    if (nm == "max_inner_iters") {
        
        libmesh_assert_greater(val, 0);
        
        _max_inner_iters = val;
    }
    else
        libMesh::out
        << "Unrecognized integer parameter: " << nm << std::endl;
}



void
MAST::GCMMAOptimizationInterface::optimize() {
#if MAST_ENABLE_GCMMA == 1

    // make sure that all processes have the same problem setup
    _feval->sanitize_parallel();
    
    int
    N                  = _feval->n_vars(),
    M                  = _feval->n_eq() + _feval->n_ineq(),
    n_rel_change_iters = _feval->n_iters_relative_change();
    
    libmesh_assert_greater(N, 0);
    
    std::vector<Real>  XVAL(N, 0.), XOLD1(N, 0.), XOLD2(N, 0.),
    XMMA(N, 0.), XMIN(N, 0.), XMAX(N, 0.), XLOW(N, 0.), XUPP(N, 0.),
    ALFA(N, 0.), BETA(N, 0.), DF0DX(N, 0.),
    A(M, 0.), B(M, 0.), C(M, 0.), Y(M, 0.), RAA(M, 0.), ULAM(M, 0.),
    FVAL(M, 0.), FAPP(M, 0.), FNEW(M, 0.), FMAX(M, 0.),
    DFDX(M*N, 0.), P(M*N, 0.), Q(M*N, 0.), P0(N, 0.), Q0(N, 0.),
    UU(M, 0.), GRADF(M, 0.), DSRCH(M, 0.), HESSF(M*(M+1)/2, 0.),
    f0_iters(n_rel_change_iters);
    
    std::vector<int> IYFREE(M, 0);
    std::vector<bool> eval_grads(M, false);
    
    Real
    ALBEFA  = 0.1,
    GHINIT  = _initial_rel_step,
    GHDECR  = _asymptote_reduction,
    GHINCR  = _asymptote_expansion,
    F0VAL   = 0.,
    F0NEW   = 0.,
    F0APP   = 0.,
    RAA0    = 0.,
    Z       = 0.,
    GEPS    =_feval->tolerance();
    
    
    /*C********+*********+*********+*********+*********+*********+*********+
     C
     C  The meaning of some of the scalars and vectors in the program:
     C
     C     N  = Complex of variables x_j in the problem.
     C     M  = Complex of constraints in the problem (not including
     C          the simple upper and lower bounds on the variables).
     C ALBEFA = Relative spacing between asymptote and mode limit. Lower value
     C          will cause the move limit (alpha,beta) to move closer to asymptote
     C          values (l, u).
     C GHINIT = Initial asymptote setting. For the first two iterations the
     C          asymptotes (l, u) are defined based on offsets from the design
     C          point as this fraction of the design variable bounds, ie.
     C              l_j   =   x_j^k  - GHINIT * (x_j^max - x_j^min)
     C              u_j   =   x_j^k  + GHINIT * (x_j^max - x_j^min)
     C GHDECR = Fraction by which the asymptote is reduced for oscillating
     C          changes in design variables based on three consecutive iterations
     C GHINCR = Fraction by which the asymptote is increased for non-oscillating
     C          changes in design variables based on three consecutive iterations
     C INNMAX = Maximal number of inner iterations within each outer iter.
     C          A reasonable choice is INNMAX=10.
     C  ITER  = Current outer iteration number ( =1 the first iteration).
     C  GEPS  = Tolerance parameter for the constraints.
     C          (Used in the termination criteria for the subproblem.)
     C
     C   XVAL(j) = Current value of the variable x_j.
     C  XOLD1(j) = Value of the variable x_j one iteration ago.
     C  XOLD2(j) = Value of the variable x_j two iterations ago.
     C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
     C   XMIN(j) = Original lower bound for the variable x_j.
     C   XMAX(j) = Original upper bound for the variable x_j.
     C   XLOW(j) = Value of the lower asymptot l_j.
     C   XUPP(j) = Value of the upper asymptot u_j.
     C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
     C   BETA(j) = Upper bound for x_j in the MMA subproblem.
     C    F0VAL  = Value of the objective function f_0(x)
     C   FVAL(i) = Value of the i:th constraint function f_i(x).
     C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
     C   FMAX(i) = Right hand side of the i:th constraint.
     C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
     C             where k = (j-1)*M + i.
     C      P(k) = Coefficient p_ij in the MMA subproblem, where
     C             k = (j-1)*M + i.
     C      Q(k) = Coefficient q_ij in the MMA subproblem, where
     C             k = (j-1)*M + i.
     C     P0(j) = Coefficient p_0j in the MMA subproblem.
     C     Q0(j) = Coefficient q_0j in the MMA subproblem.
     C      B(i) = Right hand side b_i in the MMA subproblem.
     C    F0APP  = Value of the approximating objective function
     C             at the optimal soultion of the MMA subproblem.
     C   FAPP(i) = Value of the approximating i:th constraint function
     C             at the optimal soultion of the MMA subproblem.
     C    RAA0   = Parameter raa_0 in the MMA subproblem.
     C    RAA(i) = Parameter raa_i in the MMA subproblem.
     C      Y(i) = Value of the "artificial" variable y_i.
     C      Z    = Value of the "minimax" variable z.
     C      A(i) = Coefficient a_i for the variable z.
     C      C(i) = Coefficient c_i for the variable y_i.
     C   ULAM(i) = Value of the dual variable lambda_i.
     C  GRADF(i) = Gradient component of the dual objective function.
     C  DSRCH(i) = Search direction component in the dual subproblem.
     C  HESSF(k) = Hessian matrix component of the dual function.
     C IYFREE(i) = 0 for dual variables which are fixed to zero in
     C               the current subspace of the dual subproblem,
     C           = 1 for dual variables which are "free" in
     C               the current subspace of the dual subproblem.
     C
     C********+*********+*********+*********+*********+*********+*********+*/
    
    
    /*
     *  The USER should now give values to the parameters
     *  M, N, GEPS, XVAL (starting point),
     *  XMIN, XMAX, FMAX, A and C.
     */
    // _initi(M,N,GEPS,XVAL,XMIN,XMAX,FMAX,A,C);
    // Assumed:  FMAX == A
    _feval->_init_dvar_wrapper(XVAL, XMIN, XMAX);
    // set the value of C[i] to be very large numbers
    Real max_x = 0.;
    for (unsigned int i=0; i<N; i++)
        if (max_x < fabs(XVAL[i]))
            max_x = fabs(XVAL[i]);
    
    int INNMAX=_max_inner_iters, ITER=0, ITE=0, INNER=0, ICONSE=0;
    /*
     *  The outer iterative process starts.
     */
    bool terminate = false, inner_terminate=false;
    while (!terminate) {
        
        std::fill(C.begin(), C.end(), std::max(1.e0*max_x, _constr_penalty));
        GHINIT  = _initial_rel_step,
        GHDECR  = _asymptote_reduction,
        GHINCR  = _asymptote_expansion,

        ITER=ITER+1;
        ITE=ITE+1;
        /*
         *  The USER should now calculate function values and gradients
         *  at XVAL. The result should be put in F0VAL,DF0DX,FVAL,DFDX.
         */
        std::fill(eval_grads.begin(), eval_grads.end(), true);
        _feval->_evaluate_wrapper(XVAL,
                                  F0VAL, true, DF0DX,
                                  FVAL, eval_grads, DFDX);
        if (ITER == 1)
            // output the very first iteration
            _feval->_output_wrapper(0, XVAL, F0VAL, FVAL, true);
        
        /*
         *  RAA0,RAA,XLOW,XUPP,ALFA and BETA are calculated.
         */
        raasta_(&M, &N, &RAA0, &RAA[0], &XMIN[0], &XMAX[0], &DF0DX[0], &DFDX[0]);
        asympg_(&ITER, &M, &N, &ALBEFA, &GHINIT, &GHDECR, &GHINCR,
                &XVAL[0], &XMIN[0], &XMAX[0], &XOLD1[0], &XOLD2[0],
                &XLOW[0], &XUPP[0], &ALFA[0], &BETA[0]);
        /*
         *  The inner iterative process starts.
         */
        
        // write the asymptote data for the inneriterations
        _output_iteration_data(ITER, XVAL, XMIN, XMAX, XLOW, XUPP, ALFA, BETA);

        INNER=0;
        inner_terminate = false;
        while (!inner_terminate) {
            
            /*
             *  The subproblem is generated and solved.
             */
            mmasug_(&ITER, &M, &N, &GEPS, &IYFREE[0], &XVAL[0], &XMMA[0],
                    &XMIN[0], &XMAX[0], &XLOW[0], &XUPP[0], &ALFA[0], &BETA[0],
                    &A[0], &B[0], &C[0], &Y[0], &Z, &RAA0, &RAA[0], &ULAM[0],
                    &F0VAL, &FVAL[0], &F0APP, &FAPP[0], &FMAX[0], &DF0DX[0], &DFDX[0],
                    &P[0], &Q[0], &P0[0], &Q0[0], &UU[0], &GRADF[0], &DSRCH[0], &HESSF[0]);
            /*
             *  The USER should now calculate function values at XMMA.
             *  The result should be put in F0NEW and FNEW.
             */
            std::fill(eval_grads.begin(), eval_grads.end(), false);
            _feval->_evaluate_wrapper(XMMA,
                                      F0NEW, false, DF0DX,
                                      FNEW, eval_grads, DFDX);
            
            ///////////////////////////////////////////////////////////////
            // if the solution is poor, backtrack
            std::vector<Real> XMMA_new(XMMA);
            Real frac = 0.02;
            while (FNEW[0] > 1.e2) {
                libMesh::out << "*** Backtracking: frac = "
                << frac
                << "  constr: " << FNEW[0]
                << std::endl;
                for (unsigned int i=0; i<XMMA.size(); i++)
                    XMMA_new[i] = XOLD1[i] + frac*(XMMA[i]-XOLD1[i]);
                
                _feval->_evaluate_wrapper(XMMA_new,
                                          F0NEW, false, DF0DX,
                                          FNEW, eval_grads, DFDX);
                frac *= frac;
            }
            for (unsigned int i=0; i<XMMA.size(); i++)
                XMMA[i] = XMMA_new[i];
            ///////////////////////////////////////////////////////////////

            if (INNER >= INNMAX) {
                libMesh::out
                << "** Max Inner Iter Reached: Terminating! Inner Iter = "
                << INNER << std::endl;
                inner_terminate = true;
            }
            else {
                /*
                 *  It is checked if the approximations were conservative.
                 */
                conser_( &M, &ICONSE, &GEPS, &F0NEW, &F0APP, &FNEW[0], &FAPP[0]);
                if (ICONSE == 1) {
                    libMesh::out
                    << "** Conservative Solution: Terminating! Inner Iter = "
                    << INNER << std::endl;
                    inner_terminate = true;
                }
                else {
                    /*
                     *  The approximations were not conservative, so RAA0 and RAA
                     *  are updated and one more inner iteration is started.
                     */
                    INNER=INNER+1;
                    raaupd_( &M, &N, &GEPS, &XMMA[0], &XVAL[0],
                            &XMIN[0], &XMAX[0], &XLOW[0], &XUPP[0],
                            &F0NEW, &FNEW[0], &F0APP, &FAPP[0], &RAA0, &RAA[0]);
                }
            }
        }
        
        /*
         *  The inner iterative process has terminated, which means
         *  that an outer iteration has been completed.
         *  The variables are updated so that XVAL stands for the new
         *  outer iteration point. The fuction values are also updated.
         */
        xupdat_( &N, &ITER, &XMMA[0], &XVAL[0], &XOLD1[0], &XOLD2[0]);
        fupdat_( &M, &F0NEW, &FNEW[0], &F0VAL, &FVAL[0]);
        /*
         *  The USER may now write the current solution.
         */
        _feval->_output_wrapper(ITER, XVAL, F0VAL, FVAL, true);
        f0_iters[(ITE-1)%n_rel_change_iters] = F0VAL;
        
        /*
         *  One more outer iteration is started as long as
         *  ITE is less than MAXITE:
         */
        if (ITE == _feval->max_iters()) {
            libMesh::out
            << "GCMMA: Reached maximum iterations, terminating! "
            << std::endl;
            terminate = true;
        }
        
        // relative change in objective
        bool rel_change_conv = true;
        Real f0_curr = f0_iters[n_rel_change_iters-1];
        
        for (unsigned int i=0; i<n_rel_change_iters-1; i++) {
            if (f0_curr > sqrt(GEPS))
                rel_change_conv = (rel_change_conv &&
                                   fabs(f0_iters[i]-f0_curr)/fabs(f0_curr) < GEPS);
            else
                rel_change_conv = (rel_change_conv &&
                                   fabs(f0_iters[i]-f0_curr) < GEPS);
        }
        if (rel_change_conv) {
            libMesh::out
            << "GCMMA: Converged relative change tolerance, terminating! "
            << std::endl;
            terminate = true;
        }
        
    }
    
#endif //MAST_ENABLE_GCMMA == 1
}


void
MAST::GCMMAOptimizationInterface::_output_iteration_data(unsigned int i,
                                                         const std::vector<Real>& XVAL,
                                                         const std::vector<Real>& XMIN,
                                                         const std::vector<Real>& XMAX,
                                                         const std::vector<Real>& XLOW,
                                                         const std::vector<Real>& XUPP,
                                                         const std::vector<Real>& ALFA,
                                                         const std::vector<Real>& BETA) {

    libmesh_assert(_feval);
    libmesh_assert_equal_to(XVAL.size(), _feval->n_vars());
    libmesh_assert_equal_to(XMIN.size(), _feval->n_vars());
    libmesh_assert_equal_to(XMAX.size(), _feval->n_vars());
    libmesh_assert_equal_to(XLOW.size(), _feval->n_vars());
    libmesh_assert_equal_to(XUPP.size(), _feval->n_vars());
    libmesh_assert_equal_to(ALFA.size(), _feval->n_vars());
    libmesh_assert_equal_to(BETA.size(), _feval->n_vars());

    /*libMesh::out
    << "****************************************************\n"
    << "             GCMMA: ASYMPTOTE DATA                  \n"
    << "****************************************************\n"
    << std::setw(5) << "Iter: " << i << std::endl
    << std::setw(5) << "DV"
    << std::setw(20) << "XMIN"
    << std::setw(20) << "XLOW"
    << std::setw(20) << "ALFA"
    << std::setw(20) << "X"
    << std::setw(20) << "BETA"
    << std::setw(20) << "XUP"
    << std::setw(20) << "XMAX" << std::endl;

    for (unsigned int j=0; j<_feval->n_vars(); j++)
        libMesh::out
        << std::setw(5) << j
        << std::setw(20) << XMIN[j]
        << std::setw(20) << XLOW[j]
        << std::setw(20) << ALFA[j]
        << std::setw(20) << XVAL[j]
        << std::setw(20) << BETA[j]
        << std::setw(20) << XUPP[j]
        << std::setw(20) << XMAX[j] << std::endl;
     */
}

