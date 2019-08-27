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
#include "optimization/npsol_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "base/mast_config.h"

extern "C" {
    
    extern void sninit_(int* iPrint,
                        int* iSumm,
                        const char* cw,
                        int* lencw,
                        int* iw,
                        int* leniw,
                        double* rw,
                        int* lenrw);

    extern void sninitf_(const char* printfile,
                         const char* summary_file,
                         int* iPrint,
                         int* iSumm,
                         const char* cw,
                         int* lencw,
                         int* iw,
                         int* leniw,
                         double* rw,
                         int* lenrw);

    extern void snspec_(int* iSpecs,
                        int* INFO,
                        const char* cw,
                        int* lencw,
                        int* iw,
                        int* leniw,
                        double* rw,
                        int* lenrw);

    extern void snspecf_(const char* specsfile,
                         int* INFO,
                         const char* cw,
                         int* lencw,
                         int* iw,
                         int* leniw,
                         double* rw,
                         int* lenrw);

    extern void npopt_(int*    n,
                       int*    nclin,
                       int*    ncnln,
                       int*    ldA,
                       int*    ldgg,
                       int*    ldH,
                       double* A,
                       double* bl,
                       double* bu,
                       void    (*)(int*    mode,
                                   int*    ncnln,
                                   int*    n,
                                   int*    ldJ,
                                   int*    needc,
                                   double* x,
                                   double* c,
                                   double* cJac,
                                   int*    nstate),
                       void    (*)(int*    mode,
                                   int*    n,
                                   double* x,
                                   double* f,
                                   double* g,
                                   int*    nstate),
                       int*    INFO,
                       int*    majIts,
                       int*    iState,
                       double* fCon,
                       double* gCon,
                       double* cMul,
                       double* fObj,
                       double* gObj,
                       double* Hess,
                       double* x,
                       int*    iw,
                       int*    leniw,
                       double* re,
                       int*    lenrw);
    
    extern void npoptn_(const char*, int );
}




MAST::NPSOLOptimizationInterface::NPSOLOptimizationInterface():
MAST::OptimizationInterface(),
_funobj(nullptr),
_funcon(nullptr) {
    
#if MAST_ENABLE_SNOPT == 0
    libmesh_error_msg("MAST configured without SNOPT support.");
#endif
    
    _exit_message[0] = "Finished successfully";
    _exit_message[10] = " The problem appears to be infeasible 20 The problem appears to be unbounded 30 Resource limit error";
    _exit_message[40] = " Terminated after numerical difficulties 50 Error in the user-supplied functions";
    _exit_message[60] = " Undefined user-supplied functions";
    _exit_message[70] = " User requested termination";
    _exit_message[80] = " Insufficient storage allocated";
    _exit_message[90] = " Input arguments out of range";
    _exit_message[100] = " Finished successfully (associated with SNOPT auxiliary routines) 110 Errors while processing MPS data";
    _exit_message[120] = " Errors while estimating Jacobian structure";
    _exit_message[130] = " Errors while reading the Specs file";
    _exit_message[140] = " System error";
    
    _info_message[1]   = "optimality conditions satisfied";
    _info_message[2]   = "feasible point found";
    _info_message[3]   = "requested accuracy could not be achieved";
    _info_message[11]  = "infeasible linear constraints";
    _info_message[12]  = "infeasible linear equality constraints";
    _info_message[13]  = "nonlinear infeasibilities minimized";
    _info_message[14]  = "linear infeasibilities minimized";
    _info_message[15]  = "infeasible linear constraints in QP subproblem";
    _info_message[16]  = "infeasible nonelastic constraints";
    _info_message[21]  = "unbounded objective";
    _info_message[22]  = "constraint violation limit reached";
    _info_message[31]  = "iteration limit reached";
    _info_message[32]  = "major iteration limit reached";
    _info_message[33]  = "the superbasics limit is too small";
    _info_message[41]  = "current point cannot be improved";
    _info_message[42]  = "singular basis";
    _info_message[43]  = "cannot satisfy the general constraints";
    _info_message[44]  = "ill-conditioned null-space basis";
    _info_message[51]  = "incorrect objective derivatives";
    _info_message[52]  = "incorrect constraint derivatives";
    _info_message[56]  = "irregular or badly scaled problem functions";
    _info_message[61]  = "undefined function at the first feasible point";
    _info_message[62]  = "undefined function at the initial point";
    _info_message[63]  = "unable to proceed into undefined region";
    _info_message[71]  = "terminated during function evaluation";
    _info_message[72]  = "terminated during constraint evaluation";
    _info_message[73]  = "terminated during objective evaluation";
    _info_message[74]  = "terminated from monitor routine";
    _info_message[81]  = "work arrays must have at least 500 elements";
    _info_message[82]  = "not enough character storage";
    _info_message[83]  = "not enough integer storage";
    _info_message[84]  = "not enough real storage";
    _info_message[91]  = "invalid input argument";
    _info_message[92]  = "basis file dimensions do not match this problem";
    _info_message[141] = "wrong number of basic variables";
    _info_message[142] = "error in basis package";
}



void
MAST::NPSOLOptimizationInterface::
attach_function_evaluation_object (MAST::FunctionEvaluation& feval) {
    
#if MAST_ENABLE_SNOPT == 1
    MAST::OptimizationInterface::attach_function_evaluation_object(feval);
    
    // make sure that these pointers haven't already been provided
    libmesh_assert(_funobj == nullptr);
    libmesh_assert(_funcon == nullptr);
    
    _funobj = feval.get_objective_evaluation_function();
    _funcon = feval.get_constraint_evaluation_function();
#endif
}



void
MAST::NPSOLOptimizationInterface::optimize() {
    
#if MAST_ENABLE_SNOPT == 1
    // make sure that functions have been provided
    libmesh_assert(_funobj);
    libmesh_assert(_funcon);
    
    _feval->sanitize_parallel();

    int
    iPrint = 9,
    iSpec  = 4,
    iSumm  = 6,
    lencw  = 500,
    N      =  _feval->n_vars(),
    NCLIN  =  0,
    NCNLN  =  _feval->n_eq()+_feval->n_ineq(),
    NCTOTL =  N+NCLIN+NCNLN,
    LDA    =  std::max(NCLIN, 1),
    LDJ    =  std::max(NCNLN, 1),
    LDR    =  N,
    INFORM =  0,           // on exit: Reports result of call to NPSOL
                           // < 0 either funobj or funcon has set this to -ve
                           // 0 => converged to point x
                           // 1 => x satisfies optimality conditions, but sequence of iterates has not converged
                           // 2 => Linear constraints and bounds cannot be satisfied. No feasible solution
                           // 3 => Nonlinear constraints and bounds cannot be satisfied. No feasible solution
                           // 4 => Major iter limit was reached
                           // 6 => x does not satisfy first-order optimality to required accuracy
                           // 7 => function derivatives seem to be incorrect
                           // 9 => input parameter invalid
    ITER   = 0,            // iter count
    LENIW  = std::max(1200*(NCTOTL+N) ,1000),
    LENW   = std::max(2400*(NCTOTL+N), 1000);
    
    Real
    F      =  0.;          // on exit: final objective

    std::vector<int>
    IW      (LENIW,  0),
    ISTATE  (NCTOTL, 0);    // status of constraints l <= r(x) <= u,
                            // -2 => lower bound is violated by more than delta
                            // -1 => upper bound is violated by more than delta
                            // 0  => both bounds are satisfied by more than delta
                            // 1  => lower bound is active (to within delta)
                            // 2  => upper bound is active (to within delta)
                            // 3  => boundars are equal and equality constraint is satisfied
    
    std::vector<Real>
    A       (LDA,    0.),   // this is used for liear constraints, not currently handled
    BL      (NCTOTL, 0.),
    BU      (NCTOTL, 0.),
    C       (NCNLN,  0.),   // on exit: nonlinear constraints
    CJAC    (LDJ* N, 0.),   //
                            // on exit: CJAC(i,j) is the partial derivative of ith nonlinear constraint
    CLAMBDA (NCTOTL, 0.),   // on entry: need not be initialized for cold start
                            // on exit: QP multiplier from the QP subproblem, >=0 if istate(j)=1, <0 if istate(j)=2
    G       (N,      0.),   // on exit: objective gradient
    R       (LDR*N,  0.),   // on entry: need not be initialized if called with Cold Statrt
                            // on exit: information about Hessian, if Hessian=Yes, R is upper Cholesky factor of approx H
    X       (N,      0.),   // on entry: initial point
                            // on exit: final estimate of solution
    W       (LENW,   0.),   // workspace
    xmin    (N,      0.),
    xmax    (N,      0.);
    
    
    // now setup the lower and upper limits for the variables and constraints
    _feval->_init_dvar_wrapper(X, xmin, xmax);
    for (unsigned int i=0; i<N; i++) {
        BL[i] = xmin[i];
        BU[i] = xmax[i];
    }
    
    // all constraints are assumed to be g_i(x) <= 0, so that the upper
    // bound is 0 and lower bound is -infinity
    for (unsigned int i=0; i<NCNLN; i++) {
        BL[i+N] = -1.e20;
        BU[i+N] =     0.;
    }
    
    std::string cw(lencw*8,' ');
//    nm = "List";
//    npoptn_(nm.c_str(), (int)nm.length());
//    nm = "Verify level 3";
//    npoptn_(nm.c_str(), (int)nm.length());
    
    sninit_(&iPrint,
            &iSumm,
            cw.c_str(),
            &lencw,
            &IW[0],
            &LENIW,
            &W[0],
            &LENW);
    
    snspec_(&iSpec, &INFORM, cw.c_str(), &lencw, &IW[0], &LENIW, &W[0], &LENW);
    
    if (INFORM != 101) libmesh_error();
    
    npopt_(&N,
           &NCLIN,
           &NCNLN,
           &LDA,
           &LDJ,
           &LDR,
           &A[0],
           &BL[0],
           &BU[0],
           _funcon,
           _funobj,
           &INFORM,
           &ITER,
           &ISTATE[0],
           &C[0],
           &CJAC[0],
           &CLAMBDA[0],
           &F,
           &G[0],
           &R[0],
           &X[0],
           &IW[0],
           &LENIW,
           &W[0],
           &LENW);
    
    _print_termination_message(INFORM);
    
#endif // MAST_ENABLE_SNOPT 1
}


void
MAST::NPSOLOptimizationInterface::_print_termination_message(const int INFORM) {
    
    int
    exit = (INFORM%10)*10; // remove the last signifncant digit
    
    libmesh_assert(_exit_message.count(exit));
    libmesh_assert(_info_message.count(INFORM));
    
    libMesh::out
    << "SNOPT EXIT :" << exit << std::endl
    << "   " << _exit_message[exit] << std::endl
    << "      INFO :" << INFORM << std::endl
    << "   " << _info_message[INFORM] << std::endl << std::endl;
}


