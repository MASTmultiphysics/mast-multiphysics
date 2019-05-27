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

#ifndef __mast__pseudo_arclength_continuation_solver_h__
#define __mast__pseudo_arclength_continuation_solver_h__

// MAST includes
#include "solver/continuation_solver_base.h"


namespace MAST {
    
    /*!
     *   The constraint equation is defined along the path \f$ s \f$ as
     *   \f[ g(X, p, ds) = (Y - \tilde{Y} )^T t_1 = 0 , \f]
     *   where, \f$ Y = \{ X^T  p \}^T \f$, \f$ X \f$ is the solution,
     *   \f$ p \f$ is the load parameter, \f$\tilde{Y} \f$ is the predictor based on
     *   the search direction \f$ t_1 \f$. Given that the predictor is defined
     *   as \f$ \tilde{Y} = Y_0 - ds t_1 \f$, with \f$ ds \f$ as the step size,
     *   the constraint is rewritten as
     *   \f[ (Y - Y_0 - ds * t_1)^T t_1 = 0 , \f]
     *   or, assuming \f$ \| t_1 \| = 1 \f$,
     *   \f[ Y^T t_1 - Y_0^T t_1 - ds  = 0 . \f]
     *
     *   The search direction is evaluated based on:
     *   \f[   \left[ \begin{array}{cc} df/dx & df/dp \\
     *              t_0^X &  t_0^p \end{array} \right]
     *          \left\{ \begin{array}{c} t_1^X \\ t_1^p \end{array} \right\} =
     *         \left\{ \begin{array}{c} 0 \\ 1 \end{array} \right\}
     *    \f]
     *    and, \f$ t_1 \f$ is then scaled to unit norm.
     */
    class PseudoArclengthContinuationSolver:
    public MAST::ContinuationSolverBase {
        
    public:
        
        PseudoArclengthContinuationSolver();
        
        virtual ~PseudoArclengthContinuationSolver();
        
        /*!
         *   initializes the search direction using the specified load
         *   step \p dp.
         */
        virtual void initialize(Real dp);
        
        
    protected:
        
        virtual void
        _solve_NR_iterate(libMesh::NumericVector<Real>       &X,
                          MAST::Parameter                    &p);
        
        /*!
         *  updates \f$ t_1 \f$ for the current iterate \p X and \p p, and
         *  stores these values to \p _t0_X and \p _t0_p for next computation.
         */
        void
        _update_search_direction(const libMesh::NumericVector<Real> &X,
                                 const MAST::Parameter              &p,
                                 libMesh::SparseMatrix<Real>        &jac,
                                 /*libMesh::NumericVector<Real>       &dfdp,
                                 libMesh::NumericVector<Real>       &dXdp,*/
                                 libMesh::NumericVector<Real>       &t1_X,
                                 Real                               &t1_p);
        
        
        /*!
         * \f[
         *    g(X, p, ds) =  X_{scale} * (X-X0)  (dX/ds)_{scaled} +
         *                   p_{scale} * (p-p0)  (dp/ds)_{scaled}  -  ds = 0,
         * \f]
         *  where,
         *  \f$ t_1^{X} = (dX/ds)_{scaled} \f$ and
         *  \f$ t_1^{p} = (dp/ds)_{scaled} \f$.
         */
        virtual Real
        _g(const libMesh::NumericVector<Real> &X,
           const MAST::Parameter              &p);
        
        
        /*!
         * \f{eqnarray*}{
         *    g(X, p, ds) & = & (Y-Y_0)^T t_1 - ds = 0 \\
         *    Y           & = & \{ X_{scale} X^T, p_{scale} p\}^T \\
         *    dg/dp       & = & p_{scale} t_1^p \\
         *    dg/dX       & = & X_{scale} t_1^X,
         * \f}
         *  where \f$  \f$,
         *  \f$ t_1^{X} = (dX/ds)_{scaled} \f$ and
         *  \f$ t_1^{p} = (dp/ds)_{scaled} \f$.
         */
        void
        _g(const libMesh::NumericVector<Real> &X,
           const MAST::Parameter              &p,
           libMesh::NumericVector<Real>       &t1_X,
           Real                               &t1_p,
           Real                               &g);
        
        std::unique_ptr<libMesh::NumericVector<Real>>  _t0_X;
        Real                                           _t0_p;
    };
}


#endif // __mast__pseudo_arclength_continuation_solver_h__

