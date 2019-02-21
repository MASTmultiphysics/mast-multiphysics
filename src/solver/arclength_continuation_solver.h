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

#ifndef __mast__arclength_continuation_solver_h__
#define __mast__arclength_continuation_solver_h__

// MAST includes
#include "solver/continuation_solver_base.h"


namespace MAST {
    
    /*!
     *   constraint equation is defined along the path \f$ s \f$ as
     *   \f[ g(X, p, ds) = (x-x0)  dx/ds +   (p-p0) dp/ds  -  ds = 0, \f]
     *   where, \f$ X \f$ is the solution, \f$ p \f$ is the load parameter, and
     *   \f$ ds \f$ is the chord length.
     */
    class ArclengthContinuationSolver:
    public MAST::ContinuationSolverBase {
        
    public:
        
        ArclengthContinuationSolver();
        
        virtual ~ArclengthContinuationSolver();
        
        /*!
         *   sets the arc length using a nonlinear solution using a
         *   step \p dp.
         */
        virtual void initialize(Real dp);

    protected:
        
        virtual void
        _solve_NR_iterate(libMesh::NumericVector<Real>       &X,
                          MAST::Parameter                    &p);
        
        virtual void
        _dXdp(const libMesh::NumericVector<Real> &X,
              const MAST::Parameter              &p,
              libMesh::NumericVector<Real>       &dfdp,
              libMesh::NumericVector<Real>       &dXdp);
        
        /*!
         * \f[
         *    g(X, p, ds) =  (x-x0)  dX/ds +   (p-p0) dp/ds  -  ds = 0
         * \f]
         */
        virtual Real
        _g(const libMesh::NumericVector<Real> &X,
           const MAST::Parameter              &p);
        
        /*!
         * \f{eqnarray*}{
         *    g(X, p, ds) & = & (x-x0)  dX/ds +   (p-p0) dp/ds  -  ds = 0 \\
         *    dg/dp       & = & dp/ds \\
         *    dg/dX       & = & dX/ds
         * \f}
         */
        void
        _g(const libMesh::NumericVector<Real> &X,
           const MAST::Parameter              &p,
           libMesh::NumericVector<Real>       &dfdp,
           libMesh::NumericVector<Real>       &dXdp,
           Real                               &g,
           Real                               &dgdp,
           libMesh::NumericVector<Real>       *dgdX);
    };
}

#endif // __mast__arclength_continuation_solver_h__

