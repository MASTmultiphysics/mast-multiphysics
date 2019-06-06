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

#ifndef __mast__continuation_solver_base_h__
#define __mast__continuation_solver_base_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


namespace MAST {

    // Forward declerations
    class AssemblyElemOperations;
    class AssemblyBase;
    class Parameter;
    
    /*!
     *    the equation set is:
     *    \f$
     *        \left\{ \begin{array}{c} f(x, p) \\ g(x, p) \end{array} \right\} =
     *        \left\{ \begin{array}{c} 0       \\ 0       \end{array} \right\}
     *     \f$
     *    the N-R updates are calculated such that
     *    \f[
     *         \left[ \begin{array}{cc}
     *            df/dx & df/dp \\ dg/dx  &  dg/dp \end{array}\right]
     *          \left\{ \begin{array}{c} dx \\ dp \end{array} \right\}  =
     *          - \left\{ \begin{array}{c} f(x0, p0) \\ g(x0, p0) \end{array} \right\}
     *    \f]
     *    This equation is solved using Schur-factorization so that the
     *    disciplinary linear solver can be used.
     */
    class ContinuationSolverBase {
        
    public:
        
        ContinuationSolverBase();
        
        
        virtual ~ContinuationSolverBase();
        
        /*!
         *   sets the assembly object for this solver
         */
        void
        set_assembly_and_load_parameter(MAST::AssemblyElemOperations&   elem_ops,
                                        MAST::AssemblyBase&             assembly,
                                        MAST::Parameter&                p);
        
        /*!
         *   clears the assembly object from this solver
         */
        void clear_assembly_and_load_parameters();
        
        /*!
         *   initializes the data structure based on initial load step \p dp.
         *   must be called before solve().
         */
        virtual void initialize(Real dp) =0;
        
        /*!
         *   solves for the next load step
         */
        virtual void solve();
        
        /*!
         *  Maximum number of Newton-Raphson iterations for the solver.
         *  Default is 20.
         */
        unsigned int
        max_it;
        
        /*!
         *  Absolute tolerance for the solver. Default is 1.e-8;
         */
        Real
        abs_tol;

        /*!
         *  Relative tolerance for the solver. Default is 1.e-8;
         */
        Real
        rel_tol;
    
        /*!
         *   arc length that the solver is required to satisfy for the update.
         */
        Real
        arc_length;
        
        /*!
         *   minimum step size allowed with adaptivity
         */
        Real
        min_step;
        
        /*!
         *   maximum step size allowed with adaptivity
         */
        Real
        max_step;
        
        /*!
         *   exponent used in step size update.
         */
        Real
        step_size_change_exponent;

        /*!
         *   desired N-R iterations per load-step. Step-size is chanegd
         *   if actual number of N-R iterates is different from this value
         *   using the expression \f$ \Delta s \left(\frac{N_{desired}}{N_{actual}}\right)^p \f$,
         *   where, p is the exponent. 
         */
        unsigned int
        step_desired_iters;
        
        /*!
         *   flag to use Schur-factorizaiton (default) or monolithic solver
         */
        bool schur_factorization;
        
    protected:

        virtual void
        _solve_NR_iterate(libMesh::NumericVector<Real>       &X,
                          MAST::Parameter                    &p) = 0;

        /*!
         *   solves for the linear system of equation as a monolithic system
         *    \f[
         *         \left[ \begin{array}{cc}
         *            df/dx & df/dp \\ dg/dx  &  dg/dp \end{array}\right]
         *          \left\{ \begin{array}{c} dx \\ dp \end{array} \right\}  =
         *          - \left\{ \begin{array}{c} f \\ g \end{array} \right\}
         *    \f]
         *   \p dX and \p dp are returned from the solution
         */
        void
        _solve(const libMesh::NumericVector<Real>  &X,
               const MAST::Parameter               &p,
               libMesh::NumericVector<Real>        &f,
               bool                                update_f,
               libMesh::NumericVector<Real>        &dfdp,
               bool                                update_dfdp,
               const libMesh::NumericVector<Real>  &dgdX,
               const Real                          dgdp,
               const Real                          g,
               libMesh::NumericVector<Real>        &dX,
               Real                                &dp);


        /*!
         *   solves for the linear system of equation using Schur factorization.
         *    \f[
         *         \left[ \begin{array}{cc}
         *            df/dx & df/dp \\ dg/dx  &  dg/dp \end{array}\right]
         *          \left\{ \begin{array}{c} dx \\ dp \end{array} \right\}  =
         *          - \left\{ \begin{array}{c} f \\ g \end{array} \right\}
         *    \f]
         *   \p dX and \p dp are returned from the solution
         */
        void
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
                                   Real                                &dp);

        /*!
         *   @return the norm of residual at given solution and
         *   load parameter.
         */
        Real _res_norm(const libMesh::NumericVector<Real>                  &X,
                       const MAST::Parameter                               &p);
        
        virtual Real
        _g(const libMesh::NumericVector<Real> &X,
           const MAST::Parameter              &p) = 0;

        /*!
         *   method saves any data for possible resuse if the solution step is restarted
         */
        virtual void _save_iteration_data() = 0;

        /*!
         *   method resets any data if a soltion step is restarted
         */
        virtual void _reset_iterations() = 0;
        
        bool                           _initialized;
        
        MAST::AssemblyElemOperations   *_elem_ops;
        MAST::AssemblyBase             *_assembly;
        MAST::Parameter                *_p;

        Real
        _p0,
        _X_scale,
        _p_scale;
        
        std::unique_ptr<libMesh::NumericVector<Real>>
        _X0;
    };
}


#endif // __mast__continuation_solver_base_h__
