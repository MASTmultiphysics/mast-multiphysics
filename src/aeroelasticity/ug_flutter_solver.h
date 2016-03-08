/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__ug_flutter_solver_h__
#define __mast__ug_flutter_solver_h__


// MAST includes
#include "aeroelasticity/flutter_solver_base.h"


namespace MAST {
    
    
    /*!
     *   This implements a solver for a single parameter instability
     *   problem, for example a flutter solver where reduced frequency is the
     *   primary parameter.
     */
    class UGFlutterSolver:
    public MAST::FlutterSolverBase {
        
    public:
        
        /*!
         *   defalut constructor
         */
        UGFlutterSolver();
        
        
        virtual ~UGFlutterSolver();
        
        
        
        /*!
         *   clears the solution and other data from this solver
         */
        void clear();
        
        
        /*!
         *   clears the solutions stored from a previous analysis.
         */
        virtual void clear_solutions();
        
        
        /*!
         *    initializes the data structres for a flutter solution.
         */
        void initialize(MAST::Parameter&                            kr_param,
                        MAST::Parameter&                            bref_param,
                        Real                                        rho,
                        Real                                        kr_lower,
                        Real                                        kr_upper,
                        unsigned int                                n_kr_divs,
                        std::vector<libMesh::NumericVector<Real>*>& basis);
        
        
        
        /*!
         *    finds the number of critical points already identified in the
         *    procedure.
         */
        virtual unsigned int n_roots_found() const;
        
        
        /*!
         *   returns the \par n th root in terms of ascending velocity that is
         *   found by the solver
         */
        const MAST::FlutterRootBase& get_root(const unsigned int n) const;
        
                
        /*!
         *   Looks through the list of flutter cross-over points and iteratively
         *   zooms in to find the cross-over point. This should be called only
         *   after scan_for_roots() has been called. Potential cross-over points
         *   are sorted with increasing velocity, and this method will attempt to
         *   identify the next critical root in the order.
         */
        virtual std::pair<bool, MAST::FlutterRootBase*>
        find_next_root(const Real g_tol,
                       const unsigned int n_bisection_iters);
        
        
        
        /*!
         *   This method checks if the flutter root corresponding to the
         *   lowest velocity crossover has been calculated. If not, then it
         *   attempts to find that root using an iterative approach
         */
        virtual std::pair<bool, MAST::FlutterRootBase*>
        find_critical_root(const Real g_tol,
                           const unsigned int n_bisection_iters);
        
        
        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params. If the base solution has a sensitivity
         *   with respect to the parameter, then that should be provided
         *   through \par dXdp. The sensitivity solution also requires
         *   sensitivity of the eigenvalue wrt velocity, which is
         *   defined as a parameter. Hence, the sensitivity of the
         *   static solution is also required wrt the velocity parameter.
         *   If \par dXdV is \par NULL, then zero value is assumed.
         */
        virtual void
        calculate_sensitivity(MAST::FlutterRootBase& root,
                              const libMesh::ParameterVector& params,
                              const unsigned int i,
                              libMesh::NumericVector<Real>* dXdp = NULL,
                              libMesh::NumericVector<Real>* dXdV = NULL);
        
        
        /*!
         *   Prints the sorted roots to the \par output
         */
        virtual void print_sorted_roots(std::ostream* output = NULL);
        
        
        /*!
         *   Prints the crossover points output. If no pointer to output is given
         *   then the output defined by set_output_file() is used.
         */
        virtual void print_crossover_points(std::ostream* output = NULL);
        
        
        /*!
         *   Scans for flutter roots in the analyzed points, and identified the
         *   divergence (if k_red = 0. is specified) and flutter crossover points.
         *   The roots are organized in terms of increasing velocity.
         */
        virtual void scan_for_roots();
        
        
    protected:
        
        
        /*!
         *   performs an eigensolution at the specified reference value, and
         *   sort the roots based on the provided solution pointer. If the
         *   pointer is NULL, then no sorting is performed
         */
        virtual std::auto_ptr<MAST::FlutterSolutionBase>
        _analyze(const Real kr_ref,
                 const MAST::FlutterSolutionBase* prev_sol=NULL);
        
        
        
        /*!
         *    bisection method search
         */
        virtual std::pair<bool, MAST::FlutterSolutionBase*>
        _bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                          MAST::FlutterSolutionBase*>& ref_sol_range,
                          const unsigned int root_num,
                          const Real g_tol,
                          const unsigned int max_iters);
        
        
        /*!
         *    Assembles the reduced order system structural and aerodynmaic
         *    matrices for specified reduced freq \par kr.
         */
        void _initialize_matrices(Real kr,
                                  ComplexMatrixX& A,
                                  ComplexMatrixX& B);
        
        
        /*!
         *    Assembles the reduced order system structural and aerodynmaic
         *    matrices for specified flight velocity \par U_inf.
         */
        void
        _initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                                 const unsigned int i,
                                                 const libMesh::NumericVector<Real>& dXdp,
                                                 Real U_inf,
                                                 ComplexMatrixX& A,
                                                 ComplexMatrixX& B);
        
        
        /*!
         *   identifies all cross-over and divergence points from analyzed
         *   roots
         */
        virtual void _identify_crossover_points();
        
        
        /*!
         *  Parameter that define the reduced frequency.
         */
        MAST::Parameter*                                _kr_param;
        
        
        /*!
         *    reference chord
         */
        MAST::Parameter*                                _bref_param;

        
        /*!
         *    flight density
         */
        Real                                            _rho;
        
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real>                           _kr_range;
        
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int                                    _n_kr_divs;
        
        
        /*!
         *   map of kr sorted flutter solutions
         */
        std::map<Real, MAST::FlutterSolutionBase*>      _flutter_solutions;
        
        /*!
         *   the map of flutter crossover points versus average kr of the
         *   two bounding roots
         */
        std::multimap<Real, MAST::FlutterRootCrossoverBase*> _flutter_crossovers;
        
    };
}


#endif // __mast__ug_flutter_solver_h__
