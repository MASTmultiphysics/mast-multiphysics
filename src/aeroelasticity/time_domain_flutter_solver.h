/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

#ifndef __mast__time_domain_flutter_solver_h__
#define __mast__time_domain_flutter_solver_h__

// C++ includes
#include <memory>


// MAST includes
#include "aeroelasticity/flutter_solver_base.h"



namespace MAST {
    
    // Forward declerations
    class TimeDomainFlutterSolution;
    
    
    /*!
     *   This implements a solver for a single parameter instability
     *   problem, for example a flutter solver where flight speed is the 
     *   primary parameter.
     */
    class TimeDomainFlutterSolver:
    public MAST::FlutterSolverBase {
      
    public:
        
        /*!
         *   defalut constructor
         */
        TimeDomainFlutterSolver();
        
        
        virtual ~TimeDomainFlutterSolver();
        
        
        
        /*!
         *   clears the solution and other data from this solver
         */
        virtual void clear();

        
        /*!
         *   clears the solutions stored from a previous analysis.
         */
        virtual void clear_solutions();

        
        /*!
         *    initializes the data structres for a flutter solution.
         */
        void initialize(MAST::Parameter&                            velocity_param,
                        Real                                        V_lower,
                        Real                                        V_upper,
                        unsigned int                                n_V_divs,
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
         *   This root starts with the lower velocity and increments the speed
         *   till a single unstable root is identified.
         */
        virtual std::pair<bool, MAST::FlutterRootBase*>
        analyze_and_find_critical_root_without_tracking(const Real g_tol,
                                                        const unsigned int n_iters);
        
        
        
        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params. If the base solution has a sensitivity
         *   with respect to the parameter, then that should be provided 
         *   through \par dXdp. The sensitivity solution also requires
         *   sensitivity of the eigenvalue wrt velocity, which is 
         *   defined as a parameter. Hence, the sensitivity of the 
         *   static solution is also required wrt the velocity parameter. 
         *   If \par dXdV is \par nullptr, then zero value is assumed. 
         */
        virtual void
        calculate_sensitivity(MAST::FlutterRootBase& root,
                              const libMesh::ParameterVector& params,
                              const unsigned int i,
                              libMesh::NumericVector<Real>* dXdp = nullptr,
                              libMesh::NumericVector<Real>* dXdV = nullptr);
        
        
        /*!
         *   Prints the sorted roots to the \par output
         */
        virtual void print_sorted_roots();
        
        
        /*!
         *   Prints the crossover points output. If no pointer to output is given
         *   then the output defined by set_output_file() is used.
         */
        virtual void print_crossover_points();
        
        
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
         *   pointer is nullptr, then no sorting is performed
         */
        virtual std::auto_ptr<MAST::TimeDomainFlutterSolution>
        _analyze(const Real v_ref,
                 const MAST::FlutterSolutionBase* prev_sol=nullptr);
        
        
        
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
         *    matrices for specified flight velocity \par U_inf.
         */
        void _initialize_matrices(Real U_inf,
                                  RealMatrixX& A,
                                  RealMatrixX& B);

        
        /*!
         *    Assembles the reduced order system structural and aerodynmaic
         *    matrices for specified flight velocity \par U_inf.
         */
        void
        _initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                                 const unsigned int i,
                                                 const libMesh::NumericVector<Real>& dXdp,
                                                 Real U_inf,
                                                 RealMatrixX& A,
                                                 RealMatrixX& B);

        
        /*!
         *   identifies all cross-over and divergence points from analyzed
         *   roots
         */
        virtual void _identify_crossover_points();


        /*!
         *  Parameter that define the velocity.
         */
        MAST::Parameter*                                _velocity_param;
        
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real>                           _V_range;


        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int                                    _n_V_divs;

        
        /*!
         *   map of velocity sorted flutter solutions
         */
        std::map<Real, MAST::FlutterSolutionBase*>      _flutter_solutions;
        
        /*!
         *   the map of flutter crossover points versus average velocity of the
         *   two bounding roots
         */
        std::multimap<Real, MAST::FlutterRootCrossoverBase*> _flutter_crossovers;

    };
}


#endif // __mast__time_domain_flutter_solver_h__

