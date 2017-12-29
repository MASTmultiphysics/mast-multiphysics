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

#ifndef __mast__pk_flutter_solver_h__
#define __mast__pk_flutter_solver_h__

// C++ includes
#include <memory>


// MAST includes
#include "aeroelasticity/flutter_solver_base.h"


namespace MAST {
    
    class PKFlutterSolver: public MAST::FlutterSolverBase
    {
    public:
        PKFlutterSolver();
        
        
        virtual ~PKFlutterSolver();
        

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
        void initialize(MAST::Parameter&                             V_param,
                        MAST::Parameter&                             kr_param,
                        MAST::Parameter&                             bref_param,
                        Real                                         rho,
                        Real                                         V_lower,
                        Real                                         V_upper,
                        unsigned int                                 n_V_divs,
                        Real                                         kr_lower,
                        Real                                         kr_upper,
                        unsigned int                                 n_kr_divs,
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
        
        
        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params
         */
        virtual void calculate_sensitivity(MAST::FlutterRootBase& root,
                                           const libMesh::ParameterVector& params,
                                           const unsigned int i);
        
        
    protected:
        
        
        void _insert_new_solution(const Real k_red_ref,
                                  MAST::FlutterSolutionBase* sol);
        
        
        virtual std::pair<bool, MAST::FlutterSolutionBase*>
        _bisection_search(const std::pair<MAST::FlutterSolutionBase*,
                          MAST::FlutterSolutionBase*>& ref_sol_range,
                          const unsigned int root_num,
                          const Real g_tol,
                          const unsigned int max_iters);
        
        
        
        /*!
         *    Newton method to look for cross-over point method search
         */
        /*virtual std::pair<bool, MAST::FlutterSolutionBase*>
        newton_search(const MAST::FlutterSolutionBase& init_sol,
                      const unsigned int root_num,
                      const Real tol,
                      const unsigned int max_iters);*/
        
        /*!
         *   performs an eigensolution at the specified reduced frequency, and
         *   sort the roots based on the provided solution pointer. If the
         *   pointer is nullptr, then no sorting is performed
         *   solve the eigenproblem  \f\[ L x = lambda R x \f\]
         */
        virtual std::unique_ptr<MAST::FlutterSolutionBase>
        _analyze(const Real k_red,
                 const Real v_ref,
                 const MAST::FlutterSolutionBase* prev_sol=nullptr);
        
        
        
        /*!
         *    initializes the matrices for the specified k_red. UG does not account
         *    for structural damping.
         */
        void _initialize_matrices(const Real k_red,
                                  const Real v_ref,
                                  ComplexMatrixX& L,   // stiff, aero, damp
                                  ComplexMatrixX& R,   // mass
                                  RealMatrixX& stiff); // stiffness
        
        
        /*!
         *    Assembles the reduced order system structural and aerodynmaic
         *    matrices for specified flight velocity \par U_inf.
         */
        void
        _initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                                 const unsigned int i,
                                                 const Real k_red,
                                                 const Real U_inf,
                                                 ComplexMatrixX& L,  // stiff, aero, damp
                                                 ComplexMatrixX& R); // mass

        
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
         *  Parameter that define the reduced frequency.
         */
        MAST::Parameter*                                _kred_param;
        
        /*!
         *    flight density
         */
        Real                                            _rho;

        /*!
         *    reference chord
         */
        MAST::Parameter*                                _bref_param;

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
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real>                           _kr_range;
        
        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int                                    _n_k_red_divs;
        
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


#endif // __mast__pk_flutter_solver_h__
