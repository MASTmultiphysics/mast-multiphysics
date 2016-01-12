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

#ifndef __mast__time_domain_flutter_solver_h__
#define __mast__time_domain_flutter_solver_h__

// C++ includes
#include <fstream>
#include <iomanip>


// MAST includes
#include "base/mast_data_types.h"


// libMesh includes
#include "libmesh/parameter_vector.h"
#include "libmesh/numeric_vector.h"


namespace MAST {
    
    // Forward declerations
    class Parameter;
    class FlutterModel;
    class TimeDomainFlutterRootBase;
    class FlutterSolutionBase;
    class FlutterRootCrossoverBase;
    class StructuralFluidInteractionAssembly;
    template <typename ValType> class BasisMatrix;
    
    
    /*!
     *   This implements a solver for a single parameter instability
     *   problem, for example a flutter solver where flight speed is the 
     *   primary parameter.
     */
    class TimeDomainFlutterSolver {
      
    public:
        
        /*!
         *   defalut constructor
         */
        TimeDomainFlutterSolver();
        
        
        virtual ~TimeDomainFlutterSolver();
        
        
        
        /*!
         *    initializes the data structres for a flutter solution.
         */
        void initialize(MAST::StructuralFluidInteractionAssembly&   assembly,
                        Real                                        V_lower,
                        Real                                        V_upper,
                        unsigned int                                n_V_divs,
                        std::vector<libMesh::NumericVector<Real>*>& basis);


        /*!
         *   clears the solutions stored from a previous analysis.
         */
        virtual void clear_solutions();
        
        /*!
         *    finds the number of critical points already identified in the
         *    procedure.
         */
        virtual unsigned int n_roots_found() const;
        
        
        /*!
         *   returns the \par n th root in terms of ascending velocity that is
         *   found by the solver
         */
        const MAST::TimeDomainFlutterRootBase& get_root(const unsigned int n) const;
        
        
        
        void set_output_file(std::string& nm) {
            
            _output.close();
            _output.open(nm.c_str(), std::ofstream::out);
        }
        
        
        /*!
         *   Looks through the list of flutter cross-over points and iteratively
         *   zooms in to find the cross-over point. This should be called only
         *   after scan_for_roots() has been called. Potential cross-over points
         *   are sorted with increasing velocity, and this method will attempt to
         *   identify the next critical root in the order.
         */
        virtual std::pair<bool, const MAST::TimeDomainFlutterRootBase*>
        find_next_root(const Real g_tol,
                       const unsigned int n_bisection_iters);
        
        
        
        /*!
         *   This method checks if the flutter root corresponding to the
         *   lowest velocity crossover has been calculated. If not, then it
         *   attempts to find that root using an iterative approach
         */
        virtual std::pair<bool, const MAST::TimeDomainFlutterRootBase*>
        find_critical_root(const Real g_tol,
                           const unsigned int n_bisection_iters);
        
        
        
        /*!
         *   Calculate the sensitivity of the flutter root with respect to the
         *   \par i^th parameter in params
         */
        virtual void calculate_sensitivity(MAST::TimeDomainFlutterRootBase& root,
                                           const libMesh::ParameterVector& params,
                                           const unsigned int i) {
            
            libmesh_error(); // to be implemented
        }
        
        
        /*!
         *   Prints the sorted roots to the \par output
         */
        void print_sorted_roots(std::ostream* output = NULL);
        
        
        /*!
         *   Prints the crossover points output. If no pointer to output is given
         *   then the output defined by set_output_file() is used.
         */
        void print_crossover_points(std::ostream* output = NULL);
        
        
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
        _analyze(const Real v_ref,
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
         *    matrices for specified flight velocity \par U_inf.
         *
         *    m  is the mass matrix,
         *    k  is the stiffness matrix,
         *    c  is the damping matrix,
         */
        void _initialize_matrices(Real U_inf,
                                  RealMatrixX& m,
                                  RealMatrixX& c,
                                  RealMatrixX& k);
        
        /*!
         *   identifies all cross-over and divergence points from analyzed
         *   roots
         */
        virtual void _identify_crossover_points();


        /*!
         *   structural assembly that provides the assembly of the system
         *   matrices.
         */
        MAST::StructuralFluidInteractionAssembly*       _assembly;
        
        
        /*!
         *   basis vector used to define the reduced order model
         */
        std::vector<libMesh::NumericVector<Real>*>*     _basis_vectors;
        
        
        /*!
         *   range of reference values within which to find flutter roots
         */
        std::pair<Real, Real> _V_range;


        /*!
         *    number of division in the reference value range for initial
         *    scanning
         */
        unsigned int _n_V_divs;


        
        /*!
         *    file to which the result will be written
         */
        std::ofstream _output;
        
        
        /*!
         *   map of velocity sorted flutter solutions
         */
        std::map<Real, MAST::FlutterSolutionBase*> _flutter_solutions;
        
        /*!
         *   the map of flutter crossover points versus average velocity of the
         *   two bounding roots
         */
        std::multimap<Real, MAST::FlutterRootCrossoverBase*> _flutter_crossovers;

    };
}


#endif // __mast__time_domain_flutter_solver_h__

