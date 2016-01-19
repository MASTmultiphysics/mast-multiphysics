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

#ifndef __mast__flutter_solution_base_h__
#define __mast__flutter_solution_base_h__

// C++ includes
#include <vector>


// MAST includes
#include "base/mast_data_types.h"


namespace MAST {
    
    // Forward declerations
    class TimeDomainFlutterSolver;
    class TimeDomainFlutterRootBase;
    class LAPACK_DGGEV;
    
    
    class FlutterSolutionBase {
        
    public:
        FlutterSolutionBase():
        _ref_val(0.)
        {}
        
        /*!
         *    delete the flutter root objects
         */
        virtual ~FlutterSolutionBase();


        /*!
         *   initializes the root
         */
        void init (const MAST::TimeDomainFlutterSolver& solver,
                   const Real v_ref,
                   const MAST::LAPACK_DGGEV& eig_sol);

        
        /*!
         *   the reduced frequency for this solution
         */
        Real ref_val() const {
            return _ref_val;
        }
        
        /*!
         *   number of roots in this solution
         */
        unsigned int n_roots() const{
            return (unsigned int)_roots.size();
        }

        
        /*!
         *   number of unstable roots in this solution. Only roots with damping 
         *   greater than \par tol will be considered unstable.
         */
        unsigned int n_unstable_roots_in_upper_complex_half (Real tol) const;

        
        /*!
         *    @returns the critical root at the lowest velocity
         */
        MAST::TimeDomainFlutterRootBase* get_critical_root();

        
        /*!
         *    @returns the root
         */
        const MAST::TimeDomainFlutterRootBase&
        get_root(const unsigned int i) const {
            
            libmesh_assert_less(i, _roots.size());
            return *_roots[i];
        }
        
        /*!
         *    @returns a non-const reference to the root
         */
        MAST::TimeDomainFlutterRootBase& get_root(const unsigned int i) {
            libmesh_assert_less(i, _roots.size());
            return *_roots[i];
        }
        
        /*!
         *    sort this root with respect to the given solution from a previous
         *    eigen solution. This method relies on the modal participation.
         *    Flutter roots from previous and current solutions with highest
         *    dot product of modal participation vector are considered to be
         *    similar.
         */
        void sort(const MAST::FlutterSolutionBase& sol);
        
        /*!
         *
         */
        void swap_root(MAST::FlutterSolutionBase& sol,
                       unsigned int root_num);
        
        /*!
         *    prints the data and modes from this solution
         */
        void print(std::ostream& output);
        
    protected:
        
        /*!
         *    Reference value of the sweeping parameter for which this solution
         *    was obtained. For UG solver, this is k_red, and for time domain
         *    solver this could be velocity. PK solver will need additional
         *    reference values, provided in the inherited class.
         */
        Real _ref_val;
        
        /*!
         *    Matrix used for scaling of eigenvectors, and sorting of roots
         */
        RealMatrixX _Amat, _Bmat;
        
        std::vector<MAST::TimeDomainFlutterRootBase*> _roots;
    };

}


#endif // __mast__flutter_solution_base_h__
