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

#ifndef __mast__time_domain_flutter_solution_h__
#define __mast__time_domain_flutter_solution_h__

// C++ includes
#include <vector>


// MAST includes
#include "aeroelasticity/flutter_solution_base.h"


namespace MAST {
    
    // Forward declerations
    class TimeDomainFlutterSolver;
    class LAPACK_DGGEV;

    
    class TimeDomainFlutterSolution:
    public MAST::FlutterSolutionBase {
        
    public:
        TimeDomainFlutterSolution();
        
        /*!
         *    delete the flutter root objects
         */
        virtual ~TimeDomainFlutterSolution();
        
        
        /*!
         *   initializes the root
         */
        void init (const MAST::TimeDomainFlutterSolver& solver,
                   const Real v_ref,
                   const MAST::LAPACK_DGGEV& eig_sol);
        
        /*!
         *   number of unstable roots in this solution. Only roots with damping
         *   greater than \p tol will be considered unstable.
         */
        unsigned int n_unstable_roots_in_upper_complex_half (Real tol) const;
        
        /*!
         *    sort this root with respect to the given solution from a previous
         *    eigen solution. This method relies on the modal participation.
         *    Flutter roots from previous and current solutions with highest
         *    dot product of modal participation vector are considered to be
         *    similar.
         */
        virtual void sort(const MAST::FlutterSolutionBase& sol);

        
        /*!
         *    prints the data and modes from this solution
         */
        virtual void print(std::ostream& output);

        
        /*!
         *    @returns the critical root at the lowest velocity
         */
        MAST::FlutterRootBase* get_critical_root(Real tol);

        
    protected:
        
        /*!
         *    Matrix used for scaling of eigenvectors, and sorting of roots
         */
        RealMatrixX _Amat, _Bmat;
        
    };
    
}


#endif // __mast__time_domain_flutter_solution_h__

