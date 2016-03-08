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
    class FlutterSolverBase;
    class FlutterRootBase;
    
    
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
         *    @returns the root
         */
        const MAST::FlutterRootBase&
        get_root(const unsigned int i) const {
            
            libmesh_assert_less(i, _roots.size());
            return *_roots[i];
        }
        
        
        /*!
         *    @returns a non-const reference to the root
         */
        MAST::FlutterRootBase&
        get_root(const unsigned int i) {
            
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
        virtual void sort(const MAST::FlutterSolutionBase& sol) = 0;
        
        
        /*!
         *
         */
        void swap_root(MAST::FlutterSolutionBase& sol,
                       unsigned int root_num);
                
        
        /*!
         *    prints the data and modes from this solution
         */
        virtual void print(std::ostream& output) = 0;
        
        
    protected:
        
        /*!
         *    Reference value of the sweeping parameter for which this solution
         *    was obtained. For UG solver, this is k_red, and for time domain
         *    solver this could be velocity. PK solver will need additional
         *    reference values, provided in the inherited class.
         */
        Real _ref_val;
        
        std::vector<MAST::FlutterRootBase*> _roots;
    };

}


#endif // __mast__flutter_solution_base_h__
