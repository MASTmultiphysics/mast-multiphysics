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


#ifndef __mast__filter_base_h__
#define __mast__filter_base_h__


// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/system.h"



namespace MAST {
    
    /*!
     *   Creates a geometric filter for the level-set design variables.
     */
    class FilterBase {
    
    public:
        
        FilterBase(libMesh::System& sys,
                   const Real radius,
                   const std::set<unsigned int>& dv_dof_ids);
        
        
        virtual ~FilterBase();
        
        /*!
         *   computes the filtered output from the provided input.
         */
        template <typename ValType>
        void compute_filtered_values(const ValType& input,
                                     ValType& output) const;
    
    protected:
        
        /*!
         *   initializes the algebraic data structures
         */
        void _init();
        
        /*!
         *   system on which the level set discrete function is defined
         */
        libMesh::System& _level_set_system;
        
        /*!
         *   radius of the filter.
         */
        Real _radius;
        
        /*!
         *   dof ids that are design variables. If a id is not in this set, then
         *   the dof value assumes its value from the input
         */
        const std::set<unsigned int>&   _dv_dof_ids;

        /*!
         *   Algebraic relation between filtered level set values and the
         *   design variables \f$ \tilde{\phi}_i = B_{ij} \phi_j \f$
         */
        std::map<unsigned int, std::vector<std::pair<unsigned int, Real>>> _filter_map;
    };
    
    
}


#endif // __mast__filter_base_h__
