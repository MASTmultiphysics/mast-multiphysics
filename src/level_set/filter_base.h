/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "libmesh/node.h"
#include "libmesh/elem.h"


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
        void compute_filtered_values(const libMesh::NumericVector<Real>& input,
                                     libMesh::NumericVector<Real>& output,
                                     bool close_vector = true) const;

        /*!
         *  for large problems it is more efficient to specify only the non-zero entries in the input vector in
         *  \p nonzero_vals. If \p output is of type SERIAL then all ranks in the
         *  communicator will perform the same operaitons and provide an identical \p output vector.
         *  If \p output is of type PARALLEL or GHOSTED then only the local contributions will be set.
         *  If \p close_vector is \p true then \p output.close() will be called in this
         *  routines, otherwise not.
         */
        void compute_filtered_values(std::map<unsigned int, Real>& nonzero_vals,
                                     libMesh::NumericVector<Real>& output,
                                     bool close_vector) const;
        
        void compute_filtered_values(const std::vector<Real>& input,
                                     std::vector<Real>& output) const;

        /*!
         *   function identifies if the given element is within the domain of
         *   influence of this specified level set design variable. Currently,
         *   this is identified based on the filter radius, the distance of
         *   element nodes from the specified level set design variable location
         *   and the element sizes.
         */
        bool if_elem_in_domain_of_influence(const libMesh::Elem& elem,
                                            const libMesh::Node& level_set_node) const;
        

        
        /*!
         *  prints the filter data.
         */
        virtual void print(std::ostream& o) const;
        
        
    protected:
        
        /*!
         *   initializes the algebraic data structures
         */
        void _init();
        void _init2();
        
        /*!
         *   system on which the level set discrete function is defined
         */
        libMesh::System& _level_set_system;
        
        /*!
         *   radius of the filter.
         */
        Real _radius;
        
        
        /*!
         *   largest element size in the level set mesh
         */
        Real _level_set_fe_size;
        
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
