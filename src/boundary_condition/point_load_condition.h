/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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

#ifndef __mast__point_load_condition__
#define __mast__point_load_condition__

// C++ includes
#include <set>

// MAST includes
#include "base/boundary_condition_base.h"

// libMesh includes
#include "libmesh/node.h"

namespace MAST {
    
    /*!
     *   This class allows for the specification of load associated with 
     *   specified nodes in a user-provided set. The user is responsible for
     *   maintaining consistency of the nodes during mesh-refinement.
     */
    class PointLoadCondition:
    public MAST::BoundaryConditionBase {
        
    public:
        
        PointLoadCondition(MAST::BoundaryConditionType t);
        
        virtual ~PointLoadCondition();
        
        
        /*!
         *   @returns the set of nodes on which the load is 
         *   specified as a constant reference
         */
        const std::set<libMesh::Node*>& get_nodes() const;
        

        /*!
         *   @returns the set of nodes on which the load is
         *   specified as a writable reference
         */
        std::set<libMesh::Node*>& get_nodes();

        
    protected:
        
        
        /*!
         *   set of nodes on which load is specified
         */
        std::set<libMesh::Node*> _nodes;
        
    };
}

#endif // __mast__point_load_condition__
