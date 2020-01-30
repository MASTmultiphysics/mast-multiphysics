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

#ifndef __mast_sub_elem_node_map_h__
#define __mast_sub_elem_node_map_h__

// C++ includes
#include <unordered_map>
#include <functional> // std::hash

// libMesh includes
#include "libmesh/mesh_base.h"

namespace MAST{

    
    class SubElemNodeMap {
        
        // borrowing this from libMesh::TopologyMap
        struct myhash {
        public:
            template <typename T1, typename T2>
            std::size_t operator()(const std::pair<T1, T2> & x) const
            {
                return 3 * std::hash<T1>()(x.first) + std::hash<T2>()(x.second);
            }
        };
        

        // first pair is the key, which identifies the bounding nodes on a
        // parent edge where the node is created.
        // second pair is the set of nodes used by the adjacent elements on the
        // positive and negative sides of the level set function. Specifying
        // both to be same will allow a weak discontinuity, while specifying
        // both to be different will allow a strong discontinuity
        typedef
        std::unordered_map
        <std::pair<libMesh::dof_id_type, libMesh::dof_id_type>,
        std::pair<libMesh::Node*, libMesh::Node*>, MAST::SubElemNodeMap::myhash> map_type;
        
        
    public:

        SubElemNodeMap();
        
        
        virtual ~SubElemNodeMap();
        
        
        void clear() { _map.clear(); }
        
        
        bool empty() const { return _map.empty(); }
        
        
        unsigned int
        count(libMesh::dof_id_type bracket_node1, libMesh::dof_id_type bracket_node2) const;


        std::pair<libMesh::Node*, libMesh::Node*>&
        add(libMesh::dof_id_type bracket_node1, libMesh::dof_id_type bracket_node2);
        
        
    protected:
        
        
        MAST::SubElemNodeMap::map_type          _map;
    };
    
}


#endif // __mast_sub_elem_node_map_h__
