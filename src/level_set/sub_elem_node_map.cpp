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

// MAST includes
#include "level_set/sub_elem_node_map.h"


MAST::SubElemNodeMap::SubElemNodeMap() {
    
}
    
    
MAST::SubElemNodeMap::~SubElemNodeMap() {
    
}
    
    


unsigned int
MAST::SubElemNodeMap::count(libMesh::dof_id_type bracket_node1,
                            libMesh::dof_id_type bracket_node2) const {
    
    return _map.count(std::make_pair(bracket_node1, bracket_node2));
}
    
    
std::pair<libMesh::Node*, libMesh::Node*>&
MAST::SubElemNodeMap::add(libMesh::dof_id_type bracket_node1,
                          libMesh::dof_id_type bracket_node2)  {
    
    MAST::SubElemNodeMap::map_type::iterator
    it  =  _map.find(std::make_pair(bracket_node1, bracket_node2));
    
    if (it == _map.end()) {
        
        it = _map.insert(std::make_pair(std::make_pair(bracket_node1, bracket_node2),
                                        std::make_pair(nullptr, nullptr))).first;
    }
    
    return it->second;
}

