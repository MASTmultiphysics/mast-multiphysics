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

#ifndef __mast__mesh_couplings__
#define __mast__mesh_couplings__

// C++ includes
#include <set>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/node.h"


namespace MAST {
  
    // Forward declerations
    class SystemInitialization;

    
    class MeshCouplingBase {
 
    public:
    
        MeshCouplingBase(MAST::SystemInitialization& sys_init);
        
        
        virtual ~MeshCouplingBase();
        
        
        void
        add_master_and_slave_boundary_coupling(unsigned int master_b_id,
                                               unsigned int slave_b_id,
                                               Real tol);

        
        void
        add_slave_boundary_and_master_subdomain_coupling(unsigned int master_id,
                                                         unsigned int slave_b_id,
                                                         Real tol);
        
        /*!
         *  \returns a const reference to a vector of pairs where each pair
         *  identifies the slave and the set of master nodes to which the
         *  slave is connected.
         */
        const std::vector<std::pair<const libMesh::Node*, std::set<const libMesh::Node*>>>&
        get_node_couplings() const {
            
            return _node_couplings;
        }
        
        
    protected:

        bool
        _check_if_side_on_boundary(libMesh::MeshBase& mesh,
                                   const libMesh::Elem& elem,
                                   unsigned int side,
                                   unsigned int b_id);
        
        
        MAST::SystemInitialization &_sys_init;
        
        std::vector<std::pair<const libMesh::Node*, std::set<const libMesh::Node*>>> _node_couplings;
        
    };
    
}


#endif // __mast__mesh_couplings__
