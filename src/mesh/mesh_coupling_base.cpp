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
#include "mesh/mesh_coupling_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"
#include "libmesh/point_locator_tree.h"



MAST::MeshCouplingBase::MeshCouplingBase(MAST::SystemInitialization& sys_init):
_sys_init(sys_init) {
    
    
}


MAST::MeshCouplingBase::~MeshCouplingBase() {
    
    
}




void
MAST::MeshCouplingBase::
add_master_and_slave_boundary_coupling(unsigned int master_b_id,
                                       unsigned int slave_b_id,
                                       Real tol) {
    
    // iterate on all nodes and check with boundary info about
    // current ids on it.
    libMesh::MeshBase& mesh = _sys_init.system().get_mesh();
    
    std::unique_ptr<libMesh::PointLocatorBase>
    pt_locator(mesh.sub_point_locator());
    
    libMesh::PointLocatorTree
    &locator_tree = dynamic_cast<libMesh::PointLocatorTree&>(*pt_locator);
    
    libMesh::MeshBase::const_element_iterator
    e_it   =  mesh.local_elements_begin(),
    e_end  =  mesh.local_elements_end();
    
    std::set<const libMesh::Node*> slave_nodes;
    
    for ( ; e_it != e_end; e_it++) {
        
        // iterate on sides and check if it is on specified boundary id
        for (unsigned int slave_side=0;
             slave_side < (*e_it)->n_sides();
             slave_side++) {
            
            // check if the side is on the specified boundary
            if (_check_if_side_on_boundary(mesh, **e_it, slave_side, slave_b_id)) {
                
                std::unique_ptr<const libMesh::Elem>
                slave_side_ptr((*e_it)->side_ptr(slave_side));
                
                // check for coupling of nodes on this side
                for (unsigned int slave_n_id=0;
                     slave_n_id < slave_side_ptr->n_nodes();
                     slave_n_id++) {
                    
                    const libMesh::Node*
                    slave_node = slave_side_ptr->node_ptr(slave_n_id);
                    
                    // if the node has not already been constrained, then
                    // identify the constraint on it
                    if (!slave_nodes.count(slave_node)) {
                        
                        std::set<const libMesh::Node*> master_nodes;
                        
                        std::set<const libMesh::Elem*>
                        elems = locator_tree.perform_fuzzy_linear_search(*slave_node, nullptr, tol);
                        
                        // make sure some elements were found for this node
                        libmesh_assert(elems.size());
                        
                        // now check which of these elements have sides on
                        // the boundary.
                        std::set<const libMesh::Elem*>::const_iterator
                        master_e_it   = elems.begin(),
                        master_e_end  = elems.end();
                        
                        for (; master_e_it!=master_e_end; master_e_it++) {
                            
                            
                            for (unsigned int master_side=0;
                                 master_side < (*master_e_it)->n_sides();
                                 master_side++) {
                                
                                if (_check_if_side_on_boundary(mesh,
                                                               **master_e_it,
                                                               master_side,
                                                               master_b_id)) {
                                    
                                    std::unique_ptr<const libMesh::Elem>
                                    master_side_ptr((*master_e_it)->side_ptr(master_side));
                                    
                                    // check for coupling of nodes on this side
                                    for (unsigned int master_n_id=0;
                                         master_n_id < master_side_ptr->n_nodes();
                                         master_n_id++)
                                        master_nodes.insert(master_side_ptr->node_ptr(master_n_id));
                                }
                            }
                        }
                            
                        // add this slave/master information to the data
                        _node_couplings.push_back
                        (std::pair<const libMesh::Node*, std::set<const libMesh::Node*>>
                         (slave_node, master_nodes));
                        
                        // register this slave node as having been processed
                        slave_nodes.insert(slave_node);
                    }
                }
            }
        }
    }
}


void
MAST::MeshCouplingBase::
add_slave_boundary_and_master_subdomain_coupling(unsigned int master_b_id,
                                                 unsigned int slave_b_id,
                                                 Real tol) {
    
    libmesh_assert(false); // to be implemented
}



bool
MAST::MeshCouplingBase::
_check_if_side_on_boundary(libMesh::MeshBase& mesh,
                           const libMesh::Elem& elem,
                           unsigned int side,
                           unsigned int b_id) {
    
    std::vector<libMesh::boundary_id_type> bc_ids;
    mesh.boundary_info->boundary_ids(&elem, side, bc_ids);
    
    bool
    on_side = false;
    for (unsigned int i=0; i<bc_ids.size(); i++)
        if (bc_ids[i] == b_id)
            on_side = true;
    
    return on_side;
}
