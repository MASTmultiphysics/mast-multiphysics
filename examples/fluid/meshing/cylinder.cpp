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


// MAST includes
#include "examples/fluid/meshing/cylinder.h"
#include "examples/fluid/meshing/mesh_initializer.h"

// libMesh includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/boundary_info.h"


MAST::Examples::CylinderMesh2D::CylinderMesh2D() {
    
    
}


MAST::Examples::CylinderMesh2D::~CylinderMesh2D() {
    
    
}


void
MAST::Examples::CylinderMesh2D::mesh(const Real r,
                                     const Real L,
                                     const unsigned int radial_divs,
                                     const unsigned int quarter_divs,
                                     const Real elem_size_ratio,
                                     libMesh::UnstructuredMesh& mesh,
                                     libMesh::ElemType etype,
                                     bool  add_downstream_block,
                                     const Real downstream_block_length,
                                     const unsigned int downstream_block_divs) {

    
    Real
    eta   = 0.,
    phi   = 0.,
    pi    = acos(-1.);

    libMesh::MeshTools::Generation::build_square(mesh,
                                                 radial_divs, 4 * quarter_divs,
                                                 0., 1.,
                                                 -1., 7.,
                                                 etype);

    std::unique_ptr<libMesh::ReplicatedMesh>
    block_mesh;
    
    if (add_downstream_block) {
        block_mesh.reset(new libMesh::ReplicatedMesh(mesh.comm()));
        libMesh::MeshTools::Generation::build_square(*block_mesh,
                                                     downstream_block_divs, quarter_divs,
                                                     0., 1.,
                                                     -1., 1.,
                                                     etype);
    }
    
    std::vector<Real>
    div_loc      = {0., 1.},
    dx_relative  = {1., elem_size_ratio};
    std::vector<unsigned int>
    n_dx         = {20};
    
    MAST::MeshInitializer::CoordinateDivisions divs;
    divs.init(1, div_loc, dx_relative, n_dx);
    
    // collect nodes at the beginning and end of the domain so that
    // we can close the domain by removing nodes at phi = 7 and use the
    // nodes at phi = -1 instead.
    
    libMesh::MeshSerializer serializer(mesh);

    libMesh::MeshBase::node_iterator
    n_it    = mesh.nodes_begin(),
    n_end   = mesh.nodes_end();
    
    
    Real
    tol = 1.e-8*r;

    //create a map of old to new node
    std::map<libMesh::Node*, libMesh::Node*>
    node_map;
    
    std::map<Real, libMesh::Node*>
    nodes_to_use,
    nodes_to_delete,
    nodes_for_block;
    
    {
        for ( ; n_it != n_end; n_it++) {
            
            libMesh::Node& n = **n_it;
            
            if (fabs(n(1) + 1.) < tol)
                nodes_to_use[n(0)] = &n;
            else if (fabs(n(1) - 7.) < tol)
                nodes_to_delete[n(0)] = &n;
        }
        
        libmesh_assert_equal_to(nodes_to_use.size(), nodes_to_delete.size());

        std::map<Real, libMesh::Node*>::iterator
        n1_it  = nodes_to_use.begin(),
        n2_it  = nodes_to_delete.begin();
        
        for ( ; n2_it != nodes_to_delete.end(); n2_it++) {
            
            node_map[n2_it->second] = n1_it->second;
            n1_it++;
        }
    }

    
    n_it    = mesh.nodes_begin();
    n_end   = mesh.nodes_end();

    for ( ; n_it != n_end; n_it++) {

        libMesh::Node& n = **n_it;
        
        eta  = divs(n(0));
        phi  = n(1);
        
        if (phi <= 1.) {
            
            // identify the node to be used for the block
            if (n(0)-1. < tol)
                nodes_for_block[n(1)] = &n;
            
            n(0) =   (1. - eta) * r * cos(phi * 0.25 * pi) + eta * L;
            n(1) =  ((1. - eta) * r * sin(phi * 0.25 * pi) + eta * L * phi);
        }
        else if (phi < 3.) {
            phi -= 2.;
            n(0) = -((1. - eta) * r * sin(phi * 0.25 * pi) + eta * L * phi);
            n(1) =   (1. - eta) * r * cos(phi * 0.25 * pi) + eta * L;
        }
        else if (phi < 5.) {
            phi -= 4.;
            n(0) = -((1. - eta) * r * cos(phi * 0.25 * pi) + eta * L);
            n(1) = -((1. - eta) * r * sin(phi * 0.25 * pi) + eta * L * phi);
        }
        else if (phi <= 7.) {
            phi -= 6.;
            n(0) =  ((1. - eta) * r * sin(phi * 0.25 * pi) + eta * L * phi);
            n(1) = -((1. - eta) * r * cos(phi * 0.25 * pi) + eta * L);
        }
    }

    
    // iterate over the elements and reassign the boundary nodes, to close
    // the mesh
    libMesh::MeshBase::element_iterator
    e_it  =  mesh.elements_begin(),
    e_end =  mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        
        libMesh::Elem& e = **e_it;
        
        for (unsigned int i=0; i<e.n_nodes(); i++)
            if (node_map.count(e.node_ptr(i)))
                e.set_node(i) = node_map[e.node_ptr(i)];
    }
    
    
    // delete the nodes
    std::map<Real, libMesh::Node*>::iterator
    nodes_it   = nodes_to_delete.begin();
    
    for (; nodes_it != nodes_to_delete.end(); nodes_it++)
        mesh.delete_node(nodes_it->second);

    // if the downstream block is to be added, do so now. This is done by
    // iterating over all ndoes and elements on the block mesh and replicating
    // those in the new mesh
    if (add_downstream_block) {
        
        std::map<const libMesh::Node*, libMesh::Node*>
        block_node_map;
        
        nodes_to_delete.clear();
        
        libMesh::MeshBase::node_iterator
        n_it   = block_mesh->nodes_begin(),
        n_end  = block_mesh->nodes_end();
        
        libMesh::Point p;
        
        for ( ; n_it != n_end; n_it++) {
            
            libMesh::Node* old_node = *n_it;
            
            if ((*old_node)(0) < tol) {
                // If node is on left edge of block, we will use the node from
                // the circular mesh already processed.
                nodes_to_delete[(*old_node)(1)] = old_node;
            }
            else {
                // we will create a new node for this element
                p(0) = L + (*old_node)(0) * downstream_block_length;
                p(1) = L * (*old_node)(1);
                block_node_map[*n_it] = mesh.add_point(p);
            }
        }
        
        // we will not identify which nodes to use on the left edge
        libmesh_assert_equal_to(nodes_to_delete.size(), nodes_for_block.size());

        std::map<Real, libMesh::Node*>::iterator
        n1_it  = nodes_for_block.begin(),
        n2_it  = nodes_to_delete.begin();
        
        for ( ; n2_it != nodes_to_delete.end(); n2_it++) {
            
            block_node_map[n2_it->second] = n1_it->second;
            n1_it++;
        }
        
        // now, we should have identified all nodes for the block in this map
        libmesh_assert_equal_to(block_node_map.size(), block_mesh->n_nodes());
        
        // from the old mesh, remove sides that belong to the right edge
        libMesh::MeshBase::element_iterator
        e_it    = mesh.elements_begin(),
        e_end   = mesh.elements_end();
        
        for ( ; e_it != e_end; e_it++) {
            
            libMesh::Elem* elem = *e_it;
            
            std::vector<libMesh::boundary_id_type> bids;
            
            // check if the right side has a boudnary id
            mesh.boundary_info->boundary_ids(elem, 1, bids);
            std::unique_ptr<libMesh::Elem> edge(elem->side_ptr(1));
            if (bids.size()==1 &&
                bids[0] == 1 &&
                std::fabs(edge->centroid()(0)-L) < tol) // on the right edge
                mesh.boundary_info->remove_side(elem, 1, 1);
        }

        // now add elements to the mesh and set their nodes
        e_it   = block_mesh->elements_begin();
        e_end  = block_mesh->elements_end();
        
        for ( ; e_it != e_end; e_it++) {
            
            libMesh::Elem* old_elem = *e_it;
            
            libMesh::Elem*
            new_elem = mesh.add_elem(libMesh::Elem::build(old_elem->type()).release());
            
            bool
            on_side = false;
            std::set<unsigned int> sides;
            
            for (unsigned int i=0; i<old_elem->n_nodes(); i++) {
                
                libMesh::Node& old_node = *old_elem->node_ptr(i);
                
                libmesh_assert(block_node_map.count(&old_node));
                new_elem->set_node(i) = block_node_map[&old_node];
                
                // check to see if the element side has to be added to the side
                // set
                if (std::fabs(old_node(0)-1.) < tol) {
                    // right edge
                    on_side = true;
                    sides.insert(1);
                }
                if (std::fabs(old_node(1)+1) < tol) {
                    // bottom
                    on_side = true;
                    sides.insert(0);
                }
                if (std::fabs(old_node(1)-1.) < tol) {
                    // top
                    on_side = true;
                    sides.insert(2);
                }
            }
            
            if (on_side) {
                std::set<unsigned int>::const_iterator
                s_it  = sides.begin(),
                s_end = sides.end();
                
                // add sides to far-field
                for ( ; s_it != s_end; s_it++)
                    mesh.boundary_info->add_side(new_elem, *s_it, 1);
            }
        }
    }
    
    mesh.prepare_for_use();
    mesh.boundary_info->sideset_name(1) = "farfield";
    mesh.boundary_info->sideset_name(3) = "cylinder";
    mesh.boundary_info->remove_id(0);
    mesh.boundary_info->remove_id(2);
}



