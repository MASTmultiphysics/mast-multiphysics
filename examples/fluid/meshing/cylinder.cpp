/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"


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
                                     libMesh::ElemType etype) {

    
    Real
    eta   = 0.,
    phi   = 0.,
    pi    = acos(-1.);

    libMesh::MeshTools::Generation::build_square(mesh,
                                                 radial_divs, 4 * quarter_divs,
                                                 0., 1.,
                                                 -1., 7.,
                                                 etype);

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
    
    
    libMesh::MeshBase::node_iterator
    n_it    = mesh.nodes_begin(),
    n_end   = mesh.nodes_end();
    
    
    Real
    tol = 1.e-8;

    //create a map of old to new node
    std::map<libMesh::Node*, libMesh::Node*>
    node_map;
    
    std::map<Real, libMesh::Node*>
    nodes_to_use,
    nodes_to_delete;
    
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

    mesh.prepare_for_use();
    mesh.boundary_info->sideset_name(1) = "farfield";
    mesh.boundary_info->sideset_name(3) = "cylinder";
    mesh.boundary_info->remove_id(0);
    mesh.boundary_info->remove_id(2);
    mesh.write("mesh.exo");
    eta = 0.;
}



