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
#include "examples/fluid/meshing/mesh_initializer.h"



void
MAST::MeshInitializer::init (const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                             libMesh::UnstructuredMesh& mesh,
                             libMesh::ElemType t) {
    
    unsigned int dim = divs.size();
    _mesh = &mesh;
    
    // first create the mesh
    switch (dim) {
        case 1: {
            libMesh::MeshTools::Generation::build_line (mesh,
                                                        divs[0]->total_elem_divs(),
                                                        0., 1., t);
            
        }
            break;
            
        case 2: {
            libMesh::MeshTools::Generation::build_square (mesh,
                                                          divs[0]->total_elem_divs(),
                                                          divs[1]->total_elem_divs(),
                                                          0., 1., 0., 1., t);
        }
            break;
            
        case 3: {
            libMesh::MeshTools::Generation::build_cube (mesh,
                                                        divs[0]->total_elem_divs(),
                                                        divs[1]->total_elem_divs(),
                                                        divs[2]->total_elem_divs(),
                                                        0., 1., 0., 1., 0., 1., t);
        }
            break;
            
        default:
            libmesh_error();
            break;
    }
    
    // now iterate over the nodes in this mesh and update its coordinate
    libMesh::MeshBase::node_iterator n_it = mesh.nodes_begin();
    for (; n_it != mesh.nodes_end(); n_it++) {
        libMesh::Node& n = **n_it;
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            n(i_dim) = (*divs[i_dim])(n(i_dim));
    }
    
    // move the grid points and apply the boundary conditions
    this->process_mesh();
}



