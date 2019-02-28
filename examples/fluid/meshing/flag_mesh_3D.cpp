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
#include "examples/fluid/meshing/flag_mesh_3D.h"



// libMesh includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel_mesh.h"


void
MAST::FlagMesh3D::init (const unsigned int panel_bc_id,
                        const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                        libMesh::UnstructuredMesh& mesh, libMesh::ElemType t) {
    
    libmesh_assert(divs.size() == 3);
    libmesh_assert(divs[0]->n_divs() == 3);
    libmesh_assert(divs[1]->n_divs() == 3);
    libmesh_assert(divs[2]->n_divs() == 3);
    libmesh_assert(panel_bc_id > 5);
    
    _flag_th   = divs[2]->div_location(2) - divs[2]->div_location(1);
    
    _x_in      = divs[0]->div_location(0);
    _x_out     = divs[0]->div_location(3);
    _x_le      = divs[0]->div_location(1);
    _x_te      = divs[0]->div_location(2);

    _y_dom_left    = divs[1]->div_location(0);
    _y_dom_right   = divs[1]->div_location(3);
    _y_flag_left   = divs[1]->div_location(1);
    _y_flag_right  = divs[1]->div_location(2);

    _z_lo      = divs[2]->div_location(0);
    _z_up      = divs[2]->div_location(3);
    _panel_bc_id = panel_bc_id;
    
    MeshInitializer::init(divs, mesh, t);
}



void
MAST::FlagMesh3D::process_mesh( ) {
    
    // check if the mesh is parallel
    const bool
    parallel_mesh = !_mesh->is_serial();

    {
        libMesh::MeshSerializer serializer(*_mesh);
        
        
        //march over all the elmeents and tag the sides that all lie on the panel suface
        libMesh::MeshBase::element_iterator e_it = _mesh->elements_begin();
        
        for ( ; e_it != _mesh->elements_end(); e_it++) {
            
            // if all nodes of the element are inside/on the bounadry of the
            // flag, then delete the element
            unsigned int
            n_nodes           = (*e_it)->n_nodes(),
            n_nodes_in_domain = 0;
            
            for (unsigned int i=0; i<n_nodes; i++) {
                
                const libMesh::Point&
                p = (*e_it)->point(i);
                
                if (p(0) >= _x_le            &&
                    p(0) <= _x_te            &&
                    p(1) >= _y_flag_left     &&
                    p(2) <= _y_dom_right     &&
                    p(2) >= -_flag_th        &&
                    p(2) <= _flag_th)
                    n_nodes_in_domain++;
            }
            
            if (n_nodes_in_domain == n_nodes)
                _mesh->delete_elem(*e_it);

        }
        _mesh->prepare_for_use();
    }
    
    {
        // now identify the boundary where the tag needs to be applied
        //march over all the elmeents and tag the sides that all lie on the panel suface

        libMesh::MeshSerializer serializer(*_mesh);
        
        libMesh::MeshBase::element_iterator e_it = _mesh->elements_begin();
        
        for ( ; e_it != _mesh->elements_end(); e_it++) {
            
            // iterate over the sides of each element and check if all
            // nodes satisfy the requirement
            
            for (unsigned int i_side=0; i_side<(*e_it)->n_sides(); i_side++) {
                
                std::unique_ptr<const libMesh::Elem> side_elem ((*e_it)->side_ptr(i_side).release());
                std::vector<bool> side_on_panel(side_elem->n_nodes());
                std::fill(side_on_panel.begin(), side_on_panel.end(), false);
                
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++) {
                    
                    const libMesh::Node& n = *(side_elem->node_ptr(i_node));
                    if ((n(2) >= -_flag_th) &&
                        (n(2) <= _flag_th) &&
                        (n(1) >= _y_flag_left) &&
                        (n(1) <= _y_flag_right) &&
                        (n(0) >= _x_le) &&
                        (n(0) <= _x_te))
                        side_on_panel[i_node] = true;
                }
                
                // check for side on panel
                bool if_apply_bc = true;
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                    if_apply_bc = side_on_panel[i_node] && if_apply_bc;
                if (if_apply_bc) {
                    _mesh->boundary_info->add_side(*e_it, i_side, _panel_bc_id);
                    if (parallel_mesh)
                        dynamic_cast<libMesh::DistributedMesh*>(_mesh)->add_extra_ghost_elem(*e_it);
                }
            }
        }
        
        // set the boudnary id names
        _mesh->boundary_info->sideset_name(_panel_bc_id) = "Panel";
    }
    
}



