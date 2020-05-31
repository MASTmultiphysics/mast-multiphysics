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
#include "examples/fluid/meshing/panel_mesh_3D.h"

// libMesh includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/boundary_info.h"


void
MAST::PanelMesh3D::init (const Real tc,
                         bool cos_profile,
                         const unsigned int n_maxima_x,
                         const unsigned int n_maxima_y,
                         const unsigned int panel_bc_id,
                         const unsigned int symmetry_bc_id,
                         const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
                         libMesh::UnstructuredMesh& mesh, libMesh::ElemType t) {
    
    
    libmesh_assert_equal_to(divs.size(), 3);
    libmesh_assert_equal_to(divs[0]->n_divs(), 3);
    libmesh_assert_equal_to(divs[1]->n_divs(), 3);
    libmesh_assert_greater_equal(divs[2]->n_divs(), 1);
    libmesh_assert_greater(panel_bc_id, 5);
    libmesh_assert_greater(symmetry_bc_id, 5);
    
    _tc_ratio = tc;
    _x0 = divs[0]->div_location(1);
    _x1 = divs[0]->div_location(2);
    _y0 = divs[1]->div_location(1);
    _y1 = divs[1]->div_location(2);
    _z0 = divs[2]->div_location(0);
    _z1 = divs[2]->div_location(1);
    _n_maxima_x = n_maxima_x;
    _n_maxima_y = n_maxima_y;
    _cos_profile = cos_profile;
    _panel_bc_id = panel_bc_id;
    _symmetry_bc_id = symmetry_bc_id;
    MeshInitializer::init(divs, mesh, t);
}



void
MAST::PanelMesh3D::process_mesh( ) {
    
    
    // note that we apply the boundary conditions before moving the mesh since
    // the application of boudnary conditions is contingent upon the panel surface
    // points lying on the y = _y0, and the mesh movement ends up altering that
    
    // check if the mesh is parallel
    const bool
    parallel_mesh = !_mesh->is_serial();
    
    {
        libMesh::MeshSerializer serializer(*_mesh);
        
        //march over all the elmeents and tag the sides that all lie on the panel suface
        libMesh::MeshBase::element_iterator e_it = _mesh->elements_begin();
        
        for ( ; e_it != _mesh->elements_end(); e_it++)
        {
            // iterate over the sides of each element and check if all
            // nodes satisfy the requirement
            
            for (unsigned int i_side=0; i_side<(*e_it)->n_sides(); i_side++)
            {
                std::unique_ptr<const libMesh::Elem> side_elem ((*e_it)->side_ptr(i_side).release());
                std::vector<bool> side_on_panel(side_elem->n_nodes()),
                side_on_slip_wall(side_elem->n_nodes());
                
                std::fill(    side_on_panel.begin(),     side_on_panel.end(), false);
                std::fill(side_on_slip_wall.begin(), side_on_slip_wall.end(), false);
                
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++) {
                    
                    const libMesh::Node& n = *(side_elem->node_ptr(i_node));
                    if ((n(2)==_z0) && //bottom face
                        (n(0) >= _x0-1.0e-6) && (n(0) <= _x1+1.0e-6) && // x-coord
                        (n(1) >= _y0-1.0e-6) && (n(1) <= _y1+1.0e-6))
                        side_on_panel[i_node] = true;
                    
                    if ((n(2)==_z0) &&
                        ((n(0) <= _x0+1.0e-6) || (n(0) >= _x1-1.0e-6) || // x-coord
                         (n(1) <= _y0+1.0e-6) || (n(1) >= _y1-1.0e-6)))
                        side_on_slip_wall[i_node] = true;
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
                
                // now check for the slip wall
                if_apply_bc = true;
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                    if_apply_bc = side_on_slip_wall[i_node] && if_apply_bc;

                if (if_apply_bc)
                    _mesh->boundary_info->add_side(*e_it, i_side, _symmetry_bc_id);
            }
        }
        
        // set the boudnary id names
        _mesh->boundary_info->sideset_name(_panel_bc_id) = "Panel";
        _mesh->boundary_info->sideset_name(_symmetry_bc_id) = "Symmetry";
        
    }
    
    // now move the mesh points
    libMesh::MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
    const libMesh::MeshBase::node_iterator n_end = _mesh->nodes_end();
    
    const Real pi = acos(-1.);
    Real x, y, z;
    
    for (; n_it != n_end; n_it++)
    {
        libMesh::Node& n =  **n_it;
        
        // this is for sine bump
        if ((n(0) >= _x0) && (n(0) <= _x1))
        {
            x = n(0);
            y = n(1);
            z = n(2);
            
            if (!_cos_profile)
                n(2) += 0.5*_tc_ratio * ( 1. -(z - _z0)/(_z1-_z0) ) *
                sin( _n_maxima_x * pi*(x-_x0)/(_x1-_x0)) *
                sin( _n_maxima_y * pi*(y-_y0)/(_y1-_y0));
            else
                n(2) += 0.5*_tc_ratio * ( 1. -(z-_z0)/(_z1-_z0) ) *
                (1. - cos( 2*_n_maxima_x * pi*(x-_x0)/(_x1-_x0))) *
                (1. - cos( 2*_n_maxima_y * pi*(y-_y0)/(_y1-_y0)));
        }
    }
}

