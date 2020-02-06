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
#include "examples/fluid/meshing/flag_mesh_2D.h"



// libMesh includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"


MAST::FlagMesh2D::FlagMesh2D(const unsigned int n_flags,
                             const unsigned int panel_bc_id,
                             const std::vector<MeshInitializer::CoordinateDivisions*>& divs):
_x_in(0.), _x_out(0.),
_x_le(0.),  _x_te(0.),
_y_lo(0.),  _y_up(0.),
_n_flags(n_flags),
_panel_bc_id (panel_bc_id),
_divs(divs) {
    
    libmesh_assert_greater(n_flags,                      0);
    libmesh_assert_equal_to(divs.size(),                 2);
    libmesh_assert_equal_to(divs[0]->n_divs(),           3);
    libmesh_assert_equal_to(divs[1]->n_divs(), 1+2*n_flags);
    libmesh_assert_greater(panel_bc_id,                  3);
    
    _x_in      = _divs[0]->div_location(0);
    _x_out     = _divs[0]->div_location(3);
    _x_le      = _divs[0]->div_location(1);
    _x_te      = _divs[0]->div_location(2);
    _y_lo      = _divs[1]->div_location(0);
    _y_up      = _divs[0]->div_location(3);
    
    // store the flag details
    _flag_th.resize(_n_flags);
    _flag_center.resize(_n_flags);
    
    for (unsigned int i=0; i<n_flags; i++) {
        
        _flag_th[i]     =        divs[1]->div_location(2*(i+1)) - divs[1]->div_location(2*(i+1)-1);
        _flag_center[i] = 0.5 * (divs[1]->div_location(2*(i+1)) + divs[1]->div_location(2*(i+1)-1));
    }
}



void
MAST::FlagMesh2D::init_fluid_mesh (libMesh::UnstructuredMesh& mesh,
                                   libMesh::ElemType t) {

    MeshInitializer::init(_divs, mesh, t);
}



void
MAST::FlagMesh2D::process_mesh( ) {
    
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
                
                if (p(0) >= _x_le     &&
                    p(0) <= _x_te) {
                    
                    // now identify if this is in any one of the flag domains
                    bool
                    if_inside = false;
                    
                    unsigned int
                    i = 0;
                    
                    while ( i<_n_flags  && !if_inside) {

                        if (p(1) >= _flag_center[i] - 0.5*_flag_th[i]   &&
                            p(1) <= _flag_center[i] + 0.5*_flag_th[i])
                            if_inside = true;
                        
                        i++;
                    }
                    
                    if (if_inside)
                        n_nodes_in_domain++;
                }
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
                    if ((n(0) >= _x_le) &&
                        (n(0) <= _x_te)) {
                        
                        // now identify if this is in any one of the flag domains
                        unsigned int
                        i = 0;
                        
                        while ( i<_n_flags  && !side_on_panel[i_node]) {
                            
                            if (n(1) >= _flag_center[i] - 0.5*_flag_th[i]   &&
                                n(1) <= _flag_center[i] + 0.5*_flag_th[i])
                                side_on_panel[i_node] = true;
                            
                            i++;
                        }
                    }
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




void
MAST::FlagMesh2D::init_structural_mesh (libMesh::UnstructuredMesh& mesh,
                                        libMesh::ElemType t) {
    
    std::vector<Real>
    x_div_loc        (2),
    x_relative_dx    (2);
    
    std::vector<unsigned int>
    x_divs           (1);
    
    std::unique_ptr<MAST::MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MAST::MeshInitializer::CoordinateDivisions);
    
    std::vector<MAST::MeshInitializer::CoordinateDivisions*>
    divs(1);

    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<2; i_div++) {
        
        x_div_loc[i_div]        = _divs[0]->div_location(i_div+1);
        x_relative_dx[i_div]    = _divs[0]->relative_mesh_size(i_div+1);
    }
    x_divs[0]       = _divs[0]->n_elements_in_div(0);
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(1, x_div_loc, x_relative_dx, x_divs);
    
    // initialize a local mesh and copy them for each flag
    for (unsigned int i=0; i<_n_flags; i++) {

        libMesh::SerialMesh
        local_mesh(mesh.comm());
    
        MeshInitializer().init(divs, local_mesh, t);

        std::map<libMesh::Node*, libMesh::Node*>
        added_nodes;
        
        // displace the y-coordinate for all nodes to the mid-plane for this
        // flag
        libMesh::MeshBase::node_iterator
        n_it   = local_mesh.nodes_begin(),
        n_end  = local_mesh.nodes_end();
        
        for ( ; n_it != n_end; n_it++) {
            
            (**n_it)(1) = _flag_center[i];
            added_nodes[*n_it] = mesh.add_point((**n_it));
        }
        
        // now copy all the elements
        libMesh::MeshBase::element_iterator
        e_it   = local_mesh.elements_begin(),
        e_end  = local_mesh.elements_end();
        
        libMesh::Elem
        *e_new = nullptr;
        
        unsigned int
        elem_num = 0;
        
        for ( ; e_it != e_end; e_it++) {
            
            e_new  = copy_elem(**e_it, added_nodes);
            e_new  = mesh.add_elem(e_new);
            
            // elem num
            if (elem_num == 0)
                mesh.boundary_info->add_side(e_new, 0, 0);
            
            if (elem_num == local_mesh.n_elem()-1)
                mesh.boundary_info->add_side(e_new, 1, 1);
            
            // this is used to identify the element for which boundary
            // info must be set
            elem_num++;
        }
    }
    
    mesh.prepare_for_use();
    
    // Add sideset names to boundary info
    mesh.boundary_info->sideset_name(0) = "left";
    mesh.boundary_info->sideset_name(1) = "right";
}


libMesh::Elem*
MAST::FlagMesh2D::copy_elem(libMesh::Elem& e,
                            std::map<libMesh::Node*, libMesh::Node*>& added_nodes){

    libMesh::Elem*
    e_new = nullptr;
    
    switch (e.type()) {
            
        case libMesh::EDGE2:
            e_new = new libMesh::Edge2;
            break;
            
        case libMesh::EDGE3:
            e_new = new libMesh::Edge3;
            break;

        case libMesh::EDGE4:
            e_new = new libMesh::Edge4;
            break;

        default:
            libmesh_error();// invalid elem type
    }
    
    for (unsigned int i=0; i<e.n_nodes(); i++) {
        
        std::map<libMesh::Node*, libMesh::Node*>::iterator
        n_it = added_nodes.find(e.node_ptr(i));
        
        // make sure that this node exists in the map
        libmesh_assert(n_it != added_nodes.end());
        
        e_new->set_node(i) = n_it->second;
    }
    
    return e_new;
}


