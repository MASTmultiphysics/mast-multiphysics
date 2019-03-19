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
#include "examples/structural/base/blade_stiffened_panel_mesh.h"
#include "examples/base/multilinear_interpolation.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"


// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"



void
MAST::StiffenedPanelMesh::init(const unsigned int n_stiff,
                               const unsigned int n_x_divs,
                               const unsigned int n_y_divs_between_stiffeners,
                               const Real length,
                               const Real width,
                               libMesh::MeshBase& mesh,
                               libMesh::ElemType t,
                               bool beam_stiffeners,
                               const unsigned int n_z_divs_stiff,
                               const Real h_stiff) {
    
    // shell modeling of stiffeners should be accompanied by additional data
    // modeling the stiffeners
    if (!beam_stiffeners) {
        
        libmesh_assert(n_z_divs_stiff);
        libmesh_assert(h_stiff);
    }
    
    
    mesh.set_mesh_dimension(2);
    
    
    libMesh::ElemType stiff_t = t;
    
    if (beam_stiffeners) {
        
        switch (t) {
            case libMesh::TRI3:
            case libMesh::QUAD4:
                stiff_t = libMesh::EDGE2;
                break;
                
            case libMesh::TRI6:
            case libMesh::QUAD8:
            case libMesh::QUAD9:
                stiff_t = libMesh::EDGE3;
                break;
                
            default:
                libmesh_error();
        }
    }
    
    // initialize the main mesh
    libMesh::SerialMesh panel(mesh.comm());
    
    // create the mesh for each panel
    libMesh::MeshTools::Generation::build_square(panel,
                                                 n_x_divs,
                                                 n_y_divs_between_stiffeners*(n_stiff+1),
                                                 0., length,
                                                 0., width,
                                                 t);
    
    // use subdomain id for panel as 0
    _combine_mesh(mesh, panel, MAST::StiffenedPanelMesh::PANEL, 0., 0);
    
    
    // now iterate over the stiffeners and create them
    for (unsigned int i=0; i<n_stiff; i++) {
        
        libMesh::SerialMesh stiff(mesh.comm());
        
        if (beam_stiffeners)
            libMesh::MeshTools::Generation::build_line(stiff,
                                                       n_x_divs,
                                                       0.,
                                                       length,
                                                       stiff_t);
        else
            libMesh::MeshTools::Generation::build_square(stiff,
                                                         n_x_divs,
                                                         n_z_divs_stiff,
                                                         0., length,
                                                         0., h_stiff,
                                                         stiff_t);
        
        // add the elements and nodes to the panel mesh
        // the y-coordinate for this mesh will be the x-coordinate
        // subdomain id for each x-stiffener is 1+i
        _combine_mesh(mesh,
                      stiff,
                      MAST::StiffenedPanelMesh::STIFFENER_X,
                      width*(i+1.)/(n_stiff+1.),
                      i+1);
    }
    
    mesh.prepare_for_use();
}




void
MAST::StiffenedPanelMesh::_combine_mesh(libMesh::MeshBase& panel,
                                        libMesh::MeshBase& component,
                                        MAST::StiffenedPanelMesh::Component c,
                                        Real stiff_offset,
                                        libMesh::subdomain_id_type sid) {
    
    
    libMesh::BoundaryInfo
    &panel_binfo        = *panel.boundary_info,
    &component_binfo    = *component.boundary_info;
    
    panel.reserve_nodes(component.n_nodes());
    panel.reserve_elem (component.n_elem());
    
    libMesh::MeshBase::const_element_iterator
    el_it = component.elements_begin(),
    el_end = component.elements_end();
    
    std::map<libMesh::Node*, libMesh::Node*> old_to_new;
    libMesh::Node *old_node, *new_node;
    
    for ( ; el_it != el_end; el_it++ ) {
        libMesh::Elem* old_elem = *el_it;
        
        libMesh::Elem *new_elem =
        panel.add_elem(libMesh::Elem::build(old_elem->type()).release());
        
        new_elem->subdomain_id() = sid;
        
        // add boundary condition tags for the panel boundary
        for (unsigned short int n=0; n<old_elem->n_sides(); n++) {
            
            if (component_binfo.n_boundary_ids(old_elem, n)) {
                
                // add the boundary tags to the panel mesh
                std::vector<libMesh::boundary_id_type> bc_ids;
                component_binfo.boundary_ids(old_elem, n, bc_ids);
                
                for ( unsigned int bid=0; bid < bc_ids.size(); bid++) {
                    if (c == MAST::StiffenedPanelMesh::PANEL)
                        panel_binfo.add_side(new_elem, n, bc_ids[bid]);
                    else if (c == MAST::StiffenedPanelMesh::STIFFENER_X &
                             bc_ids[bid] != 0 &
                             bc_ids[bid] != 2 )
                        panel_binfo.add_side(new_elem, n, bc_ids[bid]+4);
                }
            }
        }
        
        
        for (unsigned int n=0; n<old_elem->n_nodes(); n++) {
            old_node = old_elem->node_ptr(n);
            
            if (!old_to_new.count(old_node)) {
                libMesh::Point p;
                switch (c) {
                    case MAST::StiffenedPanelMesh::PANEL:
                        p = (*old_node);
                        break;
                        
                    case MAST::StiffenedPanelMesh::STIFFENER_X:
                        p(0) = (*old_node)(0);
                        p(1) = stiff_offset;
                        p(2) = (*old_node)(1);
                        break;
                        
                    case MAST::StiffenedPanelMesh::STIFFENER_Y:
                        p(0) = stiff_offset;
                        p(1) = (*old_node)(0);
                        p(2) = (*old_node)(1);
                        break;
                        
                    default:
                        libmesh_error();
                }
                
                // first see if the node can be found on the panel
                // this will be done only when z-coordinate is zero
                if ((c != MAST::StiffenedPanelMesh::PANEL) &&
                    p(2) == 0.) {
                    
                    libMesh::MeshBase::node_iterator
                    n_it = panel.nodes_begin(),
                    n_end = panel.nodes_end();
                    
                    for ( ; n_it != n_end; n_it++)
                        if ((**n_it)(0) == p(0) &&
                            (**n_it)(1) == p(1)) {
                            
                            old_to_new[old_node] = *n_it;
                            break;
                        }
                }
                else
                    old_to_new[old_node] = panel.add_point(p);
            }
            
            new_node = old_to_new[old_node];
            new_elem->set_node(n) = new_node;
        }
    }
}




