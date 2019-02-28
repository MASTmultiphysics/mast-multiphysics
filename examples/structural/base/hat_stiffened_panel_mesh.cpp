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
#include "examples/structural/base/hat_stiffened_panel_mesh.h"
#include "examples/base/multilinear_interpolation.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/material_property_card_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"


// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"




void
MAST::HatStiffenedPanelMesh::init(const unsigned int n_stiff,
                                  const unsigned int n_x_divs,
                                  const unsigned int n_y_divs_on_stiffeners,
                                  const unsigned int n_y_divs_between_stiffeners,
                                  const Real length,
                                  const Real width,
                                  const Real skin_dip_amplitude_by_panel_w,
                                  const Real hat_dip_amplitude_by_panel_w,
                                  const Real stiff_w_by_panel_w,
                                  const Real hat_w_by_stiff_w,
                                  const Real hat_h_by_panel_w,
                                  libMesh::MeshBase& mesh,
                                  libMesh::ElemType t) {
    
    libmesh_assert_less(stiff_w_by_panel_w*n_stiff, 1.);
    
    
    mesh.set_mesh_dimension(2);
    
    
    // initialize the main mesh
    libMesh::SerialMesh panel(mesh.comm());
    
    // create the mesh for each panel
    const unsigned int
    n_y_per_stiff = n_y_divs_between_stiffeners + n_y_divs_on_stiffeners,
    n_y_divs      = n_y_per_stiff*n_stiff + n_y_divs_between_stiffeners;
    
    _create_panel(panel,
                  n_stiff,
                  n_x_divs,
                  n_y_divs_between_stiffeners,
                  n_y_divs_on_stiffeners,
                  length,
                  width,
                  stiff_w_by_panel_w,
                  hat_w_by_stiff_w,
                  skin_dip_amplitude_by_panel_w,
                  t);
    
    // use subdomain id for panel as 0
    _combine_mesh(mesh, panel, MAST::HatStiffenedPanelMesh::PANEL, 0., 0,
                  0.);
    
    
    // now iterate over the stiffeners and create them
    for (unsigned int i=0; i<n_stiff; i++) {
        
        libMesh::SerialMesh stiff(mesh.comm());
        
        _create_hat_stiff(stiff,
                          n_x_divs,
                          n_y_divs_on_stiffeners,
                          length,
                          width,
                          stiff_w_by_panel_w,
                          hat_w_by_stiff_w,
                          hat_h_by_panel_w,
                          hat_dip_amplitude_by_panel_w,
                          t);
        
        // add the elements and nodes to the panel mesh
        // subdomain id for each x-stiffener is 1+i
        Real
        dy_between_stiff = (1.-n_stiff*stiff_w_by_panel_w)/(n_stiff+1)*width,
        stiff_y_offset = dy_between_stiff*(i+1)+i*stiff_w_by_panel_w*width;
        
        _combine_mesh(mesh,
                      stiff,
                      MAST::HatStiffenedPanelMesh::STIFFENER_X,
                      stiff_y_offset,
                      i+1,
                      length*1.e-8);
    }
    
    mesh.prepare_for_use();
}



void
MAST::HatStiffenedPanelMesh::
_create_panel(libMesh::MeshBase& panel,
              const unsigned int n_stiff,
              const unsigned int n_x_divs,
              const unsigned int n_y_divs_between_stiffeners,
              const unsigned int n_y_divs_on_stiffeners,
              const Real length,
              const Real width,
              const Real stiff_w_by_panel_w,
              const Real hat_w_by_stiff_w,
              const Real skin_dip_amplitude_by_panel_w,
              const libMesh::ElemType t) {
    
    // the geometry numbers
    const Real
    w_stiff          = stiff_w_by_panel_w * width,
    dw_between_stiff = (1-n_stiff*stiff_w_by_panel_w)*width/(n_stiff+1),
    w_hat            = w_stiff * hat_w_by_stiff_w;
    
    const unsigned int
    n_y_divs = n_y_divs_on_stiffeners*n_stiff +
    n_y_divs_between_stiffeners*(n_stiff+1);
    
    // create multilinear interpolation function for y_location vs node num
    // in the y-direction
    std::vector<MAST::Parameter*> y_vals(2+2*n_stiff);
    std::vector<MAST::ConstantFieldFunction*> y_vals_f(2+2*n_stiff);
    std::map<Real, MAST::FieldFunction<Real>*> vals;
    
    for (unsigned int i=0; i<2+2*n_stiff; i++) {
        
        MAST::Parameter
        *p = new MAST::Parameter("y", 0.);
        MAST::ConstantFieldFunction
        *f = new MAST::ConstantFieldFunction("y", *p);
        
        y_vals[i]   = p;
        y_vals_f[i] = f;
    }
    
    (*y_vals[2*n_stiff+1])() = width;   // panel edge
    
    // create map of ndv and corresponding location
    vals[0]        = y_vals_f[0];            // panel edge
    vals[n_y_divs] = y_vals_f[1+2*n_stiff];  // panel edge
    for (unsigned int i=0; i<n_stiff; i++) {
        
        (*y_vals[2*i+1])() = (i+1)*dw_between_stiff+i*w_stiff;
        (*y_vals[2*i+2])() = (i+1)*(dw_between_stiff+w_stiff);
        
        vals[(i+1)*n_y_divs_between_stiffeners +
             i*n_y_divs_on_stiffeners] = y_vals_f[2*i+1];
        vals[(i+1)*(n_y_divs_between_stiffeners +
                    n_y_divs_on_stiffeners)] = y_vals_f[2*i+2];
    }
    
    MAST::MultilinearInterpolation y_loc("y", vals);
    
    libMesh::SerialMesh& panel_m = dynamic_cast<libMesh::SerialMesh&>(panel);
    
    // the reference panel is divided into unit width for each stiffener
    // and unit width for each region between two stiffeners or between
    // stiffeners and panel boundary
    libMesh::MeshTools::Generation::build_square(panel_m,
                                                 n_x_divs,
                                                 n_y_divs,
                                                 0., length,
                                                 0., n_y_divs,
                                                 t);
    
    
    
    // now iterate over all the nodes, and move the z-location based on
    // the specified input
    libMesh::MeshBase::node_iterator
    n_it  = panel.nodes_begin(),
    n_end = panel.nodes_end();
    
    Real
    y  = 0.,
    y0 = 0.,
    y1 = 0.,
    pi = acos(-1.);
    
    libMesh::Point
    p;
    
    
    for ( ; n_it != n_end; n_it++) {
        
        libMesh::Node& n = **n_it;
        
        p(0) = n(1);
        y_loc(p, 0., y);
        n(1) = y;
        
        // move the z-location if the node belongs to a stiffener
        for (unsigned int i=0; i<n_stiff; i++) {
            
            y0 = (*y_vals[2*i+1])();// + (w_stiff-w_hat)/2.;
            y1 = (*y_vals[2*i+2])();// - (w_stiff-w_hat)/2.;
            
            if ((y >= y0) && (y <= y1))
                n(2) += skin_dip_amplitude_by_panel_w * width * 0.5 *
                (1.-cos((y-y0)/(y1-y0)*2*pi));
        }
    }
    
    // now delete the constant function pointers
    for (unsigned int i=0; i<2+2*n_stiff; i++) {
        
        delete y_vals_f[i];
        delete y_vals[i];
    }
}



void
MAST::HatStiffenedPanelMesh::_create_hat_stiff(libMesh::MeshBase& stiff,
                                               const unsigned int n_x_divs,
                                               const unsigned int n_y_divs_on_stiffeners,
                                               const Real length,
                                               const Real width,
                                               const Real stiff_w_by_panel_w,
                                               const Real hat_w_by_stiff_w,
                                               const Real hat_h_by_panel_w,
                                               const Real hat_dip_amplitude_by_panel_w,
                                               const libMesh::ElemType t) {
    
    // the geometry numbers
    const Real
    w_stiff = stiff_w_by_panel_w * width,
    hat_y0  = 0.5 * stiff_w_by_panel_w*(1.-hat_w_by_stiff_w) * width,
    hat_y1  = w_stiff - hat_y0,
    hat_z   = hat_h_by_panel_w * width;
    
    libMesh::SerialMesh& stiff_m = dynamic_cast<libMesh::SerialMesh&>(stiff);
    
    libMesh::MeshTools::Generation::build_square(stiff_m,
                                                 n_x_divs,
                                                 n_y_divs_on_stiffeners,
                                                 0., length,
                                                 0., w_stiff,
                                                 t);
    
    
    // now iterate over all the nodes, and move the z-location based on
    // the specified input
    libMesh::MeshBase::node_iterator
    n_it  = stiff.nodes_begin(),
    n_end = stiff.nodes_end();
    
    Real
    y = 0.,
    pi = acos(-1.);
    
    
    for ( ; n_it != n_end; n_it++) {
        
        libMesh::Node& n = **n_it;
        
        y = n(1);
        
        if (y <= hat_y0)
            n(2) = (y/hat_y0) * hat_z;
        else if (y <= hat_y1)
            n(2) = hat_z + hat_dip_amplitude_by_panel_w * width * 0.5 *
            (1. - cos((y-hat_y0)/(hat_y1-hat_y0)*2*pi));
        else
            n(2) = (1.-(y-hat_y1)/hat_y0) * hat_z;
    }
}



void
MAST::HatStiffenedPanelMesh::_combine_mesh(libMesh::MeshBase& panel,
                                           libMesh::MeshBase& component,
                                           MAST::HatStiffenedPanelMesh::Component c,
                                           Real stiff_offset,
                                           libMesh::subdomain_id_type sid,
                                           const Real tol) {
    
    
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
                    if (c == MAST::HatStiffenedPanelMesh::PANEL)
                        panel_binfo.add_side(new_elem, n, bc_ids[bid]);
                    else if (c == MAST::HatStiffenedPanelMesh::STIFFENER_X &
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
                    case MAST::HatStiffenedPanelMesh::PANEL:
                        p = (*old_node);
                        break;
                        
                    case MAST::HatStiffenedPanelMesh::STIFFENER_X: {
                        
                        p = (*old_node);
                        p(1) += stiff_offset;
                    }
                        break;
                        
                    case MAST::HatStiffenedPanelMesh::STIFFENER_Y:
                    default:
                        libmesh_error();
                }
                
                // first see if the node can be found on the panel
                // this will be done only when z-coordinate is zero
                if ((c != MAST::HatStiffenedPanelMesh::PANEL) &&
                    p(2) == 0.) {
                    
                    libMesh::MeshBase::node_iterator
                    n_it = panel.nodes_begin(),
                    n_end = panel.nodes_end();
                    
                    for ( ; n_it != n_end; n_it++)
                        if (fabs((**n_it)(0) - p(0)) < tol &&
                            fabs((**n_it)(1) - p(1)) < tol) {
                            
                            old_to_new[old_node] = *n_it;
                            break;
                        }
                }
                else
                    old_to_new[old_node] = panel.add_point(p);
            }
            
            libmesh_assert(old_to_new.count(old_node));
            new_node = old_to_new[old_node];
            new_elem->set_node(n) = new_node;
        }
    }
}




