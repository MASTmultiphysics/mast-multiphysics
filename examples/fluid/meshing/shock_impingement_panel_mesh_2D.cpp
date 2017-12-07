///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2017  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//
//// MAST includes
//#include "examples/fluid/meshing/shock_impingement_panel_mesh_2D.h"
//
//// libMesh includes
//#include "libmesh/mesh_serializer.h"
//#include "libmesh/parallel_mesh.h"
//
//
//void
//MAST::ShockImpingementPanel2D::init (const unsigned int n_panels,
//                                     const Real mach,
//                                     const Real panel_length,
//                                     const Real shock_location_on_panel,
//                                     const unsigned int panel_bc_id,
//                                     const unsigned int symmetry_bc_id,
//                                     const std::vector<MeshInitializer::CoordinateDivisions*>& divs,
//                                     libMesh::UnstructuredMesh& mesh, libMesh::ElemType t) {
//    
//    libmesh_assert(divs.size() == 2);
//    libmesh_assert(divs[0]->n_divs() == 2);
//    libmesh_assert(divs[1]->n_divs() >= 1);
//    libmesh_assert(ramp_bc_id > 3);
//    libmesh_assert(symmetry_bc_id > 3);
//    
//    _h_by_l = h_by_l;
//    _x0 = divs[0]->div_location(1);
//    _x1 = divs[0]->div_location(2);
//    _y0 = divs[1]->div_location(0);
//    _y1 = divs[1]->div_location(1);
//    _ramp_bc_id = ramp_bc_id;
//    _symmetry_bc_id = symmetry_bc_id;
//    MeshInitializer::init(divs, mesh, t);
//}
//
//
//
//void
//MAST::ShockImpingementPanel2D::process_mesh( ) {
//    
//    // note that we apply the boundary conditions before moving the mesh since
//    // the application of boudnary conditions is contingent upon the panel surface
//    // points lying on the y = _y0, and the mesh movement ends up altering that
//    
//    // check if the mesh is parallel
//    const bool
//    parallel_mesh = !_mesh->is_serial();
//    
//    {
//        libMesh::MeshSerializer serializer(*_mesh);
//        
//        
//        //march over all the elmeents and tag the sides that all lie on the panel suface
//        libMesh::MeshBase::element_iterator e_it = _mesh->elements_begin();
//        
//        for ( ; e_it != _mesh->elements_end(); e_it++)
//        {
//            // iterate over the sides of each element and check if all
//            // nodes satisfy the requirement
//            
//            for (unsigned int i_side=0; i_side<(*e_it)->n_sides(); i_side++)
//            {
//                libMesh::UniquePtr<libMesh::Elem> side_elem ((*e_it)->side(i_side).release());
//                std::vector<bool> side_on_ramp(side_elem->n_nodes()),
//                side_on_slip_wall(side_elem->n_nodes());
//                std::fill(side_on_ramp.begin(), side_on_ramp.end(), false);
//                std::fill(side_on_slip_wall.begin(), side_on_slip_wall.end(), false);
//                
//                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
//                {
//                    const libMesh::Node& n = *(side_elem->get_node(i_node));
//                    if ((n(1)==_y0) && (n(0) >= _x0-1.0e-6))
//                        side_on_ramp[i_node] = true;
//                    
//                    if ((n(1)==_y0) && ((n(0) <= _x0+1.0e-6)))
//                        side_on_slip_wall[i_node] = true;
//                }
//                
//                // check for side on ramp
//                bool if_apply_bc = true;
//                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
//                    if_apply_bc = side_on_ramp[i_node] && if_apply_bc;
//                if (if_apply_bc) {
//                    _mesh->boundary_info->add_side(*e_it, i_side, _ramp_bc_id);
//                    if (parallel_mesh)
//                        dynamic_cast<libMesh::DistributedMesh*>(_mesh)->add_extra_ghost_elem(*e_it);
//                }
//                
//                // now check for the slip wall
//                if_apply_bc = true;
//                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
//                    if_apply_bc = side_on_slip_wall[i_node] && if_apply_bc;
//                if (if_apply_bc)
//                    _mesh->boundary_info->add_side(*e_it, i_side, _symmetry_bc_id);
//            }
//        }
//        
//        // set the boudnary id names
//        _mesh->boundary_info->sideset_name(_ramp_bc_id)     = "Ramp";
//        _mesh->boundary_info->sideset_name(_symmetry_bc_id) = "Symmetry";
//    }
//    
//    // now move the mesh points
//    libMesh::MeshBase::node_iterator   n_it  = _mesh->nodes_begin();
//    const libMesh::MeshBase::node_iterator n_end = _mesh->nodes_end();
//    
//    Real x, y;
//    
//    for (; n_it != n_end; n_it++)
//    {
//        libMesh::Node& n =  **n_it;
//        
//        // this is for sine bump
//        if ((n(0) >= _x0) && (n(0) <= _x1))
//        {
//            x = n(0);
//            y = n(1);
//            
//            n(1) += _h_by_l * ( 1. -(y-_y0)/(_y1-_y0) ) * (x-_x0)/(_x1-_x0);
//        }
//    }
//}
//
//
//
