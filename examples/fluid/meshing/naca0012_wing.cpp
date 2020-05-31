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


// boost includes
#include "boost/tuple/tuple.hpp"
#include <iomanip>

// MAST includes
#include "examples/fluid/meshing/naca0012_wing.h"
#include "examples/fluid/meshing/mesh_initializer.h"

// libMesh includes
#include "libmesh/mesh_serializer.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/boundary_info.h"


MAST::Examples::NACA0012WingMesh3D::NACA0012WingMesh3D() {
    
    
}


MAST::Examples::NACA0012WingMesh3D::~NACA0012WingMesh3D() {
    
    
}


void
MAST::Examples::NACA0012WingMesh3D::mesh(const Real root_chord,
                                         const Real taper_ratio,
                                         const Real mid_chord_sweep,
                                         const Real far_field_radius_to_root_chord,
                                         const Real span,
                                         const Real spanwise_farfield,
                                         const unsigned int radial_divs_chord,
                                         const unsigned int radial_divs_chord_to_farfield,
                                         const unsigned int quarter_divs,
                                         const unsigned int spanwise_divs,
                                         const unsigned int span_to_farfield_divs,
                                         const Real radial_elem_size_ratio,
                                         const Real spanwise_elem_size_ratio,
                                         libMesh::UnstructuredMesh& mesh,
                                         libMesh::ElemType etype) {
    
    
    Real
    eta   = 0.,
    phi   = 0.,
    pi    = acos(-1.);
    
    libMesh::MeshTools::Generation::build_cube
    (mesh,
     radial_divs_chord+radial_divs_chord_to_farfield,
     4 * quarter_divs,
     spanwise_divs+span_to_farfield_divs,
     0., 1., // from airfoil center to outer boundary
     0., 1., // from 0 to 2*pi around the geometry
     0., 1., // from root to spanwise farfield boundary
     etype);
    
    
    //create a map of old to new node
    std::map<libMesh::Node*, libMesh::Node*>
    node_map_at_y,
    node_map_at_x0;
    
    std::map<Real, libMesh::Node*>
    nodes_at_x0;
    
    std::map<Real, std::map<Real, libMesh::Node*>>
    nodes_to_replace_at_x0,
    nodes_to_use,
    nodes_to_delete;
    
    std::set<libMesh::Elem*>
    elems_to_delete;
    std::set<libMesh::Node*>
    nodes_in_wing_to_delete;

    libMesh::MeshSerializer serializer(mesh);

    // The mesh block is created between [[0,1],[0,2],[0,1]]. This is
    // mapped using to two sets of coordinates:
    // the first set is to identify the interior of the wing, so that
    // those elements can be removed.
    // the second set is the physical location of the wing.
    {
        std::vector<Real>
        x_div_loc      = {0., 1., 2.}, // [0,1] along chord, [1, 2] from surface to farfield
        y_div_loc      = {0., .5, 1., 1.5, 2.}, // from TE to LE to TE
        z_div_loc      = {0., 1., 2.}, // [0, 1] from root to tip, [1, 2] from tip to spanwise farfield
        dx_relative    = std::vector<Real>(3, 1.),
        dy_relative    = std::vector<Real>(5, 1.),
        dz_relative    = std::vector<Real>(3, 1.);
        dx_relative  = {1., 1., radial_elem_size_ratio},
        dz_relative  = {1., 1., spanwise_elem_size_ratio};
        std::vector<unsigned int>
        n_dx         = {radial_divs_chord, radial_divs_chord_to_farfield},
        n_dy         = std::vector<unsigned int>(4, quarter_divs),
        n_dz         = {spanwise_divs, span_to_farfield_divs};
        
        MAST::MeshInitializer::CoordinateDivisions
        xdivs,
        ydivs,
        zdivs;
        xdivs.init(2, x_div_loc, dx_relative, n_dx);
        ydivs.init(4, y_div_loc, dy_relative, n_dy);
        zdivs.init(2, z_div_loc, dz_relative, n_dz);
        
        libMesh::MeshBase::node_iterator
        n_it    = mesh.nodes_begin(),
        n_end   = mesh.nodes_end();
        
        Real
        tol = 1.e-8*root_chord;
        
        {
            for ( ; n_it != n_end; n_it++) {
                
                libMesh::Node& n = **n_it;
                
                // nodes at y=2 will be replaced by nodes at y=0 so that we
                // can get a closed domain.
                if (ydivs(n(1)) < tol)
                    nodes_to_use[n(0)][n(2)] = &n;
                else if (fabs(ydivs(n(1)) - 2.) < tol)
                    nodes_to_delete[n(0)][n(2)] = &n;
                
                // nodes at x=0, y=0 will be used to replace nodes at x=0 at all
                // y values since we will replace these hex8 elements with prism6
                if (n(0) < tol && ydivs(n(1)) < tol)
                    nodes_at_x0[n(2)] = &n;
                else if (n(0) < tol)
                    nodes_to_replace_at_x0[ydivs(n(1))][n(2)] = &n;
            }
            
            
            {
                libmesh_assert_equal_to(nodes_to_use.size(), nodes_to_delete.size());
                
                // this is assuming that the sorted nodes will have the same
                // sequence in both maps, since they are both sorted based on the
                // x, z coordinates that are same on the two planes.
                std::map<Real, std::map<Real, libMesh::Node*>>::iterator
                n1_it  = nodes_to_use.begin(),
                n2_it  = nodes_to_delete.begin();
                
                for ( ; n2_it != nodes_to_delete.end(); n2_it++) {
                    
                    libmesh_assert_equal_to(n2_it->second.size(), n1_it->second.size());
                    
                    std::map<Real, libMesh::Node*>::iterator
                    n1_sub_it = n1_it->second.begin(),
                    n2_sub_it = n2_it->second.begin();
                    
                    for ( ; n2_sub_it != n2_it->second.end(); n2_sub_it++) {
                        
                        node_map_at_y[n2_sub_it->second] = n1_sub_it->second;
                        n1_sub_it++;
                    }
                    
                    n1_it++;
                }
            }
            
            // we use the same logic to identify nodes for replacement on
            // the x=0 plane
            {
                std::map<Real, std::map<Real, libMesh::Node*>>::iterator
                n2_it  = nodes_to_replace_at_x0.begin();
                
                for ( ; n2_it != nodes_to_replace_at_x0.end(); n2_it++) {
                    
                    libmesh_assert_equal_to(n2_it->second.size(), nodes_at_x0.size());
                    
                    std::map<Real, libMesh::Node*>::iterator
                    n1_sub_it  = nodes_at_x0.begin(),
                    n2_sub_it  = n2_it->second.begin();
                    
                    for ( ; n2_sub_it != n2_it->second.end(); n2_sub_it++) {
                        
                        node_map_at_x0[n2_sub_it->second] = n1_sub_it->second;
                        n1_sub_it++;
                    }
                }
                
            }
        }
        
        
        // first identify the elements that will be inside the wing or
        // the ones that need to be converted to prism
        libMesh::MeshBase::element_iterator
        e_it  =  mesh.elements_begin(),
        e_end =  mesh.elements_end();
        
        for ( ; e_it != e_end; e_it++) {
            
            libMesh::Elem& e = **e_it;
            
            libMesh::Point
            c = e.centroid();
            
            if (zdivs(c(2)) < 1. && // 0...1 in the wing, 1...2 till spanwise farfield
                xdivs(c(0)) < 1.) { // 0...1 inside wing, 1...2 will radial farfield
                
                elems_to_delete.insert(&e);
                
                // now identify the nodes to delete from here
                for (unsigned int i=0; i<e.n_nodes(); i++) {
                    
                    libMesh::Node* n = e.node_ptr(i);
                    
                    if (zdivs((*n)(2)) < 1. && xdivs((*n)(0)) < 1.)
                        nodes_in_wing_to_delete.insert(n);
                }
            }
            else if (e.on_boundary()) {
                
                bool convert      = false;
                int  span_side    = -1;
                for (unsigned int i=0; i<e.n_sides(); i++) {
                    
                    // if element has a side on boundary id 4, then we
                    // convert the element
                    if (mesh.boundary_info->has_boundary_id(&e, i, 4))
                        convert   = true;
                    else if (mesh.boundary_info->has_boundary_id(&e, i, 5))
                        // boundary boundary ids, we only need to worry about the
                        // spanwise farfield boundary id, which is side id 5 of hex8,
                        // or side id 4 of the prism6.
                        span_side = i;
                }
                
                if (convert) {
                    
                    libMesh::Elem
                    *prism  = nullptr;
                    
                    switch (e.type()) {
                            
                        case libMesh::HEX8: {
                            
                            // any element with a side on boundary id 4 will be replaced with
                            // a prism
                            /*
                             *   HEX8: 7        6
                             *         o--------z
                             *        /:       /|         zeta
                             *       / :      / |          ^   eta (into page)
                             *    4 /  :   5 /  |          | /
                             *     o--------o   |          |/
                             *     |   o....|...o 2        o---> xi
                             *     |  .3    |  /
                             *     | .      | /
                             *     |.       |/
                             *     o--------o
                             *     0        1
                             *
                             *   libMesh side numbering:
                             *    {0, 3, 2, 1}, // Side 0 : back
                             *    {0, 1, 5, 4}, // Side 1 : bottom
                             *    {1, 2, 6, 5}, // Side 2 : right
                             *    {2, 3, 7, 6}, // Side 3 : top
                             *    {3, 0, 4, 7}, // Side 4 : left
                             *    {4, 5, 6, 7}  // Side 5 : front
                             *
                             *
                             *
                             *  libMesh cube generation with Hex8 puts side 4 on boundary 4, which
                             *  is the side on left (3, 0, 4, 7). This needs to collapse into
                             *  an edge, so (0, 4) of Hex8 will be same as (2, 5) of the Prism6
                             *  and (3, 7) of Hex8 will collapse into (0, 4) of Hex8.
                             *
                             *
                             *   PRISM6:
                             *           5
                             *           o
                             *          /:\
                             *         / : \             zeta
                             *        /  o  \             ^   eta (into page)
                             *     3 o-------o 4          | /
                             *       | . 2 . |            |/
                             *       |.     .|            o---> xi
                             *       o-------o
                             *       0       1
                             *   libMesh side numbering:
                             *     {0, 2, 1, 99}, // Side 0 : back
                             *     {0, 1, 4,  3}, // Side 1 : bottom
                             *     {1, 2, 5,  4}, // Side 2 :
                             *     {2, 0, 3,  5}, // Side 3 : left
                             *     {3, 4, 5, 99}  // Side 4 : front
                             *
                             *
                             *
                             *   Therefore, we will number the prism nodes with the Hex nodes as
                             *   follows:
                             *
                             *   HEX8 -> PRISM6 node mapping:
                             *           4
                             *           o
                             *          /:\
                             *         / : \             zeta
                             *        /  o  \             ^   eta (into page)
                             *     5 o-------o 6          | /
                             *       | . 0 . |            |/
                             *       |.     .|            o---> xi
                             *       o-------o
                             *       1       2
                             *
                             *    Note that the left face of Hex8 has now been removed. The
                             *    remaining faces from Hex8 map to the Prism as
                             *    back   -> 0 (Hex8) -> 0 (Prism6)
                             *    bottom -> 1 (Hex8) -> 3 (Prism6)
                             *    right  -> 2 (Hex8) -> 1 (Prism6)
                             *    top    -> 3 (Hex8) -> 2 (Prism6)
                             *    left   -> 4 (Hex8) -> None
                             *    front  -> 5 (Hex8) -> 4 (Prism6)
                             *
                             */

                            prism = new libMesh::Prism6;
                            prism->set_id(e.id());
                            prism->set_node(0) = e.node_ptr(1);
                            prism->set_node(1) = e.node_ptr(2);
                            libMesh::Node* n0  = e.node_ptr(0);
                            if ((*n0)(1) < tol)
                                prism->set_node(2) = n0; // replace
                            else {
                                libmesh_assert(node_map_at_x0.count(n0));
                                prism->set_node(2) = node_map_at_x0[n0];
                                nodes_in_wing_to_delete.insert(n0);
                            }
                            prism->set_node(3) = e.node_ptr(5);
                            prism->set_node(4) = e.node_ptr(6);
                            libMesh::Node* n4  = e.node_ptr(4);
                            if ((*n4)(1) < tol)
                                prism->set_node(5) = n4; // replace
                            else {
                                libmesh_assert(node_map_at_x0.count(n4));
                                prism->set_node(5) = node_map_at_x0[n4];
                                nodes_in_wing_to_delete.insert(n4);
                            }
                            // nodes 3 and 7 will always be deleted
                            nodes_in_wing_to_delete.insert(e.node_ptr(3));
                            nodes_in_wing_to_delete.insert(e.node_ptr(7));
                            // now replace the elements
                            mesh.delete_elem(&e);
                            mesh.add_elem(prism);
                        }
                            break;
                            
                        case libMesh::HEX27: {
                            
                            // any element with a side on boundary id 4 will be replaced with
                            // a prism
                            /*
                             *   HEX27:     7              18             6
                             *              o--------------o--------------o
                             *             /:             /              /|
                             *            / :            /              / |
                             *           /  :           /              /  |
                             *        19/   :        25/            17/   |
                             *         o--------------o--------------o    |
                             *        /     :        /              /|    |
                             *       /    15o       /    23o       / |  14o
                             *      /       :      /              /  |   /|           zeta
                             *    4/        :   16/             5/   |  / |            ^   eta (into page)
                             *    o--------------o--------------o    | /  |            | /
                             *    |         :    |   26         |    |/   |            |/
                             *    |  24o    :    |    o         |  22o    |            o---> xi
                             *    |         :    |       10     |   /|    |
                             *    |        3o....|.........o....|../.|....o
                             *    |        .     |              | /  |   / 2
                             *    |       .    21|            13|/   |  /
                             * 12 o--------------o--------------o    | /
                             *    |     .        |              |    |/
                             *    |  11o         | 20o          |    o
                             *    |   .          |              |   / 9
                             *    |  .           |              |  /
                             *    | .            |              | /
                             *    |.             |              |/
                             *    o--------------o--------------o
                             *    0              8              1
                             *
                             *   libMesh side numbering:
                             *    {0, 3, 2, 1}, // Side 0 : back
                             *    {0, 1, 5, 4}, // Side 1 : bottom
                             *    {1, 2, 6, 5}, // Side 2 : right
                             *    {2, 3, 7, 6}, // Side 3 : top
                             *    {3, 0, 4, 7}, // Side 4 : left
                             *    {4, 5, 6, 7}  // Side 5 : front
                             *
                             *
                             *
                             *  libMesh cube generation with Hex8 puts side 4 on boundary 4, which
                             *  is the side on left (3, 0, 4, 7). This needs to collapse into
                             *  an edge, so (0, 4) of Hex8 will be same as (2, 5) of the Prism6
                             *  and (3, 7) of Hex8 will collapse into (0, 4) of Hex8.
                             *
                             *
                             *   PRISM15:
                             *               5
                             *               o
                             *              /:\
                             *             / : \
                             *            /  :  \
                             *           /   :   \
                             *       14 o    :    o 13
                             *         /     :     \
                             *        /      :      \
                             *       /       o 11    \
                             *    3 /        :        \4
                             *     o---------o---------o
                             *     |         :12       |
                             *     |         :         |
                             *     |    o    :    o    |            zeta
                             *     |   17    o   16    |             ^   eta (into page)
                             *     |        .2.        |             | /
                             *     |       .   .       |             |/
                             *   9 o      .  o  .      o 10          o---> xi
                             *     |     .  15   .     |
                             *     |  8 o         o 7  |
                             *     |   .           .   |
                             *     |  .             .  |
                             *     | .               . |
                             *     |.                 .|
                             *     o---------o---------o
                             *     0         6         1
                             *   libMesh side numbering:
                             *     {0, 2, 1, 99}, // Side 0 : back
                             *     {0, 1, 4,  3}, // Side 1 : bottom
                             *     {1, 2, 5,  4}, // Side 2 :
                             *     {2, 0, 3,  5}, // Side 3 : left
                             *     {3, 4, 5, 99}  // Side 4 : front
                             *
                             *
                             *
                             *   Therefore, we will number the prism nodes with the Hex nodes as
                             *   follows:
                             *
                             *               4
                             *               o
                             *              /:\
                             *             / : \
                             *            /  :  \
                             *           /   :   \
                             *       16 o    :    o 18
                             *         /     :     \
                             *        /      :      \
                             *       /       o 12    \
                             *    5 /        :        \6
                             *     o---------o---------o
                             *     |         :17       |
                             *     |         :         |
                             *     |    o    :    o    |            zeta
                             *     |   21    o   23    |             ^   eta (into page)
                             *     |        .0.        |             | /
                             *     |       .   .       |             |/
                             *  13 o      .  o  .      o 14          o---> xi
                             *     |     .  22   .     |
                             *     |  8 o         o 10 |
                             *     |   .           .   |
                             *     |  .             .  |
                             *     | .               . |
                             *     |.                 .|
                             *     o---------o---------o
                             *     1         9         2
                             *
                             *
                             *    Note that the left face of Hex8 has now been removed. The
                             *    remaining faces from Hex8 map to the Prism as
                             *    back   -> 0 (Hex8) -> 0 (Prism6)
                             *    bottom -> 1 (Hex8) -> 3 (Prism6)
                             *    right  -> 2 (Hex8) -> 1 (Prism6)
                             *    top    -> 3 (Hex8) -> 2 (Prism6)
                             *    left   -> 4 (Hex8) -> None
                             *    front  -> 5 (Hex8) -> 4 (Prism6)
                             *
                             */

                            prism = new libMesh::Prism18;
                            prism->set_id(e.id());
                            prism->set_node(0) = e.node_ptr( 1);
                            prism->set_node(1) = e.node_ptr( 2);
                            libMesh::Node* n0  = e.node_ptr( 0);
                            if ((*n0)(1) < tol)
                                prism->set_node(2) = n0; // replace
                            else {
                                libmesh_assert(node_map_at_x0.count(n0));
                                prism->set_node(2) = node_map_at_x0[n0];
                                nodes_in_wing_to_delete.insert(n0);
                            }
                            prism->set_node(6) = e.node_ptr( 9);
                            prism->set_node(7) = e.node_ptr(10);
                            prism->set_node(8) = e.node_ptr( 8);

                            
                            prism->set_node( 9) = e.node_ptr(13);
                            prism->set_node(10) = e.node_ptr(14);
                            libMesh::Node* n12  = e.node_ptr(12);
                            if ((*n12)(1) < tol)
                                prism->set_node(11) = n12; // replace
                            else {
                                libmesh_assert(node_map_at_x0.count(n12));
                                prism->set_node(11) = node_map_at_x0[n12];
                                nodes_in_wing_to_delete.insert(n12);
                            }
                            prism->set_node(15) = e.node_ptr(22);
                            prism->set_node(16) = e.node_ptr(23);
                            prism->set_node(17) = e.node_ptr(21);


                            prism->set_node(3) = e.node_ptr(5);
                            prism->set_node(4) = e.node_ptr(6);
                            libMesh::Node* n4  = e.node_ptr(4);
                            if ((*n4)(1) < tol)
                                prism->set_node(5) = n4; // replace
                            else {
                                libmesh_assert(node_map_at_x0.count(n4));
                                prism->set_node(5) = node_map_at_x0[n4];
                                nodes_in_wing_to_delete.insert(n4);
                            }
                            prism->set_node(12) = e.node_ptr(17);
                            prism->set_node(13) = e.node_ptr(18);
                            prism->set_node(14) = e.node_ptr(16);

                            // nodes to be deleted
                            nodes_in_wing_to_delete.insert(e.node_ptr( 3));
                            nodes_in_wing_to_delete.insert(e.node_ptr(15));
                            nodes_in_wing_to_delete.insert(e.node_ptr( 7));
                            
                            nodes_in_wing_to_delete.insert(e.node_ptr(11));
                            nodes_in_wing_to_delete.insert(e.node_ptr(24));
                            nodes_in_wing_to_delete.insert(e.node_ptr(19));
                            
                            nodes_in_wing_to_delete.insert(e.node_ptr(20));
                            nodes_in_wing_to_delete.insert(e.node_ptr(26));
                            nodes_in_wing_to_delete.insert(e.node_ptr(25));


                            // now replace the elements
                            mesh.delete_elem(&e);
                            mesh.add_elem(prism);
                        }
                            break;
                            
                        default:
                            libmesh_error(); // now yet handled yet
                            break;
                    }
                    
                    
                    // set a different subdomain id for the prism
                    // element so that it can be written to the exodus
                    // file
                    prism->subdomain_id() = 1;

                    // boundary boundary ids, we only need to worry about the
                    // spanwise farfield boundary id, which is side id 5 of hex8,
                    // or side id 4 of the prism6.
                    if (span_side > 0)
                        mesh.boundary_info->add_side(prism, 4, 5);
                }
            }
        }
    }
    
    // clear boundary ID 4 so that we can add elements with sides on
    // wing surface to it
    mesh.boundary_info->remove_id(4);
    
    // initialize the mesh divs to distribute the locations
    std::vector<Real>
    x_div_loc      = {0., 1., 2.}, // [0,1] along chord, [1, 2] from surface to farfield
    y_div_loc      = {0., .5, 1., 1.5, 2.}, // from TE to LE to TE
    z_div_loc      = {0., span, spanwise_farfield},
    dx_relative    = std::vector<Real>(3, 1.),
    dy_relative    = std::vector<Real>(5, 1.),
    dz_relative    = std::vector<Real>(3, 1.);
    dx_relative  = {1., 1., radial_elem_size_ratio},
    dz_relative  = {1., 1., spanwise_elem_size_ratio};
    std::vector<unsigned int>
    n_dx         = {radial_divs_chord, radial_divs_chord_to_farfield},
    n_dy         = std::vector<unsigned int>(4, quarter_divs),
    n_dz         = {spanwise_divs, span_to_farfield_divs};
    
    MAST::MeshInitializer::CoordinateDivisions
    xdivs,
    ydivs,
    zdivs;
    xdivs.init(2, x_div_loc, dx_relative, n_dx);
    ydivs.init(4, y_div_loc, dy_relative, n_dy);
    zdivs.init(2, z_div_loc, dz_relative, n_dz);

    
    Real
    tc               =  0.12,
    x                =  0.,
    y                =  0.,
    local_chord      =  0.,
    local_midpt      =  0.,
    z                =  0.,
    far_field_radius = root_chord * far_field_radius_to_root_chord;
    
    libMesh::MeshBase::node_iterator
    n_it    = mesh.nodes_begin(),
    n_end   = mesh.nodes_end();
    
    for ( ; n_it != n_end; n_it++) {
        
        libMesh::Node& n = **n_it;
        
        eta  = xdivs(n(0));
        phi  = ydivs(n(1));
        
        x    = (cos(phi * pi) + 1.)*0.5;
        
        // added a factor of 2 to scale to the correct thickness
        y    = 5.* (2. * tc) * (.2969 * sqrt(x) - 0.1260 * x - 0.3516 * pow(x, 2) + 0.2843 * pow(x, 3) - 0.1036 * pow(x, 4));
        
        if (phi > 1.)
            y *= -1.;

        // sweep, taper, etc. are defined along the z-axis
        z    =  zdivs(n(2));
        n(2) = z;
        
        if (z <= span) {
            local_midpt = z * sin(mid_chord_sweep);
            local_chord = root_chord * (1. - (1. - taper_ratio) * (z/span));
        }
        else {
            local_midpt = span * sin(mid_chord_sweep);
            local_chord = root_chord * taper_ratio;
        }
        
        if (eta <= 1.) {
            // inside the airfoil section
            n(0) =  local_midpt + ((1. - eta) * 0 + eta * local_chord) * cos(phi * pi);
            n(1) = (1. - eta) * 0                  + eta * local_chord * y;
        }
        else {
            // outside the airfoil section
            eta -= 1.;
            n(0) = (1. - eta) * (local_chord * cos(phi * pi) + local_midpt) + eta * far_field_radius * cos(phi * pi);
            n(1) = (1. - eta) * local_chord * y              + eta * far_field_radius * sin(phi * pi);
        }
        
    }
    
    
    // iterate over the elements and reassign the boundary nodes, to close
    // the mesh
    libMesh::MeshBase::element_iterator
    e_it  =  mesh.elements_begin(),
    e_end =  mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        
        libMesh::Elem& e = **e_it;
        
        if (!elems_to_delete.count(&e)) {
            
            for (unsigned int i=0; i<e.n_nodes(); i++)
                if (node_map_at_y.count(e.node_ptr(i)))
                    e.set_node(i) = node_map_at_y[e.node_ptr(i)];
        }
    }
    
    
    // delete the nodes
    std::map<Real, std::map<Real, libMesh::Node*>>::iterator
    nodes_it   = nodes_to_delete.begin();
    
    for (; nodes_it != nodes_to_delete.end(); nodes_it++) {
        
        std::map<Real, libMesh::Node*>::iterator
        nodes_sub_it = nodes_it->second.begin();
        
        for (; nodes_sub_it != nodes_it->second.end(); nodes_sub_it++) {
            
            libMesh::Node *n = nodes_sub_it->second;
            
            if (!nodes_in_wing_to_delete.count(n))
                mesh.delete_node(n);
        }
    }
    
    // now, delete all elements that are inside the wing
    std::set<libMesh::Elem*>::iterator
    el_del_it  = elems_to_delete.begin(),
    el_del_end = elems_to_delete.end();
    
    for ( ; el_del_it != el_del_end; el_del_it++)
        mesh.delete_elem(*el_del_it);

    std::set<libMesh::Node*>::iterator
    nd_del_it  = nodes_in_wing_to_delete.begin(),
    nd_del_end = nodes_in_wing_to_delete.end();
    
    for ( ; nd_del_it != nd_del_end; nd_del_it++)
        mesh.delete_node(*nd_del_it);

    mesh.prepare_for_use();
    
    // now iterate over all the elements and if there is any element
    // with side on boundary that is not currently in the boudnary
    // info, we assume it to be on wing surface and add it to boundary
    // id 4.
    e_it  =  mesh.elements_begin();
    e_end =  mesh.elements_end();
    
    for ( ; e_it != e_end; e_it++) {
        
        libMesh::Elem *e = *e_it;
        
        if (e->on_boundary()) {

            for (unsigned int i=0; i<e->n_sides(); i++) {
                
                if (e->neighbor_ptr(i) == nullptr &&
                    mesh.boundary_info->n_boundary_ids(e, i) == 0)
                    mesh.boundary_info->add_side(e, i, 4);
            }
        }
    }
    
    
    // libMesh numbering:
    //  0  -> back
    //  1  -> bottom
    //  2  -> right
    //  3  -> top
    //  4  -> left
    //  5  -> front
    mesh.boundary_info->sideset_name(0) = "root";     // rename the back face
    mesh.boundary_info->sideset_name(2) = "farfield"; // rename the right face
    mesh.boundary_info->sideset_name(5) = "farfield_span";     // rename the front face
    mesh.boundary_info->sideset_name(4) = "wing";
    mesh.boundary_info->remove_id(1); // remove the bottom face
    mesh.boundary_info->remove_id(3); // remove the top face
}



