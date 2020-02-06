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
#include "level_set/material_patch.h"
#include "base/field_function_base.h"


MAST::MaterialPatch::MaterialPatch():
_initialized (false) {
    
}


MAST::MaterialPatch::~MaterialPatch() {
    
}


void
MAST::MaterialPatch::init(const libMesh::Elem& elem,
                          const libMesh::Node& node,
                          const MAST::FieldFunction<Real>& phi,
                          const Real t) {

    libmesh_assert(!_initialized);
    
    // make sure that this node is in the region with negative level set
    Real
    v = 0.;
    
    phi(node, t, v);
    libmesh_assert_less_equal(v, 0);

    std::set<const libMesh::Elem*> elems;
    
    // get the list of elements that share this node
    elem.find_point_neighbors(node, elems);
    
    // now that the values of shape function are available, find the nodes
    // in void and in material
    if (_quad4_material_levels(elem, node, elems, phi, t)) {
        
        // now check if any elements are completely in void domain. If so,
        // then remove it from factorization
        std::set<const libMesh::Elem*>::const_iterator
        it   = elems.begin(),
        end  = elems.end();
        
        for ( ; it != end; it++) {
            
            const libMesh::Elem& elem = **it;
            
            bool negative_phi = true;
            for (unsigned int i=0; i<elem.n_nodes(); i++) {
                
                phi(*elem.node_ptr(i), t, v);
                
                if (v >= 0.) {
                    negative_phi = false;
                    break;
                }
            }
            
            // if all nodes on the element are in the void, remove it from the set
            if (!negative_phi)
                _elems_to_factor.insert(&elem);
        }
    }

    _initialized  = true;
}


void
MAST::MaterialPatch::clear() {
    
    _elems_to_factor.clear();
    _initialized = false;
}


bool
MAST::MaterialPatch::_quad4_material_levels(const libMesh::Elem& elem,
                                            const libMesh::Node& node,
                                            const std::set<const libMesh::Elem*>& elem_neighbors,
                                            const MAST::FieldFunction<Real>& phi,
                                            const Real t ) {
    
    // assume that node is at the origin. Then, we can identify the
    // elements to belong to individual quadrants based on the relative
    // location of their centroid
    
    /*
     *      6             5             4
     *        o-----------o-----------o
     *        |           |           |
     *        |    3      |     2     |
     *        |           |           |
     *      7 o-----------o-----------o 3
     *        |           |           |
     *        |     0     |     1     |
     *        |           |           |
     *        o-----------o-----------o
     *       0            1             2
     */
    
    
    // if fewer than 4 elements are identified as neighbors then this node
    // is on the boundary and we ignore it.
    if (elem_neighbors.size() < 4)
        return false;
    
    std::vector<const libMesh::Elem*>
    patch_elems(4, nullptr);
    
    std::set<const libMesh::Elem*>::const_iterator
    e_it   = elem_neighbors.begin(),
    e_end  = elem_neighbors.end();
    
    std::map<Real, const libMesh::Elem*>
    elem_angle_map;
    
    for ( ; e_it != e_end; e_it++) {
        
        const libMesh::Elem& e = **e_it;
        //MAST::plot_elem(e);
        // currently only implemented for quad4 and quad9 elements
        libmesh_assert((e.type() == libMesh::QUAD4) || (e.type() == libMesh::QUAD9));
        
        // position vector from node to element centroid
        const libMesh::Point
        p = e.centroid() - node;

        elem_angle_map[atan2(p(1), p(0))] = &e;
    }

    // now that we have the elements ordered in terms of ascending angle
    // with the x-axis, we should set them in the patch_elems vector
    std::map<Real, const libMesh::Elem*>::const_iterator
    el_it  = elem_angle_map.begin(),
    el_end = elem_angle_map.end();
    
    unsigned int i=0;
    for ( ; el_it != el_end; el_it++) {
        patch_elems[i] = el_it->second;
        i++;
    }
    
    // now identify the nodes on the outer periphery of the patch
    // starting with node 0 of the element in the third quadrant.
    unsigned int
    n_periphery_nodes = 8,
    n_sign_changes    = 0,
    n_material_levels = 0;

    std::vector<const libMesh::Node*>
    periphery_nodes;
    periphery_nodes.reserve(n_periphery_nodes);

    libMesh::Point
    p1, p2;
    
    // bottom-left corner
    periphery_nodes.push_back(patch_elems[0]->node_ptr(0));

    // bottom
    p1 = patch_elems[0]->point(1) - node;
    p2 = patch_elems[1]->point(0) - node;
    if (p1.norm() < p2.norm())
        periphery_nodes.push_back(patch_elems[0]->node_ptr(1));
    else
        periphery_nodes.push_back(patch_elems[1]->node_ptr(0));

    // bottom-right corner
    periphery_nodes.push_back(patch_elems[1]->node_ptr(1));

    // right
    p1 = patch_elems[1]->point(2) - node;
    p2 = patch_elems[2]->point(1) - node;
    if (p1.norm() < p2.norm())
        periphery_nodes.push_back(patch_elems[1]->node_ptr(2));
    else
        periphery_nodes.push_back(patch_elems[2]->node_ptr(1));
    
    // top-right corner
    periphery_nodes.push_back(patch_elems[2]->node_ptr(2));
    
    // top
    p1 = patch_elems[2]->point(3) - node;
    p2 = patch_elems[3]->point(2) - node;
    if (p1.norm() < p2.norm())
        periphery_nodes.push_back(patch_elems[2]->node_ptr(3));
    else
        periphery_nodes.push_back(patch_elems[3]->node_ptr(2));
    
    // top left corner
    periphery_nodes.push_back(patch_elems[3]->node_ptr(3));
    
    // left
    p1 = patch_elems[3]->point(0) - node;
    p2 = patch_elems[0]->point(3) - node;
    if (p1.norm() < p2.norm())
        periphery_nodes.push_back(patch_elems[3]->node_ptr(0));
    else
        periphery_nodes.push_back(patch_elems[0]->node_ptr(3));

    n_periphery_nodes = periphery_nodes.size();
    
    libmesh_assert_equal_to(n_periphery_nodes, 8);
    
    // now, look at the number of sign changes across this periphery nodes
    // It should be an even number nad the number of levels is half the
    // number of sign changes
    
    Real
    v1   = 0.,
    v2   = 0.;
    
    phi(*periphery_nodes[0], t, v1);
    
    for (unsigned int i=0; i<n_periphery_nodes; i++) {
        
        phi(*periphery_nodes[(i+1)%n_periphery_nodes], t, v2);
        
        if ((v1  < 0. && v2 >= 0.) || (v1 >= 0. && v2  < 0.))
            n_sign_changes++;

        v1 = v2;
    }

    // make sure that the number of sign changes is a multiple of 2
    libmesh_assert_equal_to(n_sign_changes%2, 0);
    
    // number of levels is equal to half the number of sign changes
    n_material_levels = n_sign_changes/2;
    
    // presently, we simplify the process by requiring that if a node
    // has multiple material levels, then all involved elements should factor
    // their matrices for assembly.
    
    if (n_material_levels > 1)
        return true;
    else
        return false;
}

