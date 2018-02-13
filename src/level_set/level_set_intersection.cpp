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

// C++ includes
#include <memory>
#include <map>
#include <numeric>

// MAST includes
#include "level_set/level_set_intersection.h"
#include "base/field_function_base.h"



MAST::LevelSetIntersection::LevelSetIntersection(unsigned int max_elem_id,
                                                 unsigned int max_node_id):
_tol                             (1.e-8),
_max_iters                       (10),
_max_mesh_elem_id                (max_elem_id),
_max_mesh_node_id                (max_node_id),
_max_elem_divs                   (4),
_initialized                     (false),
_if_elem_on_positive_phi         (false),
_if_elem_on_negative_phi         (false),
_mode                            (MAST::NO_INTERSECTION),
_node_num_on_boundary            (0),
_edge_num_on_boundary            (0){
    
}
    

MAST::LevelSetIntersection::~LevelSetIntersection() {

    this->clear();
}



MAST::LevelSet2DIntersectionMode
MAST::LevelSetIntersection::get_intersection_mode() const {
    
    libmesh_assert(_initialized);
    
    return _mode;
}



unsigned int
MAST::LevelSetIntersection::node_on_boundary() const {
    
    libmesh_assert(_initialized);
    libmesh_assert(_mode == MAST::THROUGH_NODE);
    
    return _node_num_on_boundary;
}



unsigned int
MAST::LevelSetIntersection::edge_on_boundary() const {
    
    libmesh_assert(_initialized);
    libmesh_assert(_mode == MAST::COLINEAR_EDGE);
    
    return _edge_num_on_boundary;
}




bool
MAST::LevelSetIntersection::if_elem_on_positive_phi() const {
    
    libmesh_assert(_initialized);

    return _if_elem_on_positive_phi;
}


bool
MAST::LevelSetIntersection::if_elem_on_negative_phi() const {
    
    libmesh_assert(_initialized);
    
    return _if_elem_on_negative_phi;
}



bool
MAST::LevelSetIntersection::if_elem_has_positive_phi_region() const {

    libmesh_assert(_initialized);
    
    return (_if_elem_on_positive_phi || !_if_elem_on_negative_phi);
}



bool
MAST::LevelSetIntersection::if_elem_has_boundary() const {
    
    return (this->if_intersection_through_elem() ||
            (_mode == MAST::COLINEAR_EDGE));
}




bool
MAST::LevelSetIntersection::if_intersection_through_elem() const {
    
    return ((_mode == MAST::OPPOSITE_EDGES) ||
            (_mode == MAST::ADJACENT_EDGES) ||
            (_mode == MAST::OPPOSITE_NODES) ||
            (_mode == MAST::NODE_AND_EDGE));
}




void
MAST::LevelSetIntersection::clear() {
    
    _tol                        = 1.e-8;
    _max_iters                  = 10;
    _initialized                = false;
    _if_elem_on_positive_phi    = false;
    _if_elem_on_negative_phi    = false;
    _mode                       = MAST::NO_INTERSECTION;
    _node_num_on_boundary       = 0;
    _edge_num_on_boundary       = 0;
    _positive_phi_elems.clear();
    _negative_phi_elems.clear();
    _elem_sides_on_interface.clear();
    _node_local_coords.clear();
    
    std::vector<libMesh::Elem*>::iterator
    e_it  = _new_elems.begin(),
    e_end = _new_elems.end();
    
    for ( ; e_it != e_end; e_it++)
        delete *e_it;
    
    
    std::vector<libMesh::Node*>::iterator
    n_it  = _new_nodes.begin(),
    n_end = _new_nodes.end();
    
    for ( ; n_it != n_end; n_it++)
        delete *n_it;

    _new_nodes.clear();
    _new_elems.clear();
}


void
MAST::LevelSetIntersection::init(const MAST::FieldFunction<Real>& phi,
                                 const libMesh::Elem& e,
                                 const Real t) {
    
    // make sure that this has not been initialized already
    libmesh_assert(!_initialized);
    libmesh_assert_equal_to(e.dim(), 2); // this is only for 2D elements
    
    // this assumes that the phi=0 interface is not fully contained inside
    // the element. That is, the interface is assumed to intersect with the
    // element edges.
    
    // check for each mode of intersection. It is assumed that an element only
    // sees a single mode of intersection, which is hopefully true for a fine
    // mesh. However, with larger high-order elements this assumption will
    // not be valid and more accurate implementations will be needed.
    
    std::map<const libMesh::Node*, std::pair<Real, bool> >
    node_phi_vals;

    unsigned int
    n_node_intersection = 0;

    Real
    val     = 0.,
    min_val =  1.e12,  // set arbitrary high values
    max_val = -1.e12;

    bool
    on_level_set = false;
    
    for (unsigned int i=0; i<e.n_nodes(); i++) {
        
        const libMesh::Node& n = e.node_ref(i);
        phi(n, t, val);
        on_level_set      = fabs(val) <= _tol;
        node_phi_vals[&n] = std::pair<Real, bool>(val, on_level_set);
        if (on_level_set) {
            _node_num_on_boundary = i;
            n_node_intersection++;
        }
        min_val = (min_val > val)? val : min_val;
        max_val = (max_val < val)? val : max_val;
    }


    /////////////////////////////////////////////////////////////////////
    //   Check to see if the whole element is on either side of the function
    /////////////////////////////////////////////////////////////////////

    // if the sign of function on all nodes is the same, then it is assumed
    // that the element is not intersected
    if (min_val > _tol &&
        max_val > _tol) {
        // element is completely on the positive side with no intersection
        
        _mode = MAST::NO_INTERSECTION;
        _positive_phi_elems.push_back(&e);
        _if_elem_on_positive_phi         = true;
        _if_elem_on_negative_phi         = false;
        _initialized = true;
        return;
    }
    else if (min_val < _tol &&
             max_val < _tol) {
        
        // element is completely on the negative side, with no intersection
        
        _mode = MAST::NO_INTERSECTION;
        _negative_phi_elems.push_back(&e);
        _if_elem_on_positive_phi         = false;
        _if_elem_on_negative_phi         = true;
        _initialized = true;
        return;
    }
    else if (min_val < _tol &&
             max_val > _tol) {
        // if it did not get caught in the previous two cases, then there is
        // an intersection.
        // If it got here, then there is an intersection in the domain and the
        // element is not on the positive side completely.
        
        _if_elem_on_positive_phi         = false;
        _if_elem_on_negative_phi         = false;
    }
        
    /////////////////////////////////////////////////////////////////////
    //   Check to see if the function passes through one node
    /////////////////////////////////////////////////////////////////////
    // if only one node intersection was found, then the level set
    // passes through a node. Node intersection happens when the level set on
    // a node is zero and the level set does not change sign on the element.
    Real
    val1   = std::fabs(min_val)<_tol ? 0.: min_val,
    val2   = std::fabs(max_val)<_tol ? 0.: max_val;
    bool
    sign_change = (val1*val2 < 0.);
    if (n_node_intersection == 1 &&
        !sign_change) {
        
        _mode = MAST::THROUGH_NODE;
        
        // this means that the whole element is going to be on either the
        // positive or the negative side of the level-set function
        // now figure out which side of the level-set function the
        // element lies on
        
        if (min_val <= _tol &&
            max_val <= _tol ) { // element is on the negative side
            _negative_phi_elems.push_back(&e);
        }
        else {
            
            libmesh_assert_greater(max_val, _tol);
            _positive_phi_elems.push_back(&e);
        }
        
        std::vector<std::pair<libMesh::Point, libMesh::Point> >
        side_nondim_points;
        _add_node_local_coords(e, side_nondim_points, _node_local_coords);
        _initialized = true;
        return;
    }
    
    
    /////////////////////////////////////////////////////////////////////
    //   Check to see if the function is conliniear with an edge
    /////////////////////////////////////////////////////////////////////
    //  Check to see if all nodes on an edge are on the level set
    for (unsigned int i=0; i<e.n_sides(); i++) {
        
        std::unique_ptr<const libMesh::Elem>
        s(e.side_ptr(i).release());
        
        unsigned int
        n_nodes_on_level_set = 0;
        
        for (unsigned int j=0; j<s->n_nodes(); j++)
            n_nodes_on_level_set += node_phi_vals[s->node_ptr(j)].second;
        
        if (n_nodes_on_level_set == s->n_nodes()) {
            
            // side is on the level set
            _mode                        = MAST::COLINEAR_EDGE;
            _edge_num_on_boundary        = i;
            
            // now that we know that one of the edges is on the function,
            // we need to identify if the element is on the positive or
            // the negative side of the function. It is on the positive side
            // if the maximum value of the function on the nodes is positive.
            // It is on the negative side if the minimum value of the function
            // is negative.
            if (max_val > _tol)
                _positive_phi_elems.push_back(&e);
            else if (min_val < _tol)
                _negative_phi_elems.push_back(&e);
            
            std::vector<std::pair<libMesh::Point, libMesh::Point> >
            side_nondim_points;
            _add_node_local_coords(e, side_nondim_points, _node_local_coords);
            _elem_sides_on_interface[&e] = i;
            _initialized                 = true;
            return;
        }
    }

    /////////////////////////////////////////////////////////////////////
    //   Identify the edges and locations where intersection occurs
    /////////////////////////////////////////////////////////////////////
    // now that we are here, this means that none of the previous modes
    // identify the intersection of the level set with this element.
    // Therefore, we will now identify the sides of the element where
    // intersectio occurs and create the sub-elements for integration.
    // This is done separately for each element type.
    switch (e.type()) {
            
        case libMesh::QUAD4:
            _find_quad4_intersections(phi, e, t, node_phi_vals);
            break;
            
        default:
            libmesh_error_msg("level-set intersections for elem type not handled.");
            break;
    }
    
    // set the IDs of the new elements. Both the subdomain ID and
    // element IDs are set here.
    //
    // Note that the subdomain IDs are needed
    // for querying the BCs, Loads and properties from the physics. We
    // use the same subdomain ID as its parent ID
    //
    // The element id is a little more tricky. Strictly speaking, this
    // should not be necessary, since we are only creating temporary
    // elements to compute element quantities. However, the current stress
    // storage implementation is storing this for element ids. So, we provide
    // a unique element ID for each new element.
    //
    
    for (unsigned int i=0; i<_new_elems.size(); i++) {
        
        _new_elems[i]->subdomain_id() = e.subdomain_id();
        _new_elems[i]->set_id(_max_mesh_elem_id+e.id()*_max_elem_divs+i);
    }
    
    _initialized = true;
}
    



const std::vector<const libMesh::Elem*>&
MAST::LevelSetIntersection::get_sub_elems_positive_phi() const {
    
    return _positive_phi_elems;
}
    


const std::vector<const libMesh::Elem*>&
MAST::LevelSetIntersection::get_sub_elems_negative_phi() const {

    return _negative_phi_elems;
}


bool
MAST::LevelSetIntersection::
has_side_on_interface(const libMesh::Elem& e) const {
    
    std::map<const libMesh::Elem*, int>::const_iterator it;
    it  = _elem_sides_on_interface.find(&e);
    
    libmesh_assert(it != _elem_sides_on_interface.end());
    
    return it->second >= 0;
}



unsigned int
MAST::LevelSetIntersection::
get_side_on_interface(const libMesh::Elem& e) const {
    
    std::map<const libMesh::Elem*, int>::const_iterator it;
    it  = _elem_sides_on_interface.find(&e);
    
    libmesh_assert(it != _elem_sides_on_interface.end());
    libmesh_assert(it->second >= 0);
    
    return it->second;
}


void
MAST::LevelSetIntersection::_add_node_local_coords
( const libMesh::Elem& e,
 std::vector<std::pair<libMesh::Point, libMesh::Point> >& side_nondim_points,
 std::map<const libMesh::Node*, libMesh::Point>& node_coord_map) {
    
    switch (e.type()) {

        case libMesh::QUAD4: {
            
            // We create a vector storing the nondimensional locations of nodes
            // on each side of the element. This is used later to identify the
            // nondimensional location of all new nodes
            side_nondim_points.resize(4);
            // side 0: eta =-1, xi=[-1, 1]
            side_nondim_points[0] =
            std::pair<libMesh::Point, libMesh::Point>
            (libMesh::Point(-1, -1, 0.),
             libMesh::Point(+1, -1, 0.));
            side_nondim_points[1] =
            std::pair<libMesh::Point, libMesh::Point>
            (libMesh::Point(+1, -1, 0.),
             libMesh::Point(+1, +1, 0.));
            side_nondim_points[2] =
            std::pair<libMesh::Point, libMesh::Point>
            (libMesh::Point(+1, +1, 0.),
             libMesh::Point(-1, +1, 0.));
            side_nondim_points[3] =
            std::pair<libMesh::Point, libMesh::Point>
            (libMesh::Point(-1, +1, 0.),
             libMesh::Point(-1, -1, 0.));
            
            // populate the node to nondimensional coordinate map with the
            // original nodes of the element
            for (unsigned int i=0; i<4; i++)
                node_coord_map[e.node_ptr(i)] = side_nondim_points[i].first;
        }
            break;
            
        default:
            libmesh_assert(false); // not handled.
    }
    
}


void
MAST::LevelSetIntersection::_find_quad4_intersections
(const MAST::FieldFunction<Real>& phi,
 const libMesh::Elem& e,
 const Real t,
 const std::map<const libMesh::Node*, std::pair<Real, bool> >&
 node_phi_vals) {

    libmesh_assert_equal_to(e.type(), libMesh::QUAD4);
    libmesh_assert(!_initialized);
    
    std::vector<std::pair<libMesh::Point, libMesh::Point> >
    side_nondim_points;
    _add_node_local_coords(e, side_nondim_points, _node_local_coords);
    
    const unsigned int
    n_nodes = 4;   // this is equal to n_sides for QUAD4
    
    std::vector<std::pair<bool, Real> >
    side_intersection(n_nodes, std::pair<bool, Real>(false, 0.));
    std::vector<bool>
    node_intersection(n_nodes, false);
    
    std::map<const libMesh::Node*, std::pair<Real, bool> >::const_iterator
    it0,
    it1,
    it_end = node_phi_vals.end();
    
    Real
    v0  = 0.,
    v1  = 0.;

    
    for (unsigned int i=0; i<n_nodes; i++) {

        it0 =  node_phi_vals.find(e.node_ptr(i));
        libmesh_assert(it0 != it_end);

        // identify intersection with nodes
        if (std::fabs(it0->second.first) < _tol)
            node_intersection[i] = true;
        else { // if there is node intersection, then there won't be edge interseciton
            
            // identify intersection with edges
            it1 =  node_phi_vals.find(e.node_ptr((i+1)%n_nodes));
            libmesh_assert(it1 != it_end);
            v0  =  it0->second.first;
            v1  =  it1->second.first;
            if (std::fabs(v0) > _tol &&
                std::fabs(v1) > _tol &&
                v0 * v1 < 0.) {
                
                // intersection lies somewhere in the interior of the edge
                side_intersection[i] = std::pair<bool, Real>
                (true,
                 _find_intersection_on_straight_edge(*it0->first,
                                                     *it1->first,
                                                     phi,
                                                     t));
            }
        }
    }
    
    
    
    // currently, we only handle two intersections
    unsigned int
    n_intersection = 0;
    
    // add up intersections on edges
    {
        std::vector<std::pair<bool, Real> >::const_iterator
        it  = side_intersection.begin(),
        end = side_intersection.end();
        
        for (; it != end; it++)
            n_intersection += it->first;
    }

    // add up intersections on nodes
    {
        std::vector<bool>::const_iterator
        it  = node_intersection.begin(),
        end = node_intersection.end();
        
        for (; it != end; it++)
            n_intersection += *it;
    }

    
    libmesh_assert_equal_to(n_intersection, 2);
    
    // now identify which approch to use for creation of the sub-elements.
    // one case is where two adjacent sides of the element have a
    unsigned int
    ref_side = 0,
    ref_node = 0;
    
    for (unsigned int i=0; i<n_nodes; i++) {
        
        if (side_intersection[i].first &&
            side_intersection[(i+1)%n_nodes].first) {
            // found adjacent mode with ref edge i
            _mode    = MAST::ADJACENT_EDGES;
            ref_side = i;
            break;
        }
        else if (side_intersection[i].first &&
                 side_intersection[(i+2)%n_nodes].first) {
            // found adjacent mode with ref edge i
            _mode    = MAST::OPPOSITE_EDGES;
            ref_side = i;
            break;
        }
        else if (node_intersection[i] &&
                 node_intersection[(i+2)%n_nodes]) {
            // found intersection on two opposite nodes
            _mode    = MAST::OPPOSITE_NODES;
            ref_node = i;
            break;
        }
        else if (node_intersection[i] &&
                 side_intersection[(i+1)%n_nodes].first) {
            // found intersection on node and edge
            _mode    = MAST::NODE_AND_EDGE;
            ref_node = i;
            ref_side = i+1;
        }
        else if (node_intersection[i] &&
                 side_intersection[(i+2)%n_nodes].first) {
            // found intersection on node and edge
            _mode    = MAST::NODE_AND_EDGE;
            ref_node = i;
            ref_side = i+2;
        }
    }
    
    // make sure that atleast one of the modes was found
    libmesh_assert(_mode == MAST::ADJACENT_EDGES ||
                   _mode == MAST::OPPOSITE_EDGES ||
                   _mode == MAST::OPPOSITE_NODES ||
                   _mode == MAST::NODE_AND_EDGE);
    
    ///////////////////////////////////////////////////////
    // now create the sub elements
    ///////////////////////////////////////////////////////
    // For intersection with opposite edges, both sides of the element
    // will be represented by QUAD4 elements.
    //
    // When the intersection on adjacent elements happens, then
    // one side will be represented with triangles and the other
    // side is represented with a combination of two quads.
    //
    // for opposite nodes, one triangle is used for each side of the
    // intersection
    //
    // for node and edge, one side is represented with a triangle,
    // and the other with a quad.
    //
    
    Real
    xi_ref   = 0.,
    xi_other = 0.;
    libMesh::Point
    p,
    side_p0,
    side_p1;
    
    // this should pass, but just to be sure, we will check.
    libmesh_assert_equal_to(_new_nodes.size(), 0);
    
    libMesh::Node
    *nd = nullptr;
    
    std::unique_ptr<const libMesh::Elem>
    s;

    bool
    node0_positive = false;

    // used this to set unique node ids, since some of libMesh's operations
    // depend on valid and unique ids for nodes.
    unsigned int
    node_id_incr   = 1;
    
    //////////////////////////////////////////////////////////////////
    // Create a new node at each intersection: first, on ref_side
    //
    // New nodes are needed only when intersection does not exist
    // opposite nodes. In this case, the current nodes are used for
    // the new elements.
    //////////////////////////////////////////////////////////////////
    if (_mode == MAST::ADJACENT_EDGES ||
        _mode == MAST::OPPOSITE_EDGES ||
        _mode == MAST::NODE_AND_EDGE) {
        
        xi_ref = side_intersection[ref_side%n_nodes].second;
        s.reset(e.side_ptr(ref_side%n_nodes).release());
        p  = s->point(0) + xi_ref * (s->point(1) - s->point(0));
        nd = new libMesh::Node(p);
        nd->set_id(_max_mesh_node_id + node_id_incr);
        node_id_incr++;
        _new_nodes.push_back(nd);
        side_p0 = side_nondim_points[ref_side%n_nodes].first;
        side_p1 = side_nondim_points[ref_side%n_nodes].second;
        _node_local_coords[nd] = side_p0 + xi_ref * (side_p1 - side_p0);

        
        //////////////////////////////////////////////////////////////////
        // before we create the new elements, we will figure out
        // which portion of the side is on the positive or negative side of
        // the level set function
        //////////////////////////////////////////////////////////////////

        it0 = node_phi_vals.find(s->node_ptr(0));
        v0  = it0->second.first;
        
        if (v0 > 0.)
            node0_positive = true;
    }
    else {
        
        it0 = node_phi_vals.find(e.node_ptr(ref_node+1));
        v0  = it0->second.first;
        
        if (v0 > 0.)
            node0_positive = true;

    }

    
    
    switch (_mode) {
            
        case MAST::OPPOSITE_EDGES: {
            
            // create a new node at each intersection: next, on ref_side + 2
            xi_other = side_intersection[(ref_side+2)%n_nodes].second;
            s.reset(e.side_ptr((ref_side+2)%n_nodes).release());
            p  = s->point(0) + xi_other * (s->point(1) - s->point(0));
            nd = new libMesh::Node(p);
            nd->set_id(_max_mesh_node_id + node_id_incr);
            node_id_incr++;
            _new_nodes.push_back(nd);
            side_p0 = side_nondim_points[(ref_side+2)%n_nodes].first;
            side_p1 = side_nondim_points[(ref_side+2)%n_nodes].second;
            _node_local_coords[nd] = side_p0 + xi_other * (side_p1 - side_p0);
            
            libMesh::Elem
            *e_p = const_cast<libMesh::Elem*>(&e),
            *e1  = libMesh::Elem::build(libMesh::QUAD4, e_p).release(),
            *e2  = libMesh::Elem::build(libMesh::QUAD4, e_p).release();
            
            _new_elems.push_back(e1);
            _new_elems.push_back(e2);
            
            // create the elements. We create the elements such that:
            //     node(0) = new node on ref_side
            //     node(1) = elem node (ref_side+1)%n_sides
            //     node(2) = elem node (ref_side+2)%n_sides
            //     node(3) = new node on ref_side+2
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident with ref_side of original element
            //     side(1) = coincident side with ref_side+1 of original element
            //     side(2) = coincident side with ref_side+2 of original element
            //     side(3) = coincident side with level-set interface
            
            e1->set_node(0) = _new_nodes[0];
            e1->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+1)%n_nodes));
            e1->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+2)%n_nodes));
            e1->set_node(3) = _new_nodes[1];
            
            _elem_sides_on_interface[e1] = 3;
            
            // create a second element and set its nodes so that:
            //     node(0) = new node on ref_side+2
            //     node(1) = elem node (ref_side+3)%n_nodes
            //     node(2) = elem node (ref_side+4)%n_nodes
            //     node(3) = new node on ref_side
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident with ref_side+2 of original element
            //     side(1) = coincident side with (ref_side+3)%n_sides of original element
            //     side(2) = coincident side with ref_side of original element
            //     side(3) = coincident side with level-set interface
            
            e2->set_node(0) = _new_nodes[1];
            e2->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+3)%n_nodes));
            e2->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr(ref_side));
            e2->set_node(3) = _new_nodes[0];
            
            _elem_sides_on_interface[e2] = 3;
            
            if (node0_positive) {
                
                _positive_phi_elems.push_back(e2);
                _negative_phi_elems.push_back(e1);
            }
            else {
                
                _positive_phi_elems.push_back(e1);
                _negative_phi_elems.push_back(e2);
            }
        }
            break;
            
        case MAST::ADJACENT_EDGES: {
            
            
            // create a new node at each intersection: second on ref_side+1
            xi_other = side_intersection[(ref_side+1)%n_nodes].second;
            s.reset(e.side_ptr((ref_side+1)%n_nodes).release());
            p  = s->point(0) + xi_other * (s->point(1) - s->point(0));
            nd = new libMesh::Node(p);
            nd->set_id(_max_mesh_node_id + node_id_incr);
            node_id_incr++;
            _new_nodes.push_back(nd);
            side_p0 = side_nondim_points[(ref_side+1)%n_nodes].first;
            side_p1 = side_nondim_points[(ref_side+1)%n_nodes].second;
            _node_local_coords[nd] = side_p0 + xi_other * (side_p1 - side_p0);
            
            
            // create a new node on the opposite side, which will be used to create
            // the QUAD4 elements: second on ref_side+2
            s.reset(e.side_ptr((ref_side+2)%n_nodes).release());
            p  = s->point(1) + xi_ref * (s->point(0) - s->point(1));
            nd = new libMesh::Node(p);
            nd->set_id(_max_mesh_node_id + node_id_incr);
            node_id_incr++;
            _new_nodes.push_back(nd);
            side_p0 = side_nondim_points[(ref_side+2)%n_nodes].first;
            side_p1 = side_nondim_points[(ref_side+2)%n_nodes].second;
            _node_local_coords[nd] = side_p1 + xi_ref * (side_p0 - side_p1);
            
            
            // the algorithm of identifying the ref_side always ensures that the
            // region with xi > x0 on ref_side will be triangular. The other side
            // will be modeled with two QUAD4 elements
            
            libMesh::Elem
            *e_p = const_cast<libMesh::Elem*>(&e),
            *e1  = libMesh::Elem::build( libMesh::TRI3, e_p).release(),
            *e2  = libMesh::Elem::build(libMesh::QUAD4, e_p).release(),
            *e3  = libMesh::Elem::build(libMesh::QUAD4, e_p).release();
            
            _new_elems.push_back(e1);
            _new_elems.push_back(e2);
            _new_elems.push_back(e3);
            
            // create the elements. We create the elements such that:
            //     node(0) = new node on ref_side
            //     node(1) = elem node (ref_side+1)%n_sides
            //     node(2) = new node on (ref_side+1)%n_sides
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident with ref_side of original element
            //     side(1) = coincident side with ref_side+1 of original element
            //     side(2) = coincident side with level-set interface
            
            e1->set_node(0) = _new_nodes[0];
            e1->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+1)%n_nodes));
            e1->set_node(2) = _new_nodes[1];
            
            _elem_sides_on_interface[e1] = 2;
            
            // create a second element and set its nodes so that:
            //     node(0) = new node on ref_side
            //     node(1) = new node on (ref_side+2)%n_nodes
            //     node(2) = elem node (ref_side+3)%n_nodes
            //     node(3) = elem node ref_side
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = new edge connecting ref_side to new node on ref_side+2
            //     side(1) = coincident side with (ref_side+2)%n_sides of original element
            //     side(2) = coincident side with (ref_side+3)%n_sides of original element
            //     side(3) = coincident with ref_side of original element
            
            e2->set_node(0) = _new_nodes[0];
            e2->set_node(1) = _new_nodes[2];
            e2->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+3)%n_nodes));
            e2->set_node(3) = const_cast<libMesh::Node*>(e.node_ptr(ref_side));
            
            // this element does not have a side on the interface
            _elem_sides_on_interface[e2] = -1;
            
            // create a third element and set its nodes so that:
            //     node(0) = new node on ref_side+1
            //     node(1) = elem node (ref_side+2)%n_nodes
            //     node(2) = new node on (ref_side+2)%n_nodes
            //     node(3) = new node on ref_side
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident side with (ref_side+1)%n_sides of original element
            //     side(1) = coincident side with (ref_side+2)%n_sides of original element
            //     side(2) = new edge connecting ref_side to new node on (ref_side+2)%n_sides
            //     side(3) = coincident side with level-set interface
            
            e3->set_node(0) = _new_nodes[1];
            e3->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+2)%n_nodes));
            e3->set_node(2) = _new_nodes[2];
            e3->set_node(3) = _new_nodes[0];
            
            _elem_sides_on_interface[e3] = 3;
            
            if (node0_positive) {
                
                _positive_phi_elems.push_back(e2);
                _positive_phi_elems.push_back(e3);
                _negative_phi_elems.push_back(e1);
            }
            else {
                
                _positive_phi_elems.push_back(e1);
                _negative_phi_elems.push_back(e2);
                _negative_phi_elems.push_back(e3);
            }
        }
            break;
            
            
            
        case MAST::NODE_AND_EDGE: {
            
            libMesh::Elem
            *e_p = const_cast<libMesh::Elem*>(&e),
            *e1  = libMesh::Elem::build(libMesh::TRI3,  e_p).release(),
            *e2  = libMesh::Elem::build(libMesh::QUAD4, e_p).release();
            
            _new_elems.push_back(e1);
            _new_elems.push_back(e2);
            
            if (ref_side-ref_node == 1) {
                
                // create the elements. We create the elements such that:
                //     node(0) = elem node (ref_node)
                //     node(1) = elem node (ref_node+1)%n_nodes
                //     node(2) = new node on ref_side%n_nodes
                // this way, the element has the same orientation as the original
                // element and the sides have the following association:
                //     side(0) = coincident with ref_node of original element
                //     side(1) = coincident side with ref_node+1 of original element
                //     side(2) = coincident side with level-set interface
                
                e1->set_node(0) = const_cast<libMesh::Node*>(e.node_ptr((ref_node)%n_nodes));
                e1->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+1)%n_nodes));
                e1->set_node(2) = _new_nodes[0];
                
                _elem_sides_on_interface[e1] = 2;
                
                // create a second element and set its nodes so that:
                //     node(0) = new node on ref_side%n_nodes
                //     node(1) = elem node (ref_side+1)%n_nodes
                //     node(2) = elem node (ref_side+2)%n_nodes
                //     node(3) = elem node (ref_node)
                // this way, the element has the same orientation as the original
                // element and the sides have the following association:
                //     side(0) = coincident with ref_side of original element
                //     side(1) = coincident side with (ref_side+1)%n_sides of original element
                //     side(2) = coincident side with (ref_side+2)%n_sides of original element
                //     side(3) = coincident side with level-set interface
                
                e2->set_node(0) = _new_nodes[0];
                e2->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+1)%n_nodes));
                e2->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_side+2)%n_nodes));
                e2->set_node(3) = const_cast<libMesh::Node*>(e.node_ptr((ref_node)%n_nodes));
                
                _elem_sides_on_interface[e2] = 3;
                
                
                // node0 is the first node of ref_side
                if (node0_positive) {
                    
                    _positive_phi_elems.push_back(e1);
                    _negative_phi_elems.push_back(e2);
                }
                else {
                    
                    _positive_phi_elems.push_back(e2);
                    _negative_phi_elems.push_back(e1);
                }
            }
            else if (ref_side-ref_node==2) {
                
                // create the elements. We create the elements such that:
                //     node(0) = elem node (ref_node)
                //     node(1) = new node on ref_side%n_nodes
                //     node(2) = elem node (ref_node+3)%n_nodes
                // this way, the element has the same orientation as the original
                // element and the sides have the following association:
                //     side(0) = coincident side with level-set interface
                //     side(1) = coincident side with ref_side%n_nodes of original element
                //     side(2) = coincident side with (ref_side+1)%n_nodes of original element
                
                e1->set_node(0) = const_cast<libMesh::Node*>(e.node_ptr((ref_node)%n_nodes));
                e1->set_node(1) = _new_nodes[0];
                e1->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+3)%n_nodes));
                
                _elem_sides_on_interface[e1] = 0;
                
                // create a second element and set its nodes so that:
                //     node(0) = elem node (ref_node)
                //     node(1) = elem node (ref_node+1)%n_nodes
                //     node(2) = elem node (ref_node+2)%n_nodes
                //     node(3) = new node on ref_side%n_nodes
                // this way, the element has the same orientation as the original
                // element and the sides have the following association:
                //     side(0) = coincident with ref_node of original element
                //     side(1) = coincident side with (ref_node+1)%n_sides of original element
                //     side(2) = coincident side with (ref_node+2)%n_sides of original element
                //     side(3) = coincident side with level-set interface
                
                e2->set_node(0) = const_cast<libMesh::Node*>(e.node_ptr((ref_node)%n_nodes));
                e2->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+1)%n_nodes));
                e2->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+2)%n_nodes));
                e2->set_node(3) = _new_nodes[0];
                
                _elem_sides_on_interface[e2] = 3;
                
                
                if (node0_positive) {
                    
                    _positive_phi_elems.push_back(e2);
                    _negative_phi_elems.push_back(e1);
                }
                else {
                    
                    _positive_phi_elems.push_back(e1);
                    _negative_phi_elems.push_back(e2);
                }
            }
            else
                libmesh_error(); // shoudl not get here
                
        }
            break;

            
        case MAST::OPPOSITE_NODES: {
            
            libMesh::Elem
            *e_p = const_cast<libMesh::Elem*>(&e),
            *e1  = libMesh::Elem::build(libMesh::TRI3, e_p).release(),
            *e2  = libMesh::Elem::build(libMesh::TRI3, e_p).release();
            
            _new_elems.push_back(e1);
            _new_elems.push_back(e2);
            
            // create the elements. We create the elements such that:
            //     node(0) = elem node (ref_node)
            //     node(1) = elem node (ref_node+1)%n_nodes
            //     node(2) = elem node (ref_node+2)%n_nodes
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident with ref_node of original element
            //     side(1) = coincident side with ref_node+1 of original element
            //     side(2) = coincident side with level-set interface
            
            e1->set_node(0) = const_cast<libMesh::Node*>(e.node_ptr((ref_node)%n_nodes));
            e1->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+1)%n_nodes));
            e1->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+2)%n_nodes));
            
            _elem_sides_on_interface[e1] = 2;
            
            // create a second element and set its nodes so that:
            //     node(0) = elem node on (ref_node+2)%n_nodes
            //     node(1) = elem node on (ref_node+3)%n_nodes
            //     node(2) = elem node on (ref_node+4)%n_nodes
            // this way, the element has the same orientation as the original
            // element and the sides have the following association:
            //     side(0) = coincident with (ref_node+2)%n_nodes of original element
            //     side(1) = coincident with (ref_node+3)%n_nodes of original element
            //     side(2) = coincident side with level-set interface
            
            e2->set_node(0) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+2)%n_nodes));
            e2->set_node(1) = const_cast<libMesh::Node*>(e.node_ptr((ref_node+3)%n_nodes));
            e2->set_node(2) = const_cast<libMesh::Node*>(e.node_ptr(ref_node));
            
            _elem_sides_on_interface[e2] = 2;
            
            if (node0_positive) {
                
                _positive_phi_elems.push_back(e1);
                _negative_phi_elems.push_back(e2);
            }
            else {
                
                _positive_phi_elems.push_back(e2);
                _negative_phi_elems.push_back(e1);
            }
        }
            break;
            
            
        default:
            // should not get here
            libmesh_error();
    }
}



Real
MAST::LevelSetIntersection::_find_intersection_on_straight_edge
(const libMesh::Point& p0,
 const libMesh::Point& p1,
 const MAST::FieldFunction<Real>& phi,
 const Real t) {

    // this uses the bisection search to find the location of intersection
    libMesh::Point
    pt,
    pt0 = p0,
    pt1 = p1;

    Real
    xi  = 0.,
    v   = 1.e10, // arbitrarily high value to get the algorithm going
    v0  = 0.,
    v1  = 0.;

    phi(pt0, t, v0);
    phi(pt1, t, v1);
    
    unsigned int
    n_iters = 0;
 
    //
    //       a0 + a1 xi = 0     (a0 = v0, a1 = (v1-v0))
    // or,   xi = -a0/a1
    //
    
    while (fabs(v) > _tol &&
           n_iters < _max_iters) {
        
        xi  = -v0 / (v1-v0);
        pt  = pt0 + (pt1 - pt0)*xi;
        
        phi(pt, t, v);
	    
        if (v*v1 < 0.) {
            
            v0  = v;
            pt0 = pt;
        }
        else {
            
            v1  = v;
            pt1 = pt;
        }
        
        n_iters++;
    }
    
	// now find the xi location based on the distance of the new
	// point from the old points
	pt  -= p0;
	xi   = pt.norm();
	pt   = p1-p0;
	xi  /= pt.norm();
	
    return xi;
}


const libMesh::Point&
MAST::LevelSetIntersection::get_nondimensional_coordinate_for_node
(const libMesh::Node& n) const {
    
    libmesh_assert(_initialized);
    
    std::map<const libMesh::Node*, libMesh::Point>::const_iterator
    it = _node_local_coords.find(&n);
    
    libmesh_assert(it != _node_local_coords.end());
    
    return it->second;
}

