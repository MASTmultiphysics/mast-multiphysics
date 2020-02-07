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

#ifndef __mast_level_set_intersection_h__
#define __mast_level_set_intersection_h__

// C++ includes
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"


namespace MAST {

    enum LevelSet2DIntersectionMode {
        
        THROUGH_NODE,           // level set passes through node
        COLINEAR_EDGE,          // level set coliniear with edge of element
        ADJACENT_EDGES,         // level set passes through two adjacent edges
        OPPOSITE_EDGES,         // level set passes through two opposite edges
        OPPOSITE_NODES,         // level set passes through diagonally opposite nodes
        NODE_AND_EDGE,          // level set passes through a node and edge
        TWO_ADJACENT_EDGES,     // level set passes through four edges
        NODE_AND_TWO_EDGES,     // level set passes through a node and two edges
        NO_INTERSECTION
    };
    

    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    
    
    class LevelSetIntersection {
        
    public:
        
        LevelSetIntersection();
        
        virtual ~LevelSetIntersection();

        void init(const MAST::FieldFunction<Real>& phi,
                  const libMesh::Elem& e,
                  const Real t,
                  unsigned int max_elem_id,
                  unsigned int max_node_id);

        /*!
         *   clears the data structures
         */
        void clear();
        
        /*!
         *   @return a reference to the element on which the intersection is
         *   defined
         */
        const libMesh::Elem& elem() const;
        
        /*!
         *   @returns mode of intersection
         */
        MAST::LevelSet2DIntersectionMode
        get_intersection_mode() const;

        /*!
         *   @returns the node number on the element if the mode is THROUDH_NODE,
         *   otherwise throws an error.
         */
        unsigned int node_on_boundary() const;

        /*!
         *   @returns the edge number on the element if the mode is COLINEAR_EDGE,
         *   otherwise throws an error.
         */
        unsigned int edge_on_boundary() const;


        /*!
         *  @returns the area/volume fraction of element on positive phi of
         *  the element
         */
        Real get_positive_phi_volume_fraction() const;
        
        /*!
         *   @returns value of phi on element node. This can only be the
         *   node in the
         */
        Real get_node_phi_value(const libMesh::Node* n) const;
        

        /*!
         * The case of two adjacent edges results in a new node on an edge that is not coincident
         * with the level set interface. This will end up as a hanging node if added to the mesh.
         * @returns true if this node \p n is the hanging node.
         */
        bool if_hanging_node(const libMesh::Node* n) const;
        
        /*!
         *   @returns true if the intersection is through the element, or
         *   has colinear edge.
         */
        bool if_elem_has_boundary() const;
        
        /*!
         *   @returns \p true if the intersection passes through the interior
         *   of the element. Note that COLINEAR_EDGE or THROUGH_NODE do not
         *   qualify for this. 
         */
        bool if_intersection_through_elem() const;

        /*!
         *   @returns \p true if the element is entirely on the positive side
         *   of the level set without any intersection. This will return
         *   \p true only if the mode is NO_INTERSECTION and all nodal phi
         *   values are positive.
         */
        bool if_elem_on_positive_phi() const;

        /*!
         *   @returns \p true if the element is entirely on the negative side
         *   of the level set without any intersection. This will return
         *   \p true only if the mode is NO_INTERSECTION and all nodal phi
         *   values are negative.
         */
        bool if_elem_on_negative_phi() const;

        /*!
         *  @returns \p true if there is any portion of the element (interior
         *  edge, or node) that is on the positive side of the level set function.
         */
        bool if_elem_has_positive_phi_region() const;

        
        /*!
         *  @returns \p true if there is any portion of the element (interior
         *  edge, or node) that is on the negative side of the level set function.
         */
        bool if_elem_has_negative_phi_region() const;

        const std::vector<const libMesh::Elem*>&
        get_sub_elems_positive_phi() const;
        
        const std::vector<const libMesh::Elem*>&
        get_sub_elems_negative_phi() const;

        void
        get_corner_nodes_on_negative_phi(std::set<const libMesh::Node*>& nodes) const;
        
        /*!
         *   @returns the id of side that is on the interface. In case the
         *   element does not have a side on the interface, then an error is
         *   thrown. 
         */
        unsigned int
        get_side_on_interface(const libMesh::Elem& e) const;

        bool
        has_side_on_interface(const libMesh::Elem& e) const;

        const libMesh::Point&
        get_nondimensional_coordinate_for_node(const libMesh::Node& n) const;


        /*!
         *   identifies if the node from the subelements is a new node or
         *   an existing node from the parent element.
         */
        bool if_node_is_new(const libMesh::Node& node) const;

        
        /*!
         *   identifies if the new node is on an edge along the level-set method
         *   in the interior of the element (as opposed to on the edges of the
         *   original element).
         */
        bool if_interior_node(const libMesh::Node& node) const;

        
        /*!
         *    for new nodes required to create the subelements this method
         *    returns the nodes on an edge that bound the given node. An
         *    error is thrown if the node is not a
         */
        std::pair<const libMesh::Node*, const libMesh::Node*>
        get_bounding_nodes_for_node(const libMesh::Node& node) const;
        
        
        /*!
         *    identifies the sides of the element that are completely on the material side without
         *    any intersection on them.
         */
        void
        get_material_sides_without_intersection(std::set<unsigned int>& sides) const;

        
        void get_nearest_intersection_point(const libMesh::Point& p,
                                            libMesh::Point& pt);
        
    protected:
        
        /*!
         *   creates a first order element from the given high-order element.
         *   For a QUAD9 a QUAD4 is obtained by only using the corner nodes.
         *   Note that this does not create any new nodes. The element can be
         *   deleted after use.
         */
        std::unique_ptr<libMesh::Elem>
        _first_order_elem(const libMesh::Elem& e);
        
        /*!
         *   initializes on a reference element that is a first-order
         *   counterpart of the given high-order element. For two-dimensional
         *   elements this is a QUAD4.
         */
        void _init_on_first_order_ref_elem(const MAST::FieldFunction<Real>& phi,
                                           const libMesh::Elem& e,
                                           const Real t);

        
        void _add_node_local_coords
        ( const libMesh::Elem& e,
         std::vector<std::pair<libMesh::Point, libMesh::Point> >& side_nondim_points,
         std::map<const libMesh::Node*, libMesh::Point>& node_coord_map);
        
        void
        _find_quad4_intersections(const MAST::FieldFunction<Real>& phi,
                                  const libMesh::Elem& e,
                                  const Real t,
                                  const std::map<const libMesh::Node*, std::pair<Real, bool> >&
                                  node_phi_vals);

        /*!
         *   @returns the value along a straight edge where the level-set
         *   function zero.
         */
        Real
        _find_intersection_on_straight_edge(const libMesh::Point& p0,
                                            const libMesh::Point& p1,
                                            const MAST::FieldFunction<Real>& phi,
                                            const Real t);
        
        Real                                         _tol;
        
        unsigned int                                 _max_iters;
        unsigned int                                 _max_mesh_elem_id;
        unsigned int                                 _max_mesh_node_id;
        const unsigned int                           _max_elem_divs;
        const libMesh::Elem*                         _elem;
        bool                                         _initialized;

        const MAST::FieldFunction<Real>*             _phi;
        
        /*!
         *   \p true if element is completely on the positive side of level set
         *   with no intersection
         */
        bool                                         _if_elem_on_positive_phi;
        
        /*!
         *   \p true if element is completely on the negative side of level set
         *   with no intersection
         */
        bool                                         _if_elem_on_negative_phi;
        
        MAST::LevelSet2DIntersectionMode             _mode;
        
        unsigned int                                 _node_num_on_boundary;
        
        unsigned int                                 _edge_num_on_boundary;
        
        std::vector<const libMesh::Elem*>            _positive_phi_elems;
        
        std::vector<const libMesh::Elem*>            _negative_phi_elems;
        
        std::map<const libMesh::Elem*, int>          _elem_sides_on_interface;
        
        std::vector<libMesh::Node*>                  _new_nodes;
        std::vector<libMesh::Elem*>                  _new_elems;
        std::map<const libMesh::Node*, libMesh::Point> _node_local_coords;
        std::map<const libMesh::Node*, std::pair<Real, bool> > _node_phi_vals;
        std::set<const libMesh::Node*>               _interior_nodes;
        std::map<const libMesh::Node*, std::pair<const libMesh::Node*, const libMesh::Node*>> _bounding_nodes;
        std::set<const libMesh::Node*>               _hanging_node;
    };
    
}




#endif // __mast_level_set_intersection_h__

