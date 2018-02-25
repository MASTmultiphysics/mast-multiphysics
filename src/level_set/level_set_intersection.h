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
        
        LevelSetIntersection(unsigned int max_elem_id,
                             unsigned int max_node_id);
        
        virtual ~LevelSetIntersection();

        void init(const MAST::FieldFunction<Real>& phi,
                  const libMesh::Elem& e,
                  const Real t);

        /*!
         *   clears the data structures
         */
        void clear();
        
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
         *  or edge) that is on the positive side of the level set function.
         */
        bool if_elem_has_positive_phi_region() const;
        
        const std::vector<const libMesh::Elem*>&
        get_sub_elems_positive_phi() const;
        
        const std::vector<const libMesh::Elem*>&
        get_sub_elems_negative_phi() const;

        /*!
         *   @returns the id of side that is on the interface. In case the
         *   element does not have a side on the interface, then a negative
         *   value is returned. 
         */
        unsigned int
        get_side_on_interface(const libMesh::Elem& e) const;

        bool
        has_side_on_interface(const libMesh::Elem& e) const;

        const libMesh::Point&
        get_nondimensional_coordinate_for_node(const libMesh::Node& n) const;
        
    protected:
        
        
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
        const unsigned int                           _max_mesh_elem_id;
        const unsigned int                           _max_mesh_node_id;
        const unsigned int                           _max_elem_divs;
        
        bool                                         _initialized;

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
    };
    
}




#endif // __mast_level_set_intersection_h__

