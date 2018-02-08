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
                  const Real t);

        /*!
         *   clears the data structures
         */
        void clear();
        
        /*!
         *   @returns \p true if the mode is something other than
         *   NO_INTERSECTION.
         */
        bool if_intersection() const;
        
        /*!
         *   @returns \p true if the mode is ADJACENT_EDGES or OPPOSITE_EDGES.
         */
        bool if_intersection_through_elem() const;

        /*!
         *   @returns true if the element is entirely on the positive side
         *   of the level set without any intersection. This will return
         *   \p true only if the mode is NO_INTERSECTION and all nodal phi
         *   values are positive.
         */
        bool if_elem_on_positive_phi() const;

        const std::vector<const libMesh::Elem*>&
        get_sub_elems_positive_phi() const;
        
        const std::vector<const libMesh::Elem*>&
        get_sub_elems_negative_phi() const;

        /*!
         *   @returns the id of side that is on the interface. In case the
         *   element does not have a side on the interface, then a negative
         *   value is returned. 
         */
        int
        get_side_on_interface(const libMesh::Elem& e) const;
        
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
        
        bool                                         _initialized;

        bool                                         _if_elem_on_positive_phi;
        
        MAST::LevelSet2DIntersectionMode             _mode;
        
        std::vector<const libMesh::Elem*>            _positive_phi_elems;
        
        std::vector<const libMesh::Elem*>            _negative_phi_elems;
        
        std::map<const libMesh::Elem*, int>          _elem_sides_on_interface;
        
        std::vector<libMesh::Node*>                  _new_nodes;
        std::vector<libMesh::Elem*>                  _new_elems;
        std::map<const libMesh::Node*, libMesh::Point> _node_local_coords;
    };
    
}




#endif // __mast_level_set_intersection_h__

