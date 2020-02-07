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

#ifndef __mast__sub_elem_mesh_refinement_h__
#define __mast__sub_elem_mesh_refinement_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/system.h"


namespace MAST {
    
    // Forward declerations
    template <typename ValType> class FieldFunction;
    class LevelSetIntersection;
    class SubElemNodeMap;
    
    class SubElemMeshRefinement:
    public libMesh::System::Constraint {
      
    public:
        
        SubElemMeshRefinement(libMesh::MeshBase& mesh,
                              libMesh::System& sys);
        
        virtual ~SubElemMeshRefinement();
        
        bool initialized() { return _initialized; }
        
        bool process_mesh(const MAST::FieldFunction<Real>& phi,
                          bool strong_discontinuity,
                          Real time,
                          unsigned int negative_level_set_subdomain_offset,
                          unsigned int inactive_subdomain_offset,
                          unsigned int level_set_boundary_id);

        bool clear_mesh();

        /*!
         *   provides implementation of the libMesh::System::Constraint::constrain()
         *   virtual method
         */
        virtual void
        constrain ();

    protected:

        void _process_sub_elements(bool strong_discontinuity,
                                   unsigned int negative_level_set_subdomain_offset,
                                   unsigned int level_set_boundary_id,
                                   libMesh::Elem& e,
                                   MAST::LevelSetIntersection& intersect,
                                   bool positive_phi,
                                   const std::vector<const libMesh::Elem*>& elems);
        
        void _process_negative_element(unsigned int negative_level_set_subdomain_offset,
                                       unsigned level_set_boundary_id,
                                       libMesh::Elem& e,
                                       MAST::LevelSetIntersection& intersect);
        
        /*!
         *  \returns a node between the bounding nodes at the specified
         *  location. If a node already exists between these bounding nodes the
         *  that node is returned. Else, a new node is created, added to the
         *  mesh and returned.
         */
        libMesh::Node*
        _add_node(const libMesh::Point& p,
                  bool strong_disontinuity,
                  bool positive_phi,
                  unsigned int processor_id,
                  const std::pair<const libMesh::Node*, const libMesh::Node*>& bounding_nodes);
        
        libMesh::Elem* _add_elem();
        
        bool _initialized;

        bool _strong_discontinuity;

        unsigned int _negative_level_set_subdomain_offset;
        
        unsigned int _inactive_subdomain_offset;
        
        unsigned int _level_set_boundary_id;
        
        libMesh::MeshBase&          _mesh;

        libMesh::System&            _system;

        // map for storing new nodes
        MAST::SubElemNodeMap       *_node_map;
        
        std::set<unsigned int>      _negative_level_set_ids;
        std::vector<libMesh::Node*> _new_nodes;
        std::vector<libMesh::Elem*> _new_elems;
        std::vector<std::pair<libMesh::Elem*, unsigned int>> _old_elems;
        std::set<std::pair<const libMesh::Node*, std::pair<const libMesh::Node*, const libMesh::Node*>>>    _hanging_node;
    };
}


#endif // __mast__sub_elem_mesh_refinement_h__
