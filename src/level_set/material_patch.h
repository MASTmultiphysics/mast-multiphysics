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

#ifndef __mast_material_patch_h__
#define __mast_material_patch_h__

// C++ includes
#include <vector>

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/elem.h"


namespace MAST {

    // Forward declerations
    template <typename ValType> class FieldFunction;
    
    /*!
     *   A patch is defines as the set of elements sharing a node.
     *   This class looks through the connectivity of the node with negative
     *   level set function value and identifies if material domains are
     *   connected to it. More importantly, it identifies if the material
     *   domains are continuous across the patch. If not, the number of
     *   noncontinuous domains are identified along with the nodes that
     *   belong to each level.
     */
    
    class MaterialPatch {
        
    public:
        MaterialPatch();
        
        virtual ~MaterialPatch();

        /*!
         *   initialize the patch around \p node of \p elem.
         */
        void init(const libMesh::Elem& elem,
                  const libMesh::Node& node,
                  const MAST::FieldFunction<Real>& phi,
                  const Real t);
        
        void clear();

        const std::set<const libMesh::Elem*>&
        get_elems_to_factor() const {
            
            libmesh_assert(_initialized);
            
            return _elems_to_factor;
        }
        
    protected:

        
        bool _quad4_material_levels(const libMesh::Elem& elem,
                                    const libMesh::Node& node,
                                    const std::set<const libMesh::Elem*>& elem_neighbors,
                                    const MAST::FieldFunction<Real>& phi,
                                    const Real t);

        bool _initialized;
        
        std::set<const libMesh::Elem*> _elems_to_factor;
    };
    
}

#endif //__mast_material_patch_h__
