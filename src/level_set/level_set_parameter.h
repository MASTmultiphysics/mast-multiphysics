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

#ifndef __mast__level_set_parameter__
#define __mast__level_set_parameter__


// MAST includes
#include "base/parameter.h"

// libMesh includes
#include "libmesh/node.h"


namespace MAST {
    
    /*!
     *    This defines a parameter that is a level set function
     *    and stores a pointer to the node in the level-set mesh
     *    whose value is defiend by this parameter.
     */
    class LevelSetParameter:
    public MAST::Parameter {
        
    public:
        
        LevelSetParameter(const std::string& nm,
                          const Real& val,
                          const libMesh::Node* node):
        MAST::Parameter(nm, val),
        _node(node) { }
        
        virtual ~LevelSetParameter() {
            
        }
        

        const libMesh::Node* level_set_node() const {
            
            return _node;
        }
        
    protected:
        
        /*!
         *    Pointer to the level set node
         */
        const libMesh::Node* _node;
    };
}

#endif // __mast__level_set_parameter__
