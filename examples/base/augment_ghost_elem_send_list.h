/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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

#ifndef __mast__augment_ghost_elem_send_list_h__
#define __mast__augment_ghost_elem_send_list_h__

// C++ includes
#include <vector>

// libMesh includes
#include "libmesh/system.h"
#include "libmesh/dof_map.h"


namespace MAST {
    
    class AugmentGhostElementSendListObj:
    public libMesh::DofMap::AugmentSendList {
        
    public:
        
        AugmentGhostElementSendListObj(libMesh::System& sys);
        
        
        virtual ~AugmentGhostElementSendListObj () {}
        
        /**
         * User-defined function to augment the send list.
         */
        virtual void augment_send_list (std::vector<libMesh::dof_id_type> & send_list);
        
    protected:
        
        libMesh::System&          _sys;
    };
}

#endif // __mast__augment_ghost_elem_send_list_h__
