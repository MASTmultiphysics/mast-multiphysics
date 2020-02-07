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
#include "examples/base/augment_ghost_elem_send_list.h"


MAST::AugmentGhostElementSendListObj::
AugmentGhostElementSendListObj(libMesh::System& sys):
_sys(sys)
{ }
        
        

void
MAST::AugmentGhostElementSendListObj::
augment_send_list (std::vector<libMesh::dof_id_type> & send_list) {
            
    libMesh::MeshBase& mesh = _sys.get_mesh();
    libMesh::DofMap&   dofs = _sys.get_dof_map();
    
    libMesh::MeshBase::const_element_iterator
    e_it   = mesh.ghost_elements_begin(),
    e_end  = mesh.ghost_elements_end();
    
    std::vector<libMesh::dof_id_type> dof_indices;
    
    for ( ; e_it != e_end; e_it++) {
        
        dofs.dof_indices(*e_it, dof_indices);
        
        for (unsigned int i=0; i <dof_indices.size(); i++)
            if (dof_indices[i] <  dofs.first_dof() ||
                dof_indices[i] >= dofs.end_dof())
                send_list.push_back(dof_indices[i]);
    }
}
        
