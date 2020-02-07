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

// C++ includes
#include <set>

// MAST includes
#include "examples/fsi/beam_flutter_solution/constrain_beam_dofs.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"



MAST::ConstrainBeamDofs::ConstrainBeamDofs(MAST::NonlinearSystem& sys):
libMesh::System::Constraint(),
_system(sys) {
    
}


MAST::ConstrainBeamDofs::~ConstrainBeamDofs() {
    
}


void
MAST::ConstrainBeamDofs::constrain() {

    // variables to constrain: u, w, tx, ty
    std::vector<unsigned int>
    c_vars;
    c_vars.push_back(0);   // u
    c_vars.push_back(2);   // w
    c_vars.push_back(3);   // tx
    c_vars.push_back(4);   // ty
    
    
    for ( unsigned int i=0; i<c_vars.size(); i++) {
        
        std::set<libMesh::dof_id_type>
        all_dofs;
        std::vector<libMesh::dof_id_type>
        dofs;
        
        
        libMesh::MeshBase
        &mesh  = _system.get_mesh();
        
        libMesh::DofMap
        &dofmap= _system.get_dof_map();
        
        // iterate over all local element and nodes and constrain them
        libMesh::MeshBase::element_iterator
        e_it    =   mesh.local_elements_begin(),
        e_end   =   mesh.local_elements_end();
        
        for (; e_it!=e_end; e_it++) {
            
            dofs.clear();
            dofmap.dof_indices(*e_it, dofs, c_vars[i]);
            all_dofs.insert(dofs.begin(), dofs.end());
        }
        
        // iterate over all local element and nodes and constrain them
        libMesh::MeshBase::node_iterator
        n_it    =   mesh.local_nodes_begin(),
        n_end   =   mesh.local_nodes_end();
        
        for (; n_it!=n_end; n_it++) {
            
            dofs.clear();
            dofmap.dof_indices(*n_it, dofs, c_vars[i]);
            all_dofs.insert(dofs.begin(), dofs.end());
        }
        
        // now that the total set of indices is available on the local
        // processor, we will add constraints for all these dofs
        std::set<libMesh::dof_id_type>::const_iterator
        dof_it  =  all_dofs.begin(),
        dof_end =  all_dofs.end();
        
        for ( ; dof_it != dof_end; dof_it++) {
            
            if(!dofmap.is_constrained_dof(*dof_it)) {
                
                libMesh::DofConstraintRow crow;
                dofmap.add_constraint_row(*dof_it, crow, 0., true);
            }
        }
    }
}


