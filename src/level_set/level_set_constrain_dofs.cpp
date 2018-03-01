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

// MAST includes
#include "level_set/level_set_constrain_dofs.h"
#include "level_set/level_set_intersection.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/mesh_field_function.h"


// libMesh includes
#include "libmesh/dof_map.h"



MAST::LevelSetConstrainDofs::
LevelSetConstrainDofs(MAST::SystemInitialization& sys,
                      MAST::FieldFunction<Real>& level_set):
_sys           (sys),
_level_set     (level_set),
_intersection  (nullptr) {
    
    _intersection = new MAST::LevelSetIntersection(sys.system().get_mesh().max_elem_id(),
                                                   sys.system().get_mesh().max_node_id());
}



MAST::LevelSetConstrainDofs::~LevelSetConstrainDofs() {
    
    delete _intersection;
}


MAST::LevelSetIntersection&
MAST::LevelSetConstrainDofs::get_intersection() {
    
    return *_intersection;
}



void
MAST::LevelSetConstrainDofs::constrain() {
    
    MAST::NonlinearSystem& nonlin_sys = _sys.system();
    
    libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    std::vector<libMesh::dof_id_type>
    dof_indices;
    std::set<libMesh::dof_id_type>
    all_dof_indices,
    connected_dof_indices;
    
    // our intent is to constrain only those dofs that belong to the
    // unintersected elements on the negative phi side of level set, AND
    // if they do not belong to any elements intersected by the level set.
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(_level_set, *elem, nonlin_sys.time);
        
        // if the element is entirely on the negative side of the level set,
        // we will constrain all dofs of the element to zero
        
        dof_indices.clear();
        dof_map.dof_indices(elem, dof_indices);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            all_dof_indices.insert(dof_indices[i]);
        
        if (_intersection->if_elem_has_positive_phi_region()) {
            // This identifies if the element is has some positive phi
            // region. If it does, then it will provide a contribution to
            // the dofs connected to it.
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                connected_dof_indices.insert(dof_indices[i]);
        }
        
        _intersection->clear();
    }
    
    // create a set so that we only deal with unique set of ids.
    dof_indices.clear();
    dof_indices.reserve(all_dof_indices.size() - connected_dof_indices.size());
    
    std::set_difference(all_dof_indices.begin(),
                        all_dof_indices.end(),
                        connected_dof_indices.begin(),
                        connected_dof_indices.end(),
                        std::inserter(dof_indices, dof_indices.begin()));
    
    all_dof_indices.clear();
    connected_dof_indices.clear();
    
    // now, constrain everythign in the set
    std::vector<libMesh::dof_id_type>::const_iterator
    dof_it  = dof_indices.begin(),
    dof_end = dof_indices.end();
    
    for ( ; dof_it != dof_end; dof_it++) {
        
        // if the dof is already Dirichlet constrained, then we do not
        // add another constraint on it
        if (!dof_map.is_constrained_dof(*dof_it)) {
            
            libMesh::DofConstraintRow c_row;
            dof_map.add_constraint_row(*dof_it, c_row, true);
        }
    }
}



