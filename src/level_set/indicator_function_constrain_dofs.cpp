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
#include "level_set/indicator_function_constrain_dofs.h"
#include "level_set/level_set_intersection.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/mesh_field_function.h"


// libMesh includes
#include "libmesh/dof_map.h"



MAST::IndicatorFunctionConstrainDofs::
IndicatorFunctionConstrainDofs(MAST::SystemInitialization&        sys,
                               MAST::FieldFunction<Real>&         level_set,
                               MAST::FieldFunction<RealVectorX>&  indicator):
MAST::LevelSetConstrainDofs(sys, level_set),
_indicator     (indicator) {
    
}



MAST::IndicatorFunctionConstrainDofs::~IndicatorFunctionConstrainDofs() {
    
}



void
MAST::IndicatorFunctionConstrainDofs::constrain() {
    
    // first constrain dofs based on level-set
    MAST::LevelSetConstrainDofs::constrain();
    
    // now, we add additional constraints based on indicator function
    MAST::NonlinearSystem& nonlin_sys = _sys.system();
    
    libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    std::vector<libMesh::dof_id_type>
    dof_indices;
    std::set<libMesh::dof_id_type>
    constrained_dof_indices;

    bool
    constrain_elem = true;

    RealVectorX
    vec = RealVectorX::Zero(1);
    
    Real
    tol = 1.e-6;
    
    // our intent is to constrain only those dofs that belong to elements
    // where all nodes have zero indicator function value
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        // if the element is entirely on the negative side of the level set,
        // we will constrain all dofs of the element to zero
        //
        // if any of the nodes have a positive value, we will not constrain
        // the element
        constrain_elem = true;
        for (unsigned int i=0; i<elem->n_nodes(); i++) {
            
            _indicator(elem->node_ref(i), nonlin_sys.time, vec);
            if (vec(0) > tol) {
                constrain_elem = false;
                break;
            }
        }
        
        if (constrain_elem) {

            dof_indices.clear();
            dof_map.dof_indices(elem, dof_indices);

            for (unsigned int i=0; i<dof_indices.size(); i++)
                constrained_dof_indices.insert(dof_indices[i]);
        }
    }
    
    // create a set so that we only deal with unique set of ids.
    dof_indices.clear();
    
    // now, constrain everythign in the set
    std::set<libMesh::dof_id_type>::const_iterator
    dof_it  = constrained_dof_indices.begin(),
    dof_end = constrained_dof_indices.end();
    
    for ( ; dof_it != dof_end; dof_it++) {
        
        // if the dof is already Dirichlet constrained, then we do not
        // add another constraint on it
        if (!dof_map.is_constrained_dof(*dof_it)) {
            
            libMesh::DofConstraintRow c_row;
            dof_map.add_constraint_row(*dof_it, c_row, true);
        }
    }
}



