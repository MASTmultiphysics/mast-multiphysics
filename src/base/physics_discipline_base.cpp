/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/parameter.h"
#include "boundary_condition/dirichlet_boundary_condition.h"

// libMesh includes
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/fe_interface.h"
#include "libmesh/dirichlet_boundaries.h"


void
MAST::PhysicsDisciplineBase::clear_loads() {
    _side_bc_map.clear();
    _vol_bc_map.clear();
}



void
MAST::PhysicsDisciplineBase::add_side_load(libMesh::boundary_id_type bid,
                                           MAST::BoundaryConditionBase& load) {
    // make sure that this boundary and load haven't already been applied
    std::pair<MAST::SideBCMapType::const_iterator, MAST::SideBCMapType::const_iterator> it =
    _side_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++)
        libmesh_assert(it.first->second != &load);
    
    // displacement boundary condition needs to be hadled separately
    _side_bc_map.insert(MAST::SideBCMapType::value_type(bid, &load));
}




void
MAST::PhysicsDisciplineBase::add_dirichlet_bc(libMesh::boundary_id_type bid,
                                              MAST::DirichletBoundaryCondition& load) {
    
    bool insert_success =
    _dirichlet_bc_map.insert(MAST::DirichletBCMapType::value_type(bid, &load)).second;
    
    libmesh_assert(insert_success);
}



void
MAST::PhysicsDisciplineBase::add_volume_load(libMesh::subdomain_id_type bid,
                                             MAST::BoundaryConditionBase& load) {
    std::pair<MAST::VolumeBCMapType::iterator, MAST::VolumeBCMapType::iterator> it =
    _vol_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++)
        libmesh_assert(it.first->second != &load);
    
    _vol_bc_map.insert(MAST::VolumeBCMapType::value_type(bid, &load));
}


void
MAST::PhysicsDisciplineBase::add_point_load(MAST::PointLoadCondition& load) {

    libmesh_assert(_point_loads.count(&load));
    
    _point_loads.insert(&load);
}



void
MAST::PhysicsDisciplineBase::clear_volume_load(libMesh::subdomain_id_type bid,
                                               MAST::BoundaryConditionBase& load) {
    std::pair<MAST::VolumeBCMapType::iterator, MAST::VolumeBCMapType::iterator> it =
    _vol_bc_map.equal_range(bid);
    
    for ( ; it.first != it.second; it.first++)
        if (it.first->second == &load) {
            _vol_bc_map.erase(it.first);
            return;
        }
    
    // should not get here
    libmesh_assert(false);
}



void
MAST::PhysicsDisciplineBase::
set_property_for_subdomain(const libMesh::subdomain_id_type sid,
                           const MAST::ElementPropertyCardBase& prop) {
    
    MAST::PropertyCardMapType::const_iterator elem_p_it = _element_property.find(sid);
    libmesh_assert(elem_p_it == _element_property.end());
    
    _element_property[sid] = &prop;
}



const MAST::ElementPropertyCardBase&
MAST::PhysicsDisciplineBase::get_property_card(const unsigned int i) const {
    
    MAST::PropertyCardMapType::const_iterator
    elem_p_it = _element_property.find(i);
    libmesh_assert(elem_p_it != _element_property.end());
    
    return *elem_p_it->second;
}



const MAST::ElementPropertyCardBase&
MAST::PhysicsDisciplineBase::get_property_card(const libMesh::Elem& elem) const {
    
    MAST::PropertyCardMapType::const_iterator
    elem_p_it = _element_property.find(elem.subdomain_id());
    libmesh_assert(elem_p_it != _element_property.end());
    
    return *elem_p_it->second;
}



void
MAST::PhysicsDisciplineBase::add_parameter(MAST::Parameter& f) {
    
    Real* par = f.ptr();
    // make sure it does not already exist in the map
    libmesh_assert(!_parameter_map.count(par));
    
    // now add this to the map
    bool insert_success = _parameter_map.insert
    (std::map<Real*, MAST::FunctionBase*>::value_type(par, &f)).second;
    
    libmesh_assert(insert_success);
}



const MAST::FunctionBase*
MAST::PhysicsDisciplineBase::get_parameter(Real* par) const {
    // make sure valid values are given
    libmesh_assert(par);
    
    std::map<Real*, const MAST::FunctionBase*>::const_iterator
    it = _parameter_map.find(par);
    
    // make sure it does not already exist in the map
    libmesh_assert(it != _parameter_map.end());
    
    return it->second;
}



template <typename SysType>
void
MAST::PhysicsDisciplineBase::
init_system_dirichlet_bc(SysType& sys) const {
    
    
    // iterate over all the dirichlet boundary conditions and add them
    // to the system
    MAST::DirichletBCMapType::const_iterator it = _dirichlet_bc_map.begin();
    
    for ( ; it != _dirichlet_bc_map.end(); it++)
        sys.get_dof_map().add_dirichlet_boundary(it->second->dirichlet_boundary());
}




template <typename SysType>
void
MAST::PhysicsDisciplineBase::
clear_system_dirichlet_bc(SysType& sys) const {
    
    // iterate over all the dirichlet boundary conditions and add them
    // to the system
    MAST::DirichletBCMapType::const_iterator it = _dirichlet_bc_map.begin();
    
    for ( ; it != _dirichlet_bc_map.end(); it++)
        sys.get_dof_map().remove_dirichlet_boundary(it->second->dirichlet_boundary());
}




template <>
void MAST::PhysicsDisciplineBase::
init_system_dirichlet_bc<libMesh::CondensedEigenSystem>(libMesh::CondensedEigenSystem& sys) const {
    
    // first prepare a map of boundary ids and the constrained vars on that
    // boundary
    std::map<libMesh::boundary_id_type, std::vector<unsigned int> >  constrained_vars_map;
    
    // iterate over all the dirichlet boundary conditions and add them
    // to the system
    MAST::DirichletBCMapType::const_iterator it = _dirichlet_bc_map.begin();
    
    for ( ; it != _dirichlet_bc_map.end(); it++) {
        
        libMesh::DirichletBoundary& dirichlet_b = it->second->dirichlet_boundary();
        sys.get_dof_map().add_dirichlet_boundary(dirichlet_b);
        constrained_vars_map[it->first] = dirichlet_b.variables;
    }
    
    
    //
    // now collect the ids that correspond to the specified boundary conditions
    //
    // Get a constant reference to the mesh object
    const libMesh::MeshBase& mesh = sys.get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a constant reference to the Finite Element type
    // for the first variable in the system.
    libMesh::FEType fe_type = sys.get_dof_map().variable_type(0);
    
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    
    // the constrained dofs needed for CondensedEigenSystem
    std::set<unsigned int> dof_ids;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    
    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        const libMesh::Elem* elem = *el;
        
        // boundary condition is applied only on sides with no neighbors
        // and if the side's boundary id has a boundary condition tag on it
        for (unsigned int s=0; s<elem->n_sides(); s++)
            if ((*el)->neighbor(s) == NULL &&
                mesh.boundary_info->n_boundary_ids(elem, s)) {
                
                std::vector<libMesh::boundary_id_type> bc_ids = mesh.boundary_info->boundary_ids(elem, s);
                
                for (unsigned int i_bid=0; i_bid<bc_ids.size(); i_bid++)
                    if (constrained_vars_map.count(bc_ids[i_bid])) {
                        
                        const std::vector<unsigned int>& vars = constrained_vars_map[bc_ids[i_bid]];
                        // now iterate over each constrained variable for this boundary
                        // and collect its dofs
                        for (unsigned int i_var=0; i_var<vars.size(); i_var++) {
                            
                            dof_indices.clear();
                            dof_map.dof_indices (*el, dof_indices, vars[i_var]);
                            
                            // All boundary dofs are Dirichlet dofs in this case
                            std::vector<unsigned int> side_dofs;
                            libMesh::FEInterface::dofs_on_side(*el, dim, fe_type,
                                                               s, side_dofs);
                            
                            for(unsigned int ii=0; ii<side_dofs.size(); ii++)
                                dof_ids.insert(dof_indices[side_dofs[ii]]);
                        }
                    } // end of boundary loop
            } // end of side loop
    }// end of element loop

    // also, it is likely that some of the bcs have been applied via the
    // DofMap API by specification of row constraints. In that case,
    // factor out the dofs that are not coupled to any other dofs
    libMesh::DofConstraints::const_iterator
    constraint_it  = dof_map.constraint_rows_begin(),
    constraint_end = dof_map.constraint_rows_end();
    
    for ( ; constraint_it != constraint_end; constraint_it++) {
        // if the dof constraint has only one entry, then add it to the
        // constrained set
        if (!constraint_it->second.size())
            dof_ids.insert(constraint_it->first);
    }
    
    
    // now that the dofs are available, tell the system to condense out
    // the constrained dofs
    sys.initialize_condensed_dofs(dof_ids);
}



// explicit instantiations
template void MAST::PhysicsDisciplineBase::
init_system_dirichlet_bc<libMesh::System>(libMesh::System& sys) const;

template void MAST::PhysicsDisciplineBase::
clear_system_dirichlet_bc<libMesh::System>(libMesh::System& sys) const;

