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

//MAST includes
#include "elasticity/kinematic_coupling.h"
#include "elasticity/kinematic_coupling_constraint.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "mesh/mesh_coupling_base.h"

// libMesh includes
#include "libmesh/dof_map.h"


MAST::KinematicCoupling::KinematicCoupling(MAST::SystemInitialization& sys_init):
MAST::DoFCouplingBase(sys_init),
_initialized   (false) {
    
}


MAST::KinematicCoupling::~KinematicCoupling() {
    
    std::map<const libMesh::Node*, const MAST::KinematicCouplingConstraint*>::const_iterator
    it   = _slave_node_constraints.begin(),
    end  = _slave_node_constraints.end();

    for ( ; it != end; it++)
        delete it->second;
}



void
MAST::KinematicCoupling::
add_master_and_slave
(std::vector<std::pair<const libMesh::Node*, std::set<const libMesh::Node*>>>& couplings,
 bool constrain_rotations) {
    
    for (unsigned int i=0; i<couplings.size(); i++)
        this->add_master_and_slave(couplings[i].second,
                                   *couplings[i].first,
                                   constrain_rotations);
}


void
MAST::KinematicCoupling::
add_master_and_slave(std::set<const libMesh::Node*>& master,
                     const libMesh::Node&            slave,
                     bool                            constrain_rotations) {
    
    libmesh_assert(!_initialized);
    libmesh_assert(!_slave_node_constraints.count(&slave));
    
    MAST::KinematicCouplingConstraint
    *constr = new MAST::KinematicCouplingConstraint(_system_init,
                                                    slave,
                                                    master,
                                                    constrain_rotations);
    
    _slave_node_constraints[&slave] = constr;
}



void
MAST::KinematicCoupling::
get_constraint_rows
(std::vector<std::tuple<libMesh::dof_id_type, libMesh::DofConstraintRow, Real>>& constrs) {

    // iterate over all the nodes in the map and initialize their constraints
    constrs.clear();

    std::map<const libMesh::Node*, const MAST::KinematicCouplingConstraint*>::const_iterator
    it   = _slave_node_constraints.begin(),
    end  = _slave_node_constraints.end();
    
    unsigned int
    idx  = 0;
    
    for ( ; it != end; it++)
        idx += it->second->n_constrain_nodes();

    constrs.resize(idx);

    it   = _slave_node_constraints.begin();
    
    for ( ; it != end; it++) {
        
        std::vector<std::tuple<libMesh::dof_id_type, libMesh::DofConstraintRow, Real>>
        constraints;
        
        it->second->get_dof_constraint_row(constraints);
        
        for (unsigned int i=0; i<constrs.size(); i++)
            constrs.push_back(constraints[i]);
    }

    libmesh_assert_equal_to(constrs.size(), idx);
    
    // set the flag so that new slave-master sets cannot be added
    _initialized = true;
}

