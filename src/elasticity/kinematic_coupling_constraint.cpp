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
#include "elasticity/kinematic_coupling_constraint.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"


MAST::KinematicCouplingConstraint::
KinematicCouplingConstraint(MAST::SystemInitialization& sys_init,
                            const libMesh::Node& slave_node,
                            const std::set<const libMesh::Node*>& master_nodes,
                            bool constrain_rotations):
_sys_init              (sys_init),
_slave                 (&slave_node),
_masters               (master_nodes),
_constrain_rotations   (constrain_rotations) {
    
}


MAST::KinematicCouplingConstraint::~KinematicCouplingConstraint() {
    
}


void
MAST::KinematicCouplingConstraint::
get_dof_constraint_row(std::vector<std::tuple
                       <libMesh::dof_id_type, libMesh::DofConstraintRow,Real>>& constraints) const {
    
    
    constraints.clear();
    
    libMesh::DofMap& dof_map = _sys_init.system().get_dof_map();

    std::vector<libMesh::dof_id_type>
    slave_dofs;

    // master nodes, their dofs and weights
    std::map<const libMesh::Node*, std::pair<std::vector<libMesh::dof_id_type>, Real>>
    master_node_data; // master dofs of all nodes

    
    // get dofs for all the nodes
    dof_map.dof_indices(_slave, slave_dofs);

    std::set<const libMesh::Node*>::const_iterator
    n_it   = _masters.begin(),
    n_end  = _masters.end();
    
    libMesh::Point dx;
    
    std::map<const libMesh::Node*, Real> weights;
    
    Real
    sum = 0.,
    val = 0.;
    
    // compute the position vectors for nodal weight
    for ( ; n_it != n_end; n_it++ ) {

        const libMesh::Node* nd = *n_it;

        dx   = *_slave - *nd;
        
        // exponential function, so that closer nodes are more influential
        // than distant nodes.
        val  = exp(-dx.norm());
        sum += val;
        
        weights[nd] = val;
    }
    
    // now compute the constraint data
    n_it   = _masters.begin();
    
    for ( ; n_it != n_end; n_it++ ) {
        
        const libMesh::Node* nd = *n_it;

        std::pair<std::vector<libMesh::dof_id_type>, Real> dofs_wt_pair;
        dof_map.dof_indices(nd, dofs_wt_pair.first);
        // divide by sum so that all weights sum to a value of 1.
        dofs_wt_pair.second = weights[nd]/sum;
        
        master_node_data[nd] = dofs_wt_pair;
    }

    // now initialize the vector with constraint data. This assumes that dof_ids
    // for all nodes have the same sequence: ux, uy, uz, tx, ty, tz
    unsigned int
    dofs_to_constrain = _constrain_rotations?6:3,
    idx               = 0;

    constraints.resize(dofs_to_constrain);
    
    // the translation and rotation dofs both share a common code, except
    // translation constraints require a contribution from the rotational dofs
    for (unsigned int i=0; i<dofs_to_constrain; i++) {
        
        // dof id of ith variable for slave node
        std::get<0>(constraints[i]) = slave_dofs[i];
        
        // rhs of the constraint is zero
        std::get<2>(constraints[i]) = 0.;
        
        libMesh::DofConstraintRow&
        c_row = std::get<1>(constraints[i]);
        
        std::map<const libMesh::Node*,
        std::pair<std::vector<libMesh::dof_id_type>, Real>>::const_iterator
        master_nd_it   = master_node_data.begin(),
        master_nd_end  = master_node_data.end();
        
        for ( ; master_nd_it != master_nd_end; master_nd_it++) {
            
            dx = *_slave - *master_nd_it->first;
            
            const std::vector<libMesh::dof_id_type>
            &master_dofs = master_nd_it->second.first;
            
            // dof id of the i^th var on the master node
            idx = master_dofs[i];
            val = master_nd_it->second.second; // weight for this node
            
            // add constraints with rotational dofs of all nodes
            // -1 is used since the constraint reads
            //     slave dof = avg of master dofs
            // or, slave dof - avg of master dofs = 0
            c_row[idx] = -val * -1.;
            
            // for translation dofs a contribution is needed from the
            // rotation of the master nodes
            if (i<3) {
                
                switch (i) {
                    case 0: {
                        // index of the rotational dof of the master node
                        c_row[master_dofs[4]]  = -dx(2) * val; //  theta_y * dz
                        c_row[master_dofs[5]]  =  dx(1) * val; // -theta_z * dy
                    }
                        break;
                        
                    case 1: {
                        // index of the rotational dof of the master node
                        c_row[master_dofs[3]]  = -dx(2) * val; //  theta_x * dz
                        c_row[master_dofs[5]]  =  dx(0) * val; // -theta_z * dx
                    }
                        break;

                    case 2: {
                        // index of the rotational dof of the master node
                        c_row[master_dofs[3]]  = -dx(1) * val; //  theta_x * dy
                        c_row[master_dofs[4]]  =  dx(0) * val; // -theta_y * dx
                    }
                        break;
                }
            }
        }
    }
}

