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

#ifndef __mast__kinematic_coupling_h__
#define __mast__kinematic_coupling_h__


// C++ includes
#include <set>
#include <map>

// MAST includes
#include "base/dof_coupling_base.h"

// libMesh includes
#include "libmesh/node.h"

namespace MAST {
    
    // Forward declerations
    class KinematicCouplingConstraint;
    class MeshCouplingBase;
    
    /*!
     *  This constrains the slave nodes to be kinematically constrained
     *  to the master node. 
     */
     
    class KinematicCoupling:
    public MAST::DoFCouplingBase {

    public:
        
        KinematicCoupling(MAST::SystemInitialization& sys_init);
        
        virtual ~KinematicCoupling();

        /*!
         *   Constrain the slave node to the specified set of master nodes.
         *   A geometric relation is used to compute the weights for each
         *   master dof. If \p constrain_rotations is \p true then the
         *   rotations of the slave nodes are constrained to those of the master
         *   nodes. Otherwise, the constraint acts as a pinned connection.
         */
        void
        add_master_and_slave(const std::set<const libMesh::Node*>& master,
                             const libMesh::Node&            slave,
                             bool                            constrain_rotations);
        

        /*!
         *   provides a vector of the pair of slave and set of master nodes.
         */
        void
        add_master_and_slave
        (const std::vector<std::pair<const libMesh::Node*, std::set<const libMesh::Node*>>>& couplings,
         bool constrain_rotations);
        
        
        /*!
         *   initializes the dof constraint rows in the DofMap object associated
         *   with the system.
         */
        void add_dof_constraints_to_dof_map();

        /*!
         * Provides the vector of constraints and right-hand-side value pairs
         * to be added to the system.
         */
        virtual void
        get_constraint_rows(std::vector<std::tuple<libMesh::dof_id_type, libMesh::DofConstraintRow, Real>>& constrs);
        
        
    protected:
        
        bool _initialized;
        
        std::map<const libMesh::Node*, const MAST::KinematicCouplingConstraint*> _slave_node_constraints;
    };
}

#endif // __mast__kinematic_coupling_h__
