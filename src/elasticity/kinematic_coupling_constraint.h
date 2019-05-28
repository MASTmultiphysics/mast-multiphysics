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

#ifndef __mast__kinematic_coupling_constraint_h__
#define __mast__kinematic_coupling_constraint_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/node.h"
#include "libmesh/dof_map.h"


namespace MAST {

    // Forward declerations
    class SystemInitialization;

    /*!
     *  This object stores the information about the coupling of nodes.
     *  A slave node can be constrained to one or more nodes and the rotations
     *  of the slave can be constrained (set equal to) the rotation of
     *  the master node. In case more than one master nodes are specified, then
     *  the coefficients in the constraint equation are weighted based on
     *  inverse relationship with the geometric distance and partition-of-unity.
     *
     *  The constraints are defined such that the displacement of the slave nodes
     *  is identified as
     *  \f[  u_s = u_m + \theta \times (X_s - X_m),   \f]
     *  where, \f$ \theta \f$ is the rotation vector
     *  \f$(\theta_x, \theta_y, \theta_z)\f$,
     *  \f$ u_s \f$ and \f$ X_s \f$ are the translation and location
     *  of the slave nodes, respectively, and likewise for the master node
     *  with subscript \f$ (\cdot)_m \f$. Note, that the rotations of the
     *  slave node are not constrained to those of the master node. The
     *  cross product is written as
     *  \f{eqnarray*}{ \theta \times dX
     *  & = & \left| \begin{array}{ccc}
     *           \hat{i} & \hat{j} & \hat{k} \\
     *           \theta_x & \theta_y & \theta_z \\
     *           dx & dy & dz
     *        \end{array} \right| \\
     *  & = &    \hat{i} (\theta_y dz - \theta_z dy)
     *         - \hat{j} (\theta_x dz - \theta_z dx)
     *         + \hat{k} (\theta_x dy - \theta_y dx)
     *  \f}
     *   Therefore, the constraint for displacement along the x-axis is
     *  \f[ u_{xs} = u_{xm} + (\theta_y dz - \theta_z dy),  \f]
     *   for displacement along the y-axis,
     *  \f[ u_{ys} = u_{ym} + (\theta_x dz - \theta_z dx),  \f]
     *   and for displacement along the z-axis,
     *  \f[ u_{zs} = u_{zm} + (\theta_x dy - \theta_y dx).  \f]
     */
    class KinematicCouplingConstraint {
      
    public:
        
        KinematicCouplingConstraint(MAST::SystemInitialization& sys_init,
                                    const libMesh::Node& slave_node,
                                    std::set<const libMesh::Node*>& master_nodes,
                                    bool constrain_rotations);
        
        virtual ~KinematicCouplingConstraint();

        /*!
         *  @return the number of constraints that will be added from this object
         */
        unsigned int
        n_constrain_nodes() const {
            
            return _constrain_rotations?6:3;
        }

        /*!
         *   initializes the vector of \p libMesh::DofConstraintRow objects
         *   and rhs values for this node
         */
        void
        get_dof_constraint_row(std::vector<std::tuple
                               <libMesh::dof_id_type,
                               libMesh::DofConstraintRow,
                               Real>>& constraints) const;
        
    protected:
        
        MAST::SystemInitialization&    _sys_init;
        const libMesh::Node*           _slave;
        std::set<const libMesh::Node*> _masters;
        bool                           _constrain_rotations;
    };
}

#endif // __mast__kinematic_coupling_constraint_h__
