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

#ifndef __mast__dof_coupling_base_h__
#define __mast__dof_coupling_base_h__

// MAST includes
#include "base/mast_data_types.h"

// libMesh includes
#include "libmesh/dof_map.h"


namespace MAST {
    
    // Forward declerations
    class SystemInitialization;
    
    /*!
     *  This provides a base class to couple degrees-of-freedom within a single
     *  system. This makes use of libMesh::DofConstraintRow as provided in
     *  libMesh::DofMap. For the \f$ k^{th}\f$ dof , \f$ u_k \f$, the
     *  libMesh::DofConstraintRow is defined as a linear equation
     *  \f[ u_k + \sum_{i=0, i\neq k}^{n-1} a_{u_i} u_i = b_k, \f]
     *  where, \f$ n \f$ is the number of degrees-of-freedom in the system.
     *  Note that a DofConstraintRow only requires non-zero \f$ a_{u_i} \f$ to
     *  be specified in the constraint.
     */
    
    class DoFCouplingBase {

    public:
        
        DoFCouplingBase(MAST::SystemInitialization& sys_init);
        
        virtual ~DoFCouplingBase();
        
        /*!.
         * Provides the vector of constraints and right-hand-side value pairs
         * to be added to the system.
         */
        virtual void
        get_constraint_rows(std::vector<std::tuple<libMesh::dof_id_type, libMesh::DofConstraintRow, Real>>& constrs) = 0;
        
    protected:
        
        MAST::SystemInitialization& _system_init;
    };
    
    
}

#endif // __mast__dof_coupling_base_h__
