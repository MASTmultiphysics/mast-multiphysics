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

// MAST includes
#include "boundary_condition/dirichlet_dof_boundary_condition.h"
#include "base/nonlinear_system.h"

// libMesh Includes
#include "libmesh/dof_map.h"


void MAST::DOFDirichletBoundaryCondition::init(std::vector<const libMesh::Node*> bc_nodes_in,
                  std::vector<libMesh::dof_id_type> bc_dofs_in,
                  std::vector<Real> bc_vals_in)
{
    _bc_nodes = bc_nodes_in;
    _bc_dofs =  bc_dofs_in;
    _bc_vals =  bc_vals_in;
}




MAST::DOFConstraint::DOFConstraint(std::set<MAST::DOFDirichletBoundaryCondition*> dof_bcs):
    libMesh::System::Constraint(), _dof_bcs(dof_bcs) {}
    
MAST::DOFConstraint::~DOFConstraint() {}

/*!
 * This method is called by libMesh prior to the solve command. It loops through
 * all the MAST::DOFDirichletBoundaryCondition objects, and applies the 
 * boundary condition (defined by nodes, DOFs, and values stores in the 
 * MAST::DOFDirichletBoundaryCondition object) to the system.
 */
void MAST::DOFConstraint::constrain()
{
    libMesh::DofConstraintRow c_row;
    libMesh::DofMap& dof_map = nonlinear_system->get_dof_map();
    std::vector<libMesh::dof_id_type> di;
    
    std::set<MAST::DOFDirichletBoundaryCondition*>::iterator it = _dof_bcs.begin();
    
    for ( ; it != _dof_bcs.end(); it++)
    {
        std::vector<const libMesh::Node*>& bc_nodes = (*it)->bc_nodes();
        std::vector<libMesh::dof_id_type>& bc_dofs = (*it)->bc_dofs();
        std::vector<Real>& bc_vals = (*it)->bc_vals();
        
        for (uint j=0; j<bc_nodes.size(); j++)
        {
            for (uint i=0; i<bc_dofs.size(); i++)
            {   
//                 const auto node = bc_nodes[j];
//                 int dof = bc_dofs[i];
//                 Real val = bc_vals[i];
//                 int global_dof = node->id()*6+dof;
//                 libMesh::out << "(" << node->id() << ",(" << dof << "," << global_dof << ")," << val << ")" << std::endl;
                di.clear();
                dof_map.dof_indices(bc_nodes[j], di, bc_dofs[i]);
                libmesh_assert_equal_to(di.size(), 1);
                dof_map.add_constraint_row(di[0], c_row, bc_vals[i], true);
            }
        }
    }
}

/*!
 * This method is used to set the MAST::NonlinearSystem object that the 
 * constraints are acting on.
 */
void MAST::DOFConstraint::setNonlinearSystem(MAST::NonlinearSystem& nonlinear_system_in)
{
    nonlinear_system = &nonlinear_system_in;
}
