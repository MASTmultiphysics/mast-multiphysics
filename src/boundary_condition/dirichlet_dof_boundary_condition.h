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

#ifndef __mast__dirichlet_dof_boundary_condition_h__
#define __mast__dirichlet_dof_boundary_condition_h__

// C++ includes
#include <memory>

// MAST includes
#include "base/boundary_condition_base.h"


// libmesh includes
#include "libmesh/system.h"
#include "libmesh/node.h"

/*!
 * CONCISE USAGE EXAMPLE
 * ---------------------
 * Define a std::vector of the nodes the boundary conditoin applies to
 * std::vector<const libMesh::Node*> bc_nodes{mesh.point_locator().locate_node(libMesh::Point(L0, W0+0.5*W, 0.)), mesh.point_locator().locate_node(libMesh::Point(L0+L, W0+0.5*W, 0.)) };
 * 
 * Define the local degrees of freedom the boundary condition applies to
 * std::vector<libMesh::dof_id_type> bc_dofs{0, 1, 2, 3};
 * 
 * Define the value of the boundary condition at each degree of freedom
 * std::vector<Real> bc_vals{0., 0., 0.};
 * 
 * Create the boundary condition object for dirichlet dof BCs
 * MAST::DOFDirichletBoundaryCondition  dof_bc;
 * 
 * Initialize the boundary condition with the vectors of nodes, dofs, and values
 * dof_bc.init(bc_nodes, bc_dofs, bc_vals);
 * 
 * Add the boundary condition to the discipline
 * discipline.add_dirichlet_dof_bc(dof_bc);
 * 
 * Repeat for another boundary condition
 * std::vector<const libMesh::Node*> bc_nodes3{mesh.point_locator().locate_node(libMesh::Point(L0, W0, 0.)) };
 * std::vector<libMesh::dof_id_type> bc_dofs3{0, 1, 2, 3, 4, 5};
 * std::vector<Real> bc_vals3{0., 0., 0., 0., 0., 0.};
 * MAST::DOFDirichletBoundaryCondition dof_bc3;
 * dof_bc3.init(bc_nodes3, bc_dofs3, bc_vals3);
 * discipline.add_dirichlet_dof_bc(dof_bc3);
 * 
 * Initialize the boundary conditions within the system
 * discipline.init_system_dirichlet_dof_bc(system);
 */


namespace MAST {
    
    // Forward declerations
    class NonlinearSystem;
    
    
    /*!
     * This is a boundary condition class which stores the a vector of the
     * pointers to libMesh::Node's that the BC is acting on, a vector of the
     * local degrees of freedom (0-5) the BC acts on, and a vector of the 
     * constrain values (which can be nonzero).
     */
    class DOFDirichletBoundaryCondition:  public MAST::BoundaryConditionBase 
    {
    public:
        DOFDirichletBoundaryCondition():
        MAST::BoundaryConditionBase(MAST::DOF_DIRICHLET)
        { }
        
        virtual ~DOFDirichletBoundaryCondition() { }


        void init(std::vector<const libMesh::Node*> bc_nodes_in,
                  std::vector<libMesh::dof_id_type> bc_dofs_in,
                  std::vector<Real> bc_vals_in);
        
        
        std::vector<const libMesh::Node*>& bc_nodes(){
            return _bc_nodes;
        }
        
        std::vector<libMesh::dof_id_type>& bc_dofs(){
            return _bc_dofs;
        }
        
        std::vector<Real>& bc_vals(){
            return _bc_vals;
        }
        
    protected:
        std::vector<const libMesh::Node*> _bc_nodes;
        std::vector<libMesh::dof_id_type> _bc_dofs;
        std::vector<Real> _bc_vals;
    };
    
    
    
    /*! 
     * This class is used by libMesh to impose the constraints by calling
     * the constrain() function before solving. The object created from this
     * class is added to the system through the attach_constraint_object 
     * method.
     */
    class DOFConstraint:  public libMesh::System::Constraint 
    {
    public:
        DOFConstraint(std::set<MAST::DOFDirichletBoundaryCondition*> dof_bcs);
        
        virtual ~DOFConstraint();
        
        virtual void constrain();
        
        virtual void setNonlinearSystem(MAST::NonlinearSystem& nonlinear_system_in);
        
    protected:
        MAST::NonlinearSystem* nonlinear_system = nullptr;
        const std::set<MAST::DOFDirichletBoundaryCondition*> _dof_bcs;
    };
    
    
}

#endif // __mast__dirichlet_dof_boundary_condition_h__
