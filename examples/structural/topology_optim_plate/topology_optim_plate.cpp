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
#include "examples/structural/topology_optim_plate/topology_optim_plate.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/physics_discipline_base.h"
#include "heat_conduction/heat_conduction_system_initialization.h"


MAST::Examples::TopologyOptimizationLevelSetPlate::
TopologyOptimizationLevelSetPlate(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::TopologyOptimizationLevelSet2D(comm_in) {
    
}



MAST::Examples::TopologyOptimizationLevelSetPlate::~TopologyOptimizationLevelSetPlate() { }



void
MAST::Examples::TopologyOptimizationLevelSetPlate::_init_dirichlet_conditions() {
    
    MAST::Examples::StructuralExample2D::_init_dirichlet_conditions();
    
    ///////////////////////////////////////////////////////////////////////
    // initialize Dirichlet conditions for structural system
    ///////////////////////////////////////////////////////////////////////
    _discipline->init_system_dirichlet_bc(*_sys);
    
    ///////////////////////////////////////////////////////////////////////
    // initialize Dirichlet conditions for indicator system
    ///////////////////////////////////////////////////////////////////////
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
    dirichlet->init(1, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(1,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    dirichlet   = new MAST::DirichletBoundaryCondition;   // left boundary
    dirichlet->init(3, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(3,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
}



void
MAST::Examples::TopologyOptimizationLevelSetPlate::_init_loads() {
    
    
}
