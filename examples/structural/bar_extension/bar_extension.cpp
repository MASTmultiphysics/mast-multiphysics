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
#include "examples/structural/bar_extension/bar_extension.h"
#include "base/physics_discipline_base.h"



MAST::Examples::BarExtension::BarExtension(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExample1D(comm_in) {
    
}


void
MAST::Examples::BarExtension::_init_dirichlet_conditions() {
    
    // constrain only the left
    this->_init_boundary_dirichlet_constraint(0, "left_constraint");
    
    _discipline->init_system_dirichlet_bc(*_sys);
}


void
MAST::Examples::BarExtension::_init_loads() {
    
    _init_pressure_load(true, 1);
}


