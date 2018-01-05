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
#include "base/assembly_elem_operation.h"
#include "base/assembly_base.h"
#include "mesh/fe_base.h"






MAST::AssemblyElemOperations::AssemblyElemOperations():
_assembly (nullptr) {
    
}


MAST::AssemblyElemOperations::~AssemblyElemOperations() {
    
}



void
MAST::AssemblyElemOperations::set_assembly(MAST::AssemblyBase &assembly) {

    libmesh_assert(!_assembly);
    _assembly = &assembly;
}


void
MAST::AssemblyElemOperations::clear_assembly() {
    
    _assembly = nullptr;
}


MAST::AssemblyBase&
MAST::AssemblyElemOperations::get_assembly() {
    
    libmesh_assert(_assembly);
    return *_assembly;
}

