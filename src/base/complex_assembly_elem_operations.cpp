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
#include "base/complex_assembly_elem_operations.h"
#include "base/elem_base.h"


MAST::ComplexAssemblyElemOperations::ComplexAssemblyElemOperations():
MAST::AssemblyElemOperations() {
    
}



MAST::ComplexAssemblyElemOperations::~ComplexAssemblyElemOperations() {
    
}


void
MAST::ComplexAssemblyElemOperations::set_elem_complex_solution(const ComplexVectorX &sol) {
    
    _physics_elem->set_complex_solution(sol);
}


void
MAST::ComplexAssemblyElemOperations::set_elem_complex_solution_sensitivity(const ComplexVectorX &sol) {
    
    _physics_elem->set_complex_solution(sol, true);
}

