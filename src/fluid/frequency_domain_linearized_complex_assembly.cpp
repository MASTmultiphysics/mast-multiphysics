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
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/frequency_domain_linearized_conservative_fluid_elem.h"
#include "base/assembly_base.h"



MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
FrequencyDomainLinearizedComplexAssemblyElemOperations():
MAST::ComplexAssemblyElemOperations(),
_frequency(nullptr) {
    
}




MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
~FrequencyDomainLinearizedComplexAssemblyElemOperations() {
    
}




void
MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
set_frequency_function(MAST::FrequencyFunction& f) {
    
    // make sure that is hasn't been set
    libmesh_assert(!_frequency);
    
    _frequency = &f;
}




void
MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
clear_frequency_function() {

    _frequency = nullptr;
}



void
MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
elem_calculations(bool if_jac,
                  ComplexVectorX& vec,
                  ComplexMatrixX& mat) {
    
    libmesh_assert(_physics_elem);

    MAST::FrequencyDomainLinearizedConservativeFluidElem& e =
    dynamic_cast<MAST::FrequencyDomainLinearizedConservativeFluidElem&>(*_physics_elem);
    
    vec.setZero();
    mat.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac, vec, mat, _discipline->side_loads());
}



void
MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
elem_sensitivity_calculations(ComplexVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::FrequencyDomainLinearizedConservativeFluidElem& e =
    dynamic_cast<MAST::FrequencyDomainLinearizedConservativeFluidElem&>(*_physics_elem);
    
    vec.setZero();
    ComplexMatrixX
    dummy = ComplexMatrixX::Zero(vec.size(), vec.size());
    
    // assembly of the flux terms
    e.internal_residual_sensitivity(false, vec, dummy);
    e.side_external_residual_sensitivity(false, vec, dummy, _discipline->side_loads());
}



void
MAST::FrequencyDomainLinearizedComplexAssemblyElemOperations::
init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline&>
    (_assembly->discipline()).flight_condition();
    
    FrequencyDomainLinearizedConservativeFluidElem
    *freq_elem =
    new MAST::FrequencyDomainLinearizedConservativeFluidElem(*_system, *_assembly, elem, p);
    freq_elem->freq   = _frequency;
    
    _physics_elem = freq_elem;
}

