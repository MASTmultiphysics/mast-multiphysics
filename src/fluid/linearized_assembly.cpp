/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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
#include "fluid/linearized_assembly.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/linearized_conservative_fluid_elem.h"
#include "base/physics_discipline_base.h"



MAST::LinearizedAssembly::
LinearizedAssembly():
MAST::LinearizedAssemblyBase() {
    
}




MAST::LinearizedAssembly::
~LinearizedAssembly() {
    
}



void
MAST::LinearizedAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {
    
    MAST::LinearizedConservativeFluidElem& e =
    dynamic_cast<MAST::LinearizedConservativeFluidElem&>(elem);
    
    vec.setZero();
    mat.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac, vec, mat, _discipline->side_loads());
}



void
MAST::LinearizedAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               bool if_jac,
                               RealVectorX& vec,
                               RealMatrixX& mat) {
    
    
    MAST::LinearizedConservativeFluidElem& e =
    dynamic_cast<MAST::LinearizedConservativeFluidElem&>(elem);
    
    vec.setZero();
    mat.setZero();
    
    // assembly of the flux terms
    e.internal_residual_sensitivity(if_jac, vec, mat);
    e.side_external_residual_sensitivity(if_jac, vec, mat, _discipline->side_loads());
}



std::auto_ptr<MAST::ElementBase>
MAST::LinearizedAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::FlightCondition& p =
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(_discipline)->flight_condition();
    
    MAST::LinearizedConservativeFluidElem* rval =
    new MAST::LinearizedConservativeFluidElem(*_system, elem, p);
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}

