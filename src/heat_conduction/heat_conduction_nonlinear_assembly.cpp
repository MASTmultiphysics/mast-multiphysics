/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "heat_conduction/heat_conduction_elem_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"



MAST::HeatConductionNonlinearAssembly::
HeatConductionNonlinearAssembly():
MAST::NonlinearImplicitAssembly() {
    
}



MAST::HeatConductionNonlinearAssembly::
~HeatConductionNonlinearAssembly() {
    
}



std::auto_ptr<MAST::ElementBase>
MAST::HeatConductionNonlinearAssembly::_build_elem(const libMesh::Elem& elem) {

    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    new MAST::HeatConductionElementBase(*_system, elem, p);
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}




void
MAST::HeatConductionNonlinearAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {
    
    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(elem);
    
    vec.setZero();
    mat.setZero();
    
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac, vec, mat, _discipline->side_loads());
    e.volume_external_residual(if_jac, vec, mat, _discipline->volume_loads());
}




void
MAST::HeatConductionNonlinearAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealVectorX& vec) {
    
    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(elem);
    
    vec.setZero();
    RealMatrixX mat; // dummy matrix
    
    e.internal_residual_sensitivity(false, vec, mat);
    e.side_external_residual_sensitivity(false, vec, mat, _discipline->side_loads());
    e.volume_external_residual_sensitivity(false, vec, mat, _discipline->volume_loads());
}


