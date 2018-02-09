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
#include "heat_conduction/heat_conduction_transient_assembly.h"
#include "heat_conduction/heat_conduction_elem_base.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "mesh/local_elem_fe.h"
#include "base/assembly_base.h"


MAST::HeatConductionTransientAssemblyElemOperations::
HeatConductionTransientAssemblyElemOperations():
MAST::TransientAssemblyElemOperations() {
    
}




MAST::HeatConductionTransientAssemblyElemOperations::
~HeatConductionTransientAssemblyElemOperations() {
    
}



void
MAST::HeatConductionTransientAssemblyElemOperations::
elem_calculations(bool if_jac,
                  RealVectorX& f_m,
                  RealVectorX& f_x,
                  RealMatrixX& f_m_jac_xdot,
                  RealMatrixX& f_m_jac,
                  RealMatrixX& f_x_jac) {
    
    libmesh_assert(_physics_elem);

    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(*_physics_elem);
    
    f_m.setZero();
    f_x.setZero();
    f_m_jac_xdot.setZero();
    f_m_jac.setZero();
    f_x_jac.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, f_x, f_x_jac);
    e.side_external_residual(if_jac, f_x, f_x_jac, _discipline->side_loads());
    e.volume_external_residual(if_jac, f_x, f_x_jac, _discipline->volume_loads());
    
    //assembly of the capacitance term
    e.velocity_residual(if_jac, f_m, f_m_jac_xdot, f_m_jac);
}



void
MAST::HeatConductionTransientAssemblyElemOperations::
linearized_jacobian_solution_product(RealVectorX& f) {
    
    libmesh_error(); // to be implemented
}



void
MAST::HeatConductionTransientAssemblyElemOperations::
elem_sensitivity_calculations(RealVectorX& f_m,
                              RealVectorX& f_x) {
    
    libmesh_error(); // to be implemented
    
}


void
MAST::HeatConductionTransientAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_error(); // to be implemented
}



void
MAST::HeatConductionTransientAssemblyElemOperations::
init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    _physics_elem =
    new MAST::HeatConductionElementBase(*_system, *_assembly, elem, p);
}


void
MAST::HeatConductionTransientAssemblyElemOperations::
set_local_fe_data(MAST::LocalElemFE& fe,
                  const libMesh::Elem& e) const {
    
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}

