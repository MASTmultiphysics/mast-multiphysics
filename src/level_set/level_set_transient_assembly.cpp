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
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_elem_base.h"
#include "property_cards/element_property_card_1D.h"
#include "level_set/level_set_discipline.h"
#include "mesh/local_elem_fe.h"
#include "base/assembly_base.h"


MAST::LevelSetTransientAssemblyElemOperations::
LevelSetTransientAssemblyElemOperations():
MAST::TransientAssemblyElemOperations() {
    
}




MAST::LevelSetTransientAssemblyElemOperations::
~LevelSetTransientAssemblyElemOperations() {
    
}


void
MAST::LevelSetTransientAssemblyElemOperations::
set_elem_reference_solution(const RealVectorX& sol) {
    
    libmesh_assert(_physics_elem);
    
    dynamic_cast<MAST::LevelSetElementBase*>(_physics_elem)->
    set_reference_solution_for_initialization(sol);
}



void
MAST::LevelSetTransientAssemblyElemOperations::
elem_calculations(bool if_jac,
                  RealVectorX& f_m,
                  RealVectorX& f_x,
                  RealMatrixX& f_m_jac_xdot,
                  RealMatrixX& f_m_jac,
                  RealMatrixX& f_x_jac) {
    
    libmesh_assert(_physics_elem);

    MAST::LevelSetElementBase& e =
    dynamic_cast<MAST::LevelSetElementBase&>(*_physics_elem);
    
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
MAST::LevelSetTransientAssemblyElemOperations::
linearized_jacobian_solution_product(RealVectorX& f) {
    
    libmesh_error(); // should not get called
}



void
MAST::LevelSetTransientAssemblyElemOperations::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& vec) {
    
    libmesh_error(); // should not get called
    
}


void
MAST::LevelSetTransientAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_error(); // should not get called
}


void
MAST::LevelSetTransientAssemblyElemOperations::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& f_m,
                              RealVectorX& f_x) {
    libmesh_error(); // should not get called
}



void
MAST::LevelSetTransientAssemblyElemOperations::
init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::LevelSetDiscipline&
    discipline = dynamic_cast<const MAST::LevelSetDiscipline&>(*_discipline);
    
    MAST::LevelSetElementBase
    *e = new MAST::LevelSetElementBase(*_system, *_assembly, elem);
    
    if (discipline.has_velocity_function())
        e->set_velocity_function(discipline.get_velocity_function());
    e->set_propagation_mode(discipline.if_level_set_propagation());
    
    _physics_elem = e;
}


void
MAST::LevelSetTransientAssemblyElemOperations::
set_local_fe_data(MAST::LocalElemFE& fe,
                  const libMesh::Elem& e) const {
        
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}


