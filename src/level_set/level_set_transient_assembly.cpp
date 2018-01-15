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
elem_calculations(MAST::ElementBase& elem,
                  bool if_jac,
                  RealVectorX& f_m,
                  RealVectorX& f_x,
                  RealMatrixX& f_m_jac_xdot,
                  RealMatrixX& f_m_jac,
                  RealMatrixX& f_x_jac) {
    
    MAST::LevelSetElementBase& e =
    dynamic_cast<MAST::LevelSetElementBase&>(elem);
    
    f_m.setZero();
    f_x.setZero();
    f_m_jac_xdot.setZero();
    f_m_jac.setZero();
    f_x_jac.setZero();
    
    // assembly of the flux terms
    e.internal_residual(if_jac, f_x, f_x_jac);
    e.side_external_residual(if_jac, f_x, f_x_jac, _assembly->discipline().side_loads());
    e.volume_external_residual(if_jac, f_x, f_x_jac, _assembly->discipline().volume_loads());
    
    //assembly of the capacitance term
    e.velocity_residual(if_jac, f_m, f_m_jac_xdot, f_m_jac);
}



void
MAST::LevelSetTransientAssemblyElemOperations::
linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                     RealVectorX& f) {
    
    libmesh_error(); // to be implemented
}



void
MAST::LevelSetTransientAssemblyElemOperations::
elem_sensitivity_calculations(MAST::ElementBase& elem,
                              RealVectorX& vec) {
    
    libmesh_error(); // to be implemented
    
}


void
MAST::LevelSetTransientAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(MAST::ElementBase& elem,
                                             RealMatrixX& m) {
    
    libmesh_error(); // to be implemented
}



std::unique_ptr<MAST::ElementBase>
MAST::LevelSetTransientAssemblyElemOperations::
build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::FieldFunction<RealVectorX>& vel =
    (dynamic_cast<const MAST::LevelSetDiscipline&>
    (_assembly->discipline())).get_velocity_function();
    
    MAST::ElementBase* rval =
    new MAST::LevelSetElementBase
    (_assembly->system_init(), *_assembly, elem, vel);
    
    return std::unique_ptr<MAST::ElementBase>(rval);
}


void
MAST::LevelSetTransientAssemblyElemOperations::
set_local_fe_data(const libMesh::Elem& e,
                  MAST::LocalElemFE& fe) const {
    
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_assembly->discipline().get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}


