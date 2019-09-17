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
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/assembly_base.h"
#include "base/physics_discipline_base.h"
#include "heat_conduction/heat_conduction_elem_base.h"
#include "property_cards/element_property_card_1D.h"
#include "level_set/level_set_intersected_elem.h"


MAST::HeatConductionNonlinearAssemblyElemOperations::
HeatConductionNonlinearAssemblyElemOperations():
MAST::NonlinearImplicitAssemblyElemOperations() {
    
}



MAST::HeatConductionNonlinearAssemblyElemOperations::
~HeatConductionNonlinearAssemblyElemOperations() {
    
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
set_elem_data(unsigned int dim,
              const libMesh::Elem& ref_elem,
              MAST::GeomElem& elem) const {
    
    libmesh_assert(!_physics_elem);
    
    if (dim == 1) {
        
        const MAST::ElementPropertyCard1D& p =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_discipline->get_property_card(ref_elem));
        
        elem.set_local_y_vector(p.y_vector());
    }
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
init(const MAST::GeomElem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>
    (_discipline->get_property_card(elem));
    
    _physics_elem =
    new MAST::HeatConductionElementBase(*_system, *_assembly, elem, p);
}




void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_calculations(bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {

    libmesh_assert(_physics_elem);

    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(*_physics_elem);
    
    vec.setZero();
    mat.setZero();
    
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac, vec, mat, _discipline->side_loads());
    e.volume_external_residual(if_jac, vec, mat, _discipline->volume_loads());
}




void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::HeatConductionElementBase& e =
    dynamic_cast<MAST::HeatConductionElementBase&>(*_physics_elem);
    
    vec.setZero();
    
    e.internal_residual_sensitivity(f, vec);
    e.side_external_residual_sensitivity(f, vec, _discipline->side_loads());
    e.volume_external_residual_sensitivity(f, vec, _discipline->volume_loads());
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());
    
    std::pair<const MAST::FieldFunction<RealVectorX>*, unsigned int>
    val = this->get_elem_boundary_velocity_data();
    
    if (val.first) {
        
        MAST::HeatConductionElementBase& e =
        dynamic_cast<MAST::HeatConductionElementBase&>(*_physics_elem);
        
        vec.setZero();
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.internal_residual_boundary_velocity(f, vec,
                                              val.second,
                                              *val.first);
        e.volume_external_residual_boundary_velocity(f, vec,
                                                     val.second,
                                                     *val.first,
                                                     _discipline->volume_loads());
        /*e.side_external_residual_sensitivity(f, false,
         vec,
         dummy,
         dummy,
         _discipline->side_loads());*/
    }
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       const MAST::FieldFunction<RealVectorX>& vel,
                                       RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    const MAST::LevelSetIntersectedElem
    &elem = dynamic_cast<const MAST::LevelSetIntersectedElem&>(_physics_elem->elem());

    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with a side.
    if (elem.if_elem_has_level_set_boundary() &&
        elem.if_subelem_has_side_on_level_set_boundary()) {
        
        MAST::HeatConductionElementBase& e =
        dynamic_cast<MAST::HeatConductionElementBase&>(*_physics_elem);
        
        vec.setZero();
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.internal_residual_boundary_velocity(f, vec,
                                              elem.get_subelem_side_on_level_set_boundary(),
                                              vel);
        e.volume_external_residual_boundary_velocity(f, vec,
                                                     elem.get_subelem_side_on_level_set_boundary(),
                                                     vel,
                                                     _discipline->volume_loads());
        /*e.side_external_residual_sensitivity(f, false,
         vec,
         dummy,
         dummy,
         _discipline->side_loads());*/
    }
}


void
MAST::HeatConductionNonlinearAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_error(); // to be implemented
}

