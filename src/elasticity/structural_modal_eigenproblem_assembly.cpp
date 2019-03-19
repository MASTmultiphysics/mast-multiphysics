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
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/function_base.h"
#include "base/assembly_base.h"
#include "level_set/level_set_intersected_elem.h"



MAST::StructuralModalEigenproblemAssemblyElemOperations::
StructuralModalEigenproblemAssemblyElemOperations():
MAST::EigenproblemAssemblyElemOperations() { }




MAST::StructuralModalEigenproblemAssemblyElemOperations::
~StructuralModalEigenproblemAssemblyElemOperations()
{ }




void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
set_elem_solution(const RealVectorX& sol) {
    
    unsigned int
    n = (unsigned int)sol.size();
    
    RealVectorX
    zero = RealVectorX::Zero(n);
    
    _physics_elem->set_solution    (sol);
    _physics_elem->set_velocity    (zero); // set to zero vector for a quasi-steady analysis
    _physics_elem->set_acceleration(zero); // set to zero vector for a quasi-steady analysis
    
    
    if (_incompatible_sol_assembly) {
        
        // set the incompatible mode solution if required by the
        // element
        
        MAST::StructuralElementBase& s_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        if (s_elem.if_incompatible_modes())
            _incompatible_sol_assembly->set_elem_incompatible_sol(s_elem);
    }
}




void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
set_elem_solution_sensitivity(const RealVectorX& sol) {
    
    unsigned int
    n = (unsigned int)sol.size();
    
    RealVectorX
    zero = RealVectorX::Zero(n);
    
    _physics_elem->set_solution    (sol, true);
}



void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
elem_calculations(RealMatrixX& mat_A,
                  RealMatrixX& mat_B) {
    
    libmesh_assert(_physics_elem);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealVectorX vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX mat = RealMatrixX::Zero(mat_A.rows(), mat_A.cols()); // dummy matrix
    mat_A.setZero();
    mat_B.setZero();
    
    // calculate the Jacobian components
    e.internal_residual(true, vec, mat_A);
    e.side_external_residual(true, vec, mat, mat_A, _discipline->side_loads());
    e.volume_external_residual(true, vec, mat, mat_A, _discipline->volume_loads());
    
    // calculate the mass matrix components
    e.inertial_residual(true, vec, mat_B, mat, mat_A);
}




void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              bool base_sol,
                              RealMatrixX& mat_A,
                              RealMatrixX& mat_B) {
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealVectorX vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX mat = RealMatrixX::Zero(mat_A.rows(), mat_A.cols()); // dummy matrix
    mat_A.setZero();
    mat_B.setZero();
    
    // calculate the Jacobian components
    e.internal_residual_sensitivity(f, true, vec, mat_A);
    e.side_external_residual_sensitivity(f, true, vec, mat, mat_A, _discipline->side_loads());
    e.volume_external_residual_sensitivity(f, true, vec, mat, mat_A, _discipline->volume_loads());
    
    // calculate the mass matrix components
    e.inertial_residual_sensitivity(f, true, vec, mat_B, mat, mat_A);
    
    // if the linearization is about a base state, then the sensitivity of
    // the base state will influence the sensitivity of the Jacobian
    if (base_sol)
        e.internal_residual_jac_dot_state_sensitivity(mat_A);
}




void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       bool base_sol,
                                       const MAST::FieldFunction<RealVectorX>& vel,
                                       RealMatrixX& mat_A,
                                       RealMatrixX& mat_B) {
    
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());
    
    const MAST::LevelSetIntersectedElem
    &elem = dynamic_cast<const MAST::LevelSetIntersectedElem&>(_physics_elem->elem());

    RealVectorX vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX mat = RealMatrixX::Zero(mat_A.rows(), mat_A.cols()); // dummy matrix
    mat_A.setZero();
    mat_B.setZero();

    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with a side.
    if (elem.if_elem_has_level_set_boundary() &&
        elem.if_subelem_has_side_on_level_set_boundary()) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        e.internal_residual_boundary_velocity(f,
                                              elem.get_subelem_side_on_level_set_boundary(),
                                              vel,
                                              true,
                                              vec,
                                              mat_A);
        e.volume_external_residual_boundary_velocity(f,
                                                     elem.get_subelem_side_on_level_set_boundary(),
                                                     vel,
                                                     _discipline->volume_loads(),
                                                     true,
                                                     vec,
                                                     mat_A);
        
        e.inertial_residual_boundary_velocity(f,
                                              elem.get_subelem_side_on_level_set_boundary(),
                                              vel,
                                              true,
                                              vec,
                                              mat_B,
                                              mat,
                                              mat_A);
        
        // if the linearization is about a base state, then the sensitivity of
        // the base state will influence the sensitivity of the Jacobian
        if (base_sol)
            libmesh_assert(false); // to be implemented
        //e.internal_residual_jac_dot_state_sensitivity(mat_A);
    }
}



void
MAST::StructuralModalEigenproblemAssemblyElemOperations::set_elem_data(MAST::GeomElem& elem) const {
    
    libmesh_assert(!_physics_elem);
    
    if (elem.dim() == 1) {
        
        const MAST::ElementPropertyCard1D& p =
        dynamic_cast<const MAST::ElementPropertyCard1D&>(_discipline->get_property_card(elem));
        
        elem.set_local_y_vector(p.y_vector());
    }
}


void
MAST::StructuralModalEigenproblemAssemblyElemOperations::
init(const MAST::GeomElem& elem) {

    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>
    (_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}

