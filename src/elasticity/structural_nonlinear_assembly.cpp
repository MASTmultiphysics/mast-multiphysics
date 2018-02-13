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
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "mesh/local_elem_fe.h"
#include "level_set/level_set_intersection.h"



MAST::StructuralNonlinearAssemblyElemOperations::StructuralNonlinearAssemblyElemOperations():
MAST::NonlinearImplicitAssemblyElemOperations(),
_incompatible_sol_assembly(nullptr) {
    
}



MAST::StructuralNonlinearAssemblyElemOperations::
~StructuralNonlinearAssemblyElemOperations() {
    
}



void
MAST::StructuralNonlinearAssemblyElemOperations::set_elem_solution(const RealVectorX& sol) {
    
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
MAST::StructuralNonlinearAssemblyElemOperations::
init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);

    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>
    (_discipline->get_property_card(elem));
    
    _physics_elem =
    MAST::build_structural_element(*_system, *_assembly, elem, p).release();
}




void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_calculations(bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    
    libmesh_assert(_physics_elem);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    vec.setZero();
    mat.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(mat.rows(), mat.cols());
    
    e.internal_residual(if_jac, vec, mat);
    e.side_external_residual(if_jac,
                             vec,
                             dummy,
                             mat,
                             _discipline->side_loads());
    e.volume_external_residual(if_jac,
                               vec,
                               dummy,
                               mat,
                               _discipline->volume_loads());
}



void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_linearized_jacobian_solution_product(RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    vec.setZero();
    RealMatrixX
    mat = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.linearized_internal_residual(false, vec, mat);
    e.linearized_side_external_residual(false,
                                        vec,
                                        mat,
                                        mat,
                                        _discipline->side_loads());
    e.linearized_volume_external_residual(false,
                                          vec,
                                          mat,
                                          mat,
                                          _discipline->volume_loads());
}




void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    vec.setZero();
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    e.internal_residual_sensitivity(f, false, vec, dummy);
    e.side_external_residual_sensitivity(f, false,
                                         vec,
                                         dummy,
                                         dummy,
                                         _discipline->side_loads());
    e.volume_external_residual_sensitivity(f, false,
                                           vec,
                                           dummy,
                                           dummy,
                                           _discipline->volume_loads());
}



void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       const MAST::LevelSetIntersection& intersect,
                                       const MAST::FieldFunction<RealVectorX>& vel,
                                       RealVectorX& vec) {
    
    libmesh_assert(_physics_elem);
    libmesh_assert(f.is_topology_parameter());

    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with a side.
    if (intersect.if_elem_has_boundary()) {
        
        MAST::StructuralElementBase& e =
        dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
        
        vec.setZero();
        RealMatrixX
        dummy = RealMatrixX::Zero(vec.size(), vec.size());
        
        e.internal_residual_boundary_velocity(f, vec,
                                              intersect.get_side_on_interface(_physics_elem->elem()),
                                              vel);
        /*e.side_external_residual_sensitivity(f, false,
         vec,
         dummy,
         dummy,
         _discipline->side_loads());
         e.volume_external_residual_sensitivity(f, false,
         vec,
         dummy,
         dummy,
         _discipline->volume_loads());*/
    }
}



void
MAST::StructuralNonlinearAssemblyElemOperations::
elem_second_derivative_dot_solution_assembly(RealMatrixX& m) {
    
    libmesh_assert(_physics_elem);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    m.setZero();
    
    e.internal_residual_jac_dot_state_sensitivity(m);
}




void
MAST::StructuralNonlinearAssemblyElemOperations::
set_local_fe_data(MAST::LocalElemFE& fe,
                  const libMesh::Elem& e) const {
     
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}







