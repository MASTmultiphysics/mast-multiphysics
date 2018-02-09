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
#include "elasticity/structural_buckling_eigenproblem_elem_operations.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/assembly_base.h"
#include "mesh/local_elem_fe.h"

MAST::StructuralBucklingEigenproblemElemOperations::
StructuralBucklingEigenproblemElemOperations():
MAST::EigenproblemAssemblyElemOperations() {
    
}



MAST::StructuralBucklingEigenproblemElemOperations::
~StructuralBucklingEigenproblemElemOperations() {
    
}


void
MAST::StructuralBucklingEigenproblemElemOperations::init(const libMesh::Elem& elem) {
    
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
MAST::StructuralBucklingEigenproblemElemOperations::
elem_calculations(RealMatrixX& mat_A,
                  RealMatrixX& mat_B) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealVectorX
    vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    mat_A.setZero();
    
    // calculate the Jacobian components
    e.internal_residual(true, vec, mat_A);
    e.prestress_residual(true, vec, mat_A);
    e.side_external_residual(true,
                             vec,
                             dummy,
                             mat_A,
                             _discipline->side_loads());
    e.volume_external_residual(true,
                               vec,
                               dummy,
                               mat_A,
                               _discipline->volume_loads());
}




void
MAST::StructuralBucklingEigenproblemElemOperations::
elem_sensitivity_calculations(bool base_sol,
                              RealMatrixX& mat_A,
                              RealMatrixX& mat_B) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(*_physics_elem);
    
    RealVectorX
    vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    mat_A.setZero();
    
    // calculate the Jacobian components
    e.internal_residual_sensitivity(true, vec, mat_A);
    e.prestress_residual_sensitivity(true, vec, mat_A);
    e.side_external_residual_sensitivity(true,
                                         vec,
                                         dummy,
                                         mat_A,
                                         _discipline->side_loads());
    e.volume_external_residual_sensitivity(true,
                                           vec,
                                           dummy,
                                           mat_A,
                                           _discipline->volume_loads());
}



void
MAST::StructuralBucklingEigenproblemElemOperations::
set_local_fe_data(MAST::LocalElemFE& fe,
                  const libMesh::Elem& e) const {
    
    if (e.dim() == 1) {
        
        const MAST::ElementPropertyCard1D&
        p_card = dynamic_cast<const MAST::ElementPropertyCard1D&>
        (_discipline->get_property_card(e));
        
        fe.set_1d_y_vector(p_card.y_vector());
    }
}


