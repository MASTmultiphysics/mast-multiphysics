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
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "numerics/utility.h"
#include "base/system_initialization.h"
#include "mesh/fe_base.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"


MAST::StructuralModalEigenproblemAssembly::
StructuralModalEigenproblemAssembly():
MAST::EigenproblemAssembly() { }




MAST::StructuralModalEigenproblemAssembly::
~StructuralModalEigenproblemAssembly()
{ }



std::unique_ptr<MAST::FEBase>
MAST::StructuralModalEigenproblemAssembly::build_fe(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    return std::unique_ptr<MAST::FEBase>(MAST::build_structural_fe(*_system, elem, p));
}


std::unique_ptr<MAST::ElementBase>
MAST::StructuralModalEigenproblemAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, *this, elem, p).release();
    
    return std::unique_ptr<MAST::ElementBase>(rval);
}



void
MAST::StructuralModalEigenproblemAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   RealMatrixX& mat_A,
                   RealMatrixX& mat_B) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
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
MAST::StructuralModalEigenproblemAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealMatrixX& mat_A,
                               RealMatrixX& mat_B) {
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealVectorX vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX mat = RealMatrixX::Zero(mat_A.rows(), mat_A.cols()); // dummy matrix
    mat_A.setZero();
    mat_B.setZero();
    
    // calculate the Jacobian components
    e.internal_residual_sensitivity(true, vec, mat_A);
    e.side_external_residual_sensitivity(true, vec, mat, mat_A, _discipline->side_loads());
    e.volume_external_residual_sensitivity(true, vec, mat, mat_A, _discipline->volume_loads());
    
    // calculate the mass matrix components
    e.inertial_residual_sensitivity(true, vec, mat_B, mat, mat_A);
    
    // if the linearization is about a base state, then the sensitivity of
    // the base state will influence the sensitivity of the Jacobian
    if (_base_sol)
        e.internal_residual_jac_dot_state_sensitivity(mat_A);
}


void
MAST::StructuralModalEigenproblemAssembly::_set_elem_sol(MAST::ElementBase& elem,
                                                         const RealVectorX& sol) {
    
    unsigned int
    n = (unsigned int)sol.size();
    
    RealVectorX
    zero = RealVectorX::Zero(n);
    
    elem.set_solution    (sol);
    elem.set_velocity    (zero); // set to zero vector for a quasi-steady analysis
    elem.set_acceleration(zero); // set to zero vector for a quasi-steady analysis
    
    
    // set the incompatible mode solution if required by the
    // element
    MAST::StructuralElementBase& s_elem =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    if (s_elem.if_incompatible_modes()) {
        
        // check if the vector exists in the map
        if (!_incompatible_sol.count(&elem.elem()))
            _incompatible_sol[&elem.elem()] = RealVectorX::Zero(s_elem.incompatible_mode_size());
        s_elem.set_incompatible_mode_solution(_incompatible_sol[&elem.elem()]);
    }
}


