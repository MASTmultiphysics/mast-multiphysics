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



void
MAST::StructuralModalEigenproblemAssembly::
eigenproblem_assemble(libMesh::SparseMatrix<Real> *A,
                      libMesh::SparseMatrix<Real> *B)  {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *A,
    &matrix_B = *B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // build localized solutions if needed
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution;
    
    if (_base_sol)
        localized_solution.reset(_build_localized_vector(eigen_sys,
                                                         *_base_sol).release());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dummy;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::unique_ptr<MAST::ElementBase> physics_elem;
    
    libMesh::MeshBase::const_element_iterator       el     =
    eigen_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    eigen_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;

        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dummy.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(dummy);
        physics_elem->set_acceleration(dummy);
        
        
        // set the incompatible mode solution if required by the
        // element
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }

        _elem_calculations(*physics_elem, mat_A, mat_B);
        
        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);
        
        // constrain the element matrices.
        dof_map.constrain_element_matrix(A, dof_indices);
        dof_map.constrain_element_matrix(B, dof_indices);
        
        // add to the global matrices
        matrix_A.add_matrix (A, dof_indices); // load independent
        matrix_B.add_matrix (B, dof_indices); // load dependent
    }
    
    // finalize the data structures
    A->close();
    B->close();
}





bool
MAST::StructuralModalEigenproblemAssembly::
eigenproblem_sensitivity_assemble (const libMesh::ParameterVector& parameters,
                                   const unsigned int i,
                                   libMesh::SparseMatrix<Real>* sensitivity_A,
                                   libMesh::SparseMatrix<Real>* sensitivity_B)  {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());
    
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *sensitivity_A,
    &matrix_B = *sensitivity_B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // build localized solutions if needed
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    
    if (_base_sol) {
        localized_solution.reset(_build_localized_vector(eigen_sys,
                                                         *_base_sol).release());
        
        // make sure that the sensitivity was also provided
        libmesh_assert(_base_sol_sensitivity);
        localized_solution_sens.reset(_build_localized_vector(eigen_sys,
                                                              *_base_sol_sensitivity).release());
    }

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dummy;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::unique_ptr<MAST::ElementBase> physics_elem;
    
    libMesh::MeshBase::const_element_iterator       el     =
    eigen_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    eigen_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dummy.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol) {
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(dummy);
        physics_elem->set_acceleration(dummy);
        
        // set the element's base solution sensitivity
        if (_base_sol) {
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution_sens)(dof_indices[i]);
        }
        physics_elem->set_solution(sol, true);
        
        // set the incompatible mode solution if required by the
        // element
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }
        
        _elem_sensitivity_calculations(*physics_elem, mat_A, mat_B);
        
        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);
        
        // constrain the element matrices.
        dof_map.constrain_element_matrix(A, dof_indices);
        dof_map.constrain_element_matrix(B, dof_indices);
        
        // add to the global matrices
        matrix_A.add_matrix (A, dof_indices); // load independent
        matrix_B.add_matrix (B, dof_indices); // load dependent
    }
    
    // finalize the data structures
    sensitivity_A->close();
    sensitivity_B->close();

    return true;
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


