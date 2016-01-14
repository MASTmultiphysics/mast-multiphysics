/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "elasticity/structural_buckling_eigenproblem_assembly.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "numerics/utility.h"
#include "base/system_initialization.h"
#include "base/parameter.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"


MAST::StructuralBucklingEigenproblemAssembly::
StructuralBucklingEigenproblemAssembly():
MAST::EigenproblemAssembly(),
_use_linearized_formulation(true),
_load_param(NULL),
_lambda1(0.),
_lambda2(0.),
_sol1(NULL),
_sol2(NULL) {
    
}




MAST::StructuralBucklingEigenproblemAssembly::
~StructuralBucklingEigenproblemAssembly()
{ }




void
MAST::StructuralBucklingEigenproblemAssembly::
set_buckling_data (bool use_linearized_form,
                   MAST::Parameter& p,
                   const Real lambda1,
                   const Real lambda2,
                   libMesh::NumericVector<Real>& x1,
                   libMesh::NumericVector<Real>& x2) {
    
    _use_linearized_formulation = use_linearized_form;
    _load_param                 = &p;
    _lambda1                    = lambda1;
    _lambda2                    = lambda2;
    _sol1                       = &x1;
    _sol2                       = &x2;
}


void
MAST::StructuralBucklingEigenproblemAssembly::clear_discipline_and_system() {

    _use_linearized_formulation = true;
    _load_param                 = NULL;
    _lambda1                    = 0.;
    _lambda2                    = 0.;
    _sol1                       = NULL;
    _sol2                       = NULL;
    
    MAST::EigenproblemAssembly::clear_discipline_and_system();
}



void
MAST::StructuralBucklingEigenproblemAssembly::
eigenproblem_assemble(libMesh::SparseMatrix<Real> *A,
                      libMesh::SparseMatrix<Real> *B)  {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());

    // create the localized solution vectors
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution1,
    localized_solution2;
    
    localized_solution1.reset(_build_localized_vector(eigen_sys,
                                                      *_sol1).release());
    localized_solution2.reset(_build_localized_vector(eigen_sys,
                                                      *_sol2).release());
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *A,
    &matrix_B = *B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dummy;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
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
        
        ////////////////////////////////////////////////////////////////
        // first calculate the tangent stiffness matrix about the
        // first load factor and solution
        ////////////////////////////////////////////////////////////////
        (*_load_param) = _lambda1;
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution1)(dof_indices[i]);
    
        physics_elem->set_solution(sol);
        
        
        // set the incompatible mode solution if required by the
        // element
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
        if (p_elem.if_incompatible_modes()) {
            libmesh_error(); // this needs to be carefully implemented
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }

        DenseRealMatrix AA, BB;
        _elem_calculations(*physics_elem, mat_A);
        
        MAST::copy(AA, mat_A); // copy to the libMesh matrix for further processing
        eigen_sys.get_dof_map().constrain_element_matrix(AA, dof_indices); // constrain the element matrices.
        matrix_A.add_matrix (AA, dof_indices); // add to the global matrices
        
        
        
        
        ////////////////////////////////////////////////////////////////
        // next, calculate the tangent stiffness matrix about the
        // second load factor and solution
        ////////////////////////////////////////////////////////////////
        (*_load_param) = _lambda2;
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution2)(dof_indices[i]);
        
        physics_elem->set_solution(sol);
        
        
        // set the incompatible mode solution if required by the
        // element
        if (p_elem.if_incompatible_modes()) {
            
            libmesh_error(); // this needs to be carefully implemented
                             // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }

        _elem_calculations(*physics_elem, mat_B);
        mat_B  *= -1.;
        if (_use_linearized_formulation)
            mat_B  += mat_A;
        
        MAST::copy(BB, mat_B); // copy to the libMesh matrix for further processing
        eigen_sys.get_dof_map().constrain_element_matrix(BB, dof_indices); // constrain the element matrices.
        matrix_B.add_matrix (BB, dof_indices); // add to the global matrices
    }
    
    // finalize the matrices for futher use.
    A->close();
    B->close();
}





bool
MAST::StructuralBucklingEigenproblemAssembly::
eigenproblem_sensitivity_assemble (const libMesh::ParameterVector& parameters,
                                   const unsigned int i,
                                   libMesh::SparseMatrix<Real>* sensitivity_A,
                                   libMesh::SparseMatrix<Real>* sensitivity_B)  {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());
    
    // zero the solution since it is not needed for eigenproblem
    eigen_sys.solution->zero();
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *sensitivity_A,
    &matrix_B = *sensitivity_B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dummy;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
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
        if (_base_sol.get()) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*_base_sol)(dof_indices[i]);
        }
        
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution(sol);
        
        
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
        eigen_sys.get_dof_map().constrain_element_matrix(A, dof_indices);
        eigen_sys.get_dof_map().constrain_element_matrix(B, dof_indices);
        
        // add to the global matrices
        matrix_A.add_matrix (A, dof_indices); // load independent
        matrix_B.add_matrix (B, dof_indices); // load dependent
    }
    
    // finalize the matrices for futher use.
    sensitivity_A->close();
    sensitivity_B->close();
    
    return true;
}



std::auto_ptr<MAST::ElementBase>
MAST::StructuralBucklingEigenproblemAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}



void
MAST::StructuralBucklingEigenproblemAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   RealMatrixX& mat_A) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealVectorX
    vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());

    mat_A.setZero();

    // calculate the Jacobian components
    e.internal_residual(true, vec, mat_A, false);
    e.prestress_residual(true, vec, mat_A);
    e.side_external_residual<Real>(true,
                                   vec,
                                   dummy,
                                   mat_A,
                                   _discipline->side_loads());
    e.volume_external_residual<Real>(true,
                                     vec,
                                     dummy,
                                     mat_A,
                                     _discipline->volume_loads());
}




void
MAST::StructuralBucklingEigenproblemAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealMatrixX& mat_A) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealVectorX
    vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    
    mat_A.setZero();
    
    // calculate the Jacobian components
    e.internal_residual_sensitivity(true, vec, mat_A, false);
    e.prestress_residual_sensitivity(true, vec, mat_A);
    e.side_external_residual_sensitivity<Real>(true,
                                               vec,
                                               dummy,
                                               mat_A,
                                               _discipline->side_loads());
    e.volume_external_residual_sensitivity<Real>(true,
                                                 vec,
                                                 dummy,
                                                 mat_A,
                                                 _discipline->volume_loads());
}




Real
MAST::StructuralBucklingEigenproblemAssembly::
critical_point_estimate_from_eigenproblem(Real v) const {

    Real rval = 0.;
    
    if (_use_linearized_formulation)
        rval = _lambda1 + v * (_lambda2 - _lambda1);
    else
        rval = _lambda1 + v/(1.+v) * (_lambda2 - _lambda1);
    
    return rval;
}

