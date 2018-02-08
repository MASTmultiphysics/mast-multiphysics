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
#include "elasticity/structural_buckling_eigenproblem_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_buckling_eigenproblem_elem_operations.h"
#include "property_cards/element_property_card_1D.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "mesh/local_elem_fe.h"
#include "base/assembly_base.h"
#include "base/nonlinear_system.h"
#include "base/parameter.h"
#include "numerics/utility.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"


MAST::StructuralBucklingEigenproblemAssembly::
StructuralBucklingEigenproblemAssembly():
MAST::EigenproblemAssembly(),
_use_linearized_formulation(true),
_load_param(nullptr),
_lambda1(0.),
_lambda2(0.),
_sol1(nullptr),
_sol2(nullptr) {
    
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
    _load_param                 = nullptr;
    _lambda1                    = 0.;
    _lambda2                    = 0.;
    _sol1                       = nullptr;
    _sol2                       = nullptr;
    
    MAST::EigenproblemAssembly::clear_discipline_and_system();
}



void
MAST::StructuralBucklingEigenproblemAssembly::
eigenproblem_assemble(libMesh::SparseMatrix<Real> *A,
                      libMesh::SparseMatrix<Real> *B)  {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());

    // create the localized solution vectors
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution1,
    localized_solution2;
    
    localized_solution1.reset(build_localized_vector(eigen_sys,
                                                      *_sol1).release());
    localized_solution2.reset(build_localized_vector(eigen_sys,
                                                      *_sol2).release());
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *A,
    &matrix_B = *B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    RealMatrixX dummy, mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    eigen_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    eigen_sys.get_mesh().active_local_elements_end();
    
    MAST::EigenproblemAssemblyElemOperations
    &ops = dynamic_cast<MAST::EigenproblemAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dummy.setZero(ndofs, ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        ////////////////////////////////////////////////////////////////
        // first calculate the tangent stiffness matrix about the
        // first load factor and solution
        ////////////////////////////////////////////////////////////////
        (*_load_param) = _lambda1;
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution1)(dof_indices[i]);
    
        ops.init(*elem);
        ops.set_elem_solution(sol);
        ops.elem_calculations(mat_A, dummy);

        DenseRealMatrix AA, BB;
        MAST::copy(AA, mat_A); // copy to the libMesh matrix for further processing
        dof_map.constrain_element_matrix(AA, dof_indices); // constrain the element matrices.
        matrix_A.add_matrix (AA, dof_indices); // add to the global matrices
        
        
        ////////////////////////////////////////////////////////////////
        // next, calculate the tangent stiffness matrix about the
        // second load factor and solution
        ////////////////////////////////////////////////////////////////
        (*_load_param) = _lambda2;
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution2)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        ops.elem_calculations(mat_B, dummy);
        mat_B  *= -1.;
        if (_use_linearized_formulation)
            mat_B  += mat_A;
        
        MAST::copy(BB, mat_B); // copy to the libMesh matrix for further processing
        dof_map.constrain_element_matrix(BB, dof_indices); // constrain the element matrices.
        matrix_B.add_matrix (BB, dof_indices); // add to the global matrices
    }
    
    // finalize the matrices for futher use.
    A->close();
    B->close();
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




