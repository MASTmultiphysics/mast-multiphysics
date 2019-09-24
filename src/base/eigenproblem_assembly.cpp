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
#include "base/eigenproblem_assembly.h"
#include "base/nonlinear_system.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "base/eigenproblem_assembly_elem_operations.h"
#include "mesh/geom_elem.h"
#include "numerics/utility.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::EigenproblemAssembly::EigenproblemAssembly():
MAST::AssemblyBase(),
_base_sol                (nullptr),
_base_sol_sensitivity    (nullptr) {
    
}



MAST::EigenproblemAssembly::~EigenproblemAssembly() {
    
}



void
MAST::EigenproblemAssembly::set_base_solution(const libMesh::NumericVector<Real>& sol,
                                              bool if_sens) {
    
    if (!if_sens) {
        
        // make sure the pointer has been cleared
        libmesh_assert(!_base_sol);
        
        _base_sol             = &sol;
    }
    else {
        
        // make sure the pointer has been cleared
        libmesh_assert(!_base_sol_sensitivity);

        _base_sol_sensitivity = &sol;
    }
}




void
MAST::EigenproblemAssembly::clear_base_solution(bool if_sens) {
    
    if (!if_sens)
        _base_sol             = nullptr;
    else
        _base_sol_sensitivity = nullptr;
}




void
MAST::EigenproblemAssembly::
eigenproblem_assemble(libMesh::SparseMatrix<Real>* A,
                      libMesh::SparseMatrix<Real>* B) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

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
        localized_solution.reset(build_localized_vector(eigen_sys,
                                                         *_base_sol).release());

    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    RealMatrixX mat_A, mat_B;
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
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        ops.elem_calculations(mat_A, mat_B);
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);

        // constrain the element matrices.
        dof_map.constrain_element_matrix(A, dof_indices);
        dof_map.constrain_element_matrix(B, dof_indices);
        
        matrix_A.add_matrix (A, dof_indices); // load independent
        matrix_B.add_matrix (B, dof_indices); // load dependent
    }
    
    // finalize the data structures
    A->close();
    B->close();
}





bool
MAST::EigenproblemAssembly::
eigenproblem_sensitivity_assemble(const MAST::FunctionBase& f,
                                  libMesh::SparseMatrix<Real>* sensitivity_A,
                                  libMesh::SparseMatrix<Real>* sensitivity_B) {

    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());
    
    libMesh::SparseMatrix<Real>&  matrix_A = *sensitivity_A;
    libMesh::SparseMatrix<Real>&  matrix_B = *sensitivity_B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // build localized solutions if needed
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    
    if (_base_sol) {
        
        localized_solution.reset(build_localized_vector(eigen_sys,
                                                         *_base_sol).release());
        
        // make sure that the sensitivity was also provided
        libmesh_assert(_base_sol_sensitivity);
        localized_solution_sens.reset(build_localized_vector(eigen_sys,
                                                              *_base_sol_sensitivity).release());
    }
    

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    RealMatrixX mat_A, mat_B;
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
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, tell the element about it
        if (_base_sol) {
            
            // set the element's base solution
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        
        // set the element's base solution sensitivity
        if (_base_sol) {
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution_sens)(dof_indices[i]);
        }
        
        ops.set_elem_solution_sensitivity(sol);
        ops.elem_sensitivity_calculations(f,
                                          _base_sol!=nullptr,
                                          mat_A,
                                          mat_B);
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);
        
        // constrain the element matrices.
        dof_map.constrain_element_matrix(A, dof_indices);
        dof_map.constrain_element_matrix(B, dof_indices);
        
        matrix_A.add_matrix (A, dof_indices);
        matrix_B.add_matrix (B, dof_indices);
    }
    
    // finalize the data structures
    sensitivity_A->close();
    sensitivity_B->close();
    
    return true;
}



