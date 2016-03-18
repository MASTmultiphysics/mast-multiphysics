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
#include "base/eigenproblem_assembly.h"
#include "base/nonlinear_system.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::EigenproblemAssembly::EigenproblemAssembly():
MAST::AssemblyBase(),
_base_sol(NULL),
_base_sol_sensitivity(NULL) {
    
}



MAST::EigenproblemAssembly::~EigenproblemAssembly() {
    
}




void
MAST::EigenproblemAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
    
    // now attach this to the system
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(system.system());
    
    eigen_sys.attach_eigenproblem_assemble_object(*this);
}




void
MAST::EigenproblemAssembly::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {

        MAST::NonlinearSystem& sys =
        dynamic_cast<MAST::NonlinearSystem&>(_system->system());
        
        sys.reset_eigenproblem_assemble_object();
    }
    
    _discipline           = NULL;
    _system               = NULL;
    _base_sol             = NULL;
    _base_sol_sensitivity = NULL;
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
        _base_sol             = NULL;
    else
        _base_sol_sensitivity = NULL;
}




void
MAST::EigenproblemAssembly::
eigenproblem_assemble(libMesh::SparseMatrix<Real>* A,
                      libMesh::SparseMatrix<Real>* B) {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());
    
    libMesh::SparseMatrix<Real>
    &matrix_A = *A,
    &matrix_B = *B;
    
    matrix_A.zero();
    matrix_B.zero();
    

    // build localized solutions if needed
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution;
    
    if (_base_sol)
        localized_solution.reset(_build_localized_vector(eigen_sys,
                                                         *_base_sol).release());

    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
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
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        
        _elem_calculations(*physics_elem, mat_A, mat_B);

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
    
    
}





bool
MAST::EigenproblemAssembly::
eigenproblem_sensitivity_assemble(const libMesh::ParameterVector& parameters,
                                  const unsigned int i,
                                  libMesh::SparseMatrix<Real>* sensitivity_A,
                                  libMesh::SparseMatrix<Real>* sensitivity_B) {
    
    MAST::NonlinearSystem& eigen_sys =
    dynamic_cast<MAST::NonlinearSystem&>(_system->system());

    // zero the solution since it is not needed for eigenproblem
    eigen_sys.solution->zero();
    
    libMesh::SparseMatrix<Real>&  matrix_A = *sensitivity_A;
    libMesh::SparseMatrix<Real>&  matrix_B = *sensitivity_B;
    
    matrix_A.zero();
    matrix_B.zero();
    
    // build localized solutions if needed
    std::auto_ptr<libMesh::NumericVector<Real> >
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
    RealVectorX sol;
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
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, tell the element about it
        if (_base_sol) {
            
            // set the element's base solution
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        
        // set the element's base solution sensitivity
        if (_base_sol) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution_sens)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol, true);
        
        // tell the element about the sensitivity parameter
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        
        _elem_sensitivity_calculations(*physics_elem, mat_A, mat_B);

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
    
    return true;
}

