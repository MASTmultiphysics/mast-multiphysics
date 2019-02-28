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
#include "base/nonlinear_implicit_assembly.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "numerics/utility.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::NonlinearImplicitAssembly::
NonlinearImplicitAssembly():MAST::AssemblyBase(),
_post_assembly           (nullptr),
_res_l2_norm             (0.),
_first_iter_res_l2_norm  (-1.) {
    
}



MAST::NonlinearImplicitAssembly::~NonlinearImplicitAssembly() {
    
}



void
MAST::NonlinearImplicitAssembly::
set_post_assembly_operation(MAST::NonlinearImplicitAssembly::PostAssemblyOperation& post) {
    
    _post_assembly = &post;
}



void
MAST::NonlinearImplicitAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    if (R) R->zero();
    if (J) J->zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        ops.init(*elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        ops.elem_calculations(J!=nullptr?true:false,
                                              vec, mat);
        
//        physics_elem->detach_active_solution_function();

        ops.clear_elem();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        DenseRealMatrix m;
        if (R)
            MAST::copy(v, vec);
        if (J)
            MAST::copy(m, mat);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        if (R && J)
            dof_map.constrain_element_matrix_and_vector(m, v, dof_indices);
        else if (R)
            dof_map.constrain_element_vector(v, dof_indices);
        else
            dof_map.constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        if (R) R->add_vector(v, dof_indices);
        if (J) J->add_matrix(m, dof_indices);
    }
    
    // call the post assembly object, if provided by user
    if (_post_assembly)
        _post_assembly->post_assembly(X, R, J, S);
    

    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) {
        
        R->close();
        _res_l2_norm = R->l2_norm();
        if (_first_iter_res_l2_norm < 0.)
            _first_iter_res_l2_norm = _res_l2_norm;
    }
    if (J) J->close();
}




void
MAST::NonlinearImplicitAssembly::
linearized_jacobian_solution_product (const libMesh::NumericVector<Real>& X,
                                      const libMesh::NumericVector<Real>& dX,
                                      libMesh::NumericVector<Real>& JdX,
                                      libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    // zero the solution vector
    JdX.zero();
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, dsol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_perturbed_solution;
    
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     X).release());
    localized_perturbed_solution.reset(build_localized_vector(nonlin_sys,
                                                               dX).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        ops.init(*elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dsol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol (i) = (*localized_solution)          (dof_indices[i]);
            dsol(i) = (*localized_perturbed_solution)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        ops.set_elem_perturbed_solution(dsol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        ops.elem_linearized_jacobian_solution_product(vec);
        
        //physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        JdX.add_vector(v, dof_indices);
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    JdX.close();
}



void
MAST::NonlinearImplicitAssembly::
second_derivative_dot_solution_assembly (const libMesh::NumericVector<Real>& X,
                                         const libMesh::NumericVector<Real>& dX,
                                         libMesh::SparseMatrix<Real>& d_JdX_dX,
                                         libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    // zero the matrix
    d_JdX_dX.zero();
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dsol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_perturbed_solution;
    
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     X).release());
    localized_perturbed_solution.reset(build_localized_vector(nonlin_sys,
                                                               dX).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        ops.init(*elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dsol.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol (i) = (*localized_solution)          (dof_indices[i]);
            dsol(i) = (*localized_perturbed_solution)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        ops.set_elem_solution_sensitivity(dsol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        ops.elem_second_derivative_dot_solution_assembly(mat);
        
//        physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealMatrix m;
        MAST::copy(m, mat);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        d_JdX_dX.add_matrix(m, dof_indices);
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    d_JdX_dX.close();
}





bool
MAST::NonlinearImplicitAssembly::
sensitivity_assemble (const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     *nonlin_sys.solution).release());
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( *nonlin_sys.solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        ops.init(*elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        ops.elem_sensitivity_calculations(f, vec);
        
//        physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
    }
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();
    
    return true;
}


