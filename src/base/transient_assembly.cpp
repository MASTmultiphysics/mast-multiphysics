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
#include "base/transient_assembly.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly_elem_operations.h"
#include "solver/transient_solver_base.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::TransientAssembly::TransientAssembly():
MAST::AssemblyBase       (),
_post_assembly           (nullptr) {

}



MAST::TransientAssembly::~TransientAssembly() {

}



void
MAST::TransientAssembly::
set_post_assembly_operation(MAST::TransientAssembly::PostAssemblyOperation& post) {
    
    _post_assembly = &post;
}



void
MAST::TransientAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    
    MAST::TransientSolverBase
    &solver = dynamic_cast<MAST::TransientSolverBase&>(*_elem_ops);
    MAST::NonlinearSystem&
    transient_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &transient_sys);
    
    if (R) R->zero();
    if (J) J->zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    
    
    
    // stores the localized solution, velocity, acceleration, etc. vectors.
    // These pointers will have to be deleted
    std::vector<libMesh::NumericVector<Real>*>
    local_qtys;
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    // ask the solver to localize the relevant solutions
    solver.build_local_quantities(X, local_qtys);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        solver.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        solver.init(geom_elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        solver.set_element_data(dof_indices, local_qtys);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        solver.elem_calculations(J!=nullptr?true:false,
                                             vec, mat);
        solver.clear_elem();

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

    // delete pointers to the local solutions
    for (unsigned int i=0; i<local_qtys.size(); i++)
        delete local_qtys[i];
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J && close_matrix) J->close();
}




void
MAST::TransientAssembly::
linearized_jacobian_solution_product (const libMesh::NumericVector<Real>& X,
                                      const libMesh::NumericVector<Real>& dX,
                                      libMesh::NumericVector<Real>& JdX,
                                      libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::TransientSolverBase
    &solver = dynamic_cast<MAST::TransientSolverBase&>(*_elem_ops);
    MAST::NonlinearSystem
    &transient_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &transient_sys);
    
    JdX.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    
    
    // stores the localized solution, velocity, acceleration, etc. vectors.
    // These pointers will have to be deleted
    std::vector<libMesh::NumericVector<Real>*>
    local_qtys,
    local_perturbed_qtys;
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);

    // ask the solver to localize the relevant solutions
    solver.build_local_quantities(X, local_qtys);
    solver.build_perturbed_local_quantities(dX, local_perturbed_qtys);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);

        MAST::GeomElem geom_elem;
        _elem_ops->set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        _elem_ops->init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        vec.setZero(ndofs);

        
        solver.set_element_data(dof_indices, local_qtys);
        solver.set_element_perturbed_data(dof_indices, local_perturbed_qtys);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);

        // perform the element level calculations
        solver.elem_linearized_jacobian_solution_product(vec);
        
        solver.clear_elem();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        JdX.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    // delete pointers to the local solutions
    for (unsigned int i=0; i<local_qtys.size(); i++) {
        
        delete local_qtys[i];
        delete local_perturbed_qtys[i];
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();

    JdX.close();
}



bool
MAST::TransientAssembly::
sensitivity_assemble (const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::TransientSolverBase
    &solver = dynamic_cast<MAST::TransientSolverBase&>(*_elem_ops);
    MAST::NonlinearSystem
    &nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, vec2;
    std::vector<RealVectorX> prev_local_sols;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::vector<libMesh::NumericVector<Real>*>
    local_qtys,
    prev_local_qtys;
    solver.build_local_quantities(*nonlin_sys.solution, local_qtys);
    solver.build_sensitivity_local_quantities(1, prev_local_qtys);

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( *nonlin_sys.solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        _elem_ops->set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        _elem_ops->init(geom_elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        vec.setZero(ndofs);
        vec2.setZero(ndofs);
        
        solver.set_element_data(dof_indices, local_qtys);
        solver.extract_element_sensitivity_data(dof_indices, prev_local_qtys, prev_local_sols);
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        solver.elem_sensitivity_contribution_previous_timestep(prev_local_sols, vec2);
        solver.elem_sensitivity_calculations(f, vec);
        
        vec += vec2;
        
        //        physics_elem->detach_active_solution_function();
        solver.clear_elem();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    for (unsigned int i=0; i<local_qtys.size(); i++) {
        
        delete local_qtys[i];
        delete prev_local_qtys[i];
    }

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();
    
    return true;
}


