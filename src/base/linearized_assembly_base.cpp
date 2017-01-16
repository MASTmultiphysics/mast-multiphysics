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
#include "base/linearized_assembly_base.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"
#include "base/mesh_field_function.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::LinearizedAssemblyBase::
LinearizedAssemblyBase():
MAST::NonlinearImplicitAssembly(),
_base_sol(nullptr),
_base_sol_sensitivity(nullptr) {
    
}



MAST::LinearizedAssemblyBase::~LinearizedAssemblyBase() {
    
}




void
MAST::LinearizedAssemblyBase::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    MAST::NonlinearImplicitAssembly::attach_discipline_and_system(discipline,
                                                                  system);
    
    _base_sol             = nullptr;
    _base_sol_sensitivity = nullptr;
   
}





void
MAST::LinearizedAssemblyBase::
clear_discipline_and_system( ) {
    
    MAST::NonlinearImplicitAssembly::clear_discipline_and_system();
    
    _base_sol             = nullptr;
    _base_sol_sensitivity = nullptr;
}




void
MAST::LinearizedAssemblyBase::set_base_solution(const libMesh::NumericVector<Real>& sol,
                                             bool if_sens) {
    
    if (!if_sens) {
        
        // make sure that the pointer has been cleared
        libmesh_assert(!_base_sol);
        
        _base_sol             = &sol;
    }
    else {
        
        // make sure that the pointer has been cleared
        libmesh_assert(!_base_sol_sensitivity);
        
        _base_sol_sensitivity = &sol;
    }
}




void
MAST::LinearizedAssemblyBase::clear_base_solution(bool if_sens) {
    
    if (!if_sens)
        _base_sol             = nullptr;
    else
        _base_sol_sensitivity = nullptr;
}




const libMesh::NumericVector<Real>&
MAST::LinearizedAssemblyBase::base_sol(bool if_sens) const {
    
    if (!if_sens)
        return *_base_sol;
    else
        return *_base_sol_sensitivity;
}





void
MAST::LinearizedAssemblyBase::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    if (R) R->zero();
    if (J) J->zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol, vec, delta_sol;
    RealMatrixX    mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_base_solution,
    localized_perturbed_solution;
    
    
    // localize the base solution, if it was provided
    if (_base_sol)
        localized_base_solution.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());
    
    
    // localize sol to real vector
    localized_perturbed_solution.reset(_build_localized_vector(nonlin_sys,
                                                               X).release());
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        delta_sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        // first set the velocity to be zero
        physics_elem->set_velocity(sol);
        
        // next, set the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        
        physics_elem->set_solution(sol);
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            delta_sol(i) = (*localized_perturbed_solution)(dof_indices[i]);
        
        physics_elem->set_perturbed_solution(delta_sol);
        
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        _elem_calculations(*physics_elem,
                           J!=nullptr?true:false,
                           vec, mat);
        
        physics_elem->detach_active_solution_function();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        DenseRealMatrix m;
        if (R) MAST::copy(v, vec);
        if (J) MAST::copy(m, mat);
        
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
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J) J->close();
}






bool
MAST::LinearizedAssemblyBase::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libmesh_error(); // not implemented. Call the blocked assembly instead.
    
    return false;
}
