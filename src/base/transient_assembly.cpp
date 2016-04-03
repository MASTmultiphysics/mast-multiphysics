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
#include "base/transient_assembly.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "solver/transient_solver_base.h"
#include "numerics/utility.h"
#include "base/mesh_field_function.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::TransientAssembly::
TransientAssembly():
MAST::AssemblyBase(),
_transient_solver(NULL) {

}



MAST::TransientAssembly::~TransientAssembly() {

}





void
MAST::TransientAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase& discipline,
                             MAST::TransientSolverBase& solver,
                             MAST::SystemInitialization& sys) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline        = &discipline;
    _transient_solver  = &solver;
    _system            = &sys;
    
    // now attach this to the system
    libMesh::NonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());

    // attach this assembly object to the transient solver
    solver.set_assembly(*this);
    
    transient_sys.nonlinear_solver->residual_and_jacobian_object = this;
    transient_sys.attach_sensitivity_assemble_object(*this);
}



void
MAST::TransientAssembly::reattach_to_system() {
    
    libmesh_assert(_system);
    
    // now attach this to the system
    libMesh::NonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    transient_sys.nonlinear_solver->residual_and_jacobian_object = this;
    transient_sys.attach_sensitivity_assemble_object(*this);
}



void
MAST::TransientAssembly::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {

        libMesh::NonlinearImplicitSystem& transient_sys =
        dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
        
        transient_sys.nonlinear_solver->residual_and_jacobian_object = NULL;
        transient_sys.reset_sensitivity_assembly();
    }
    
    // clear the association of this assembly object to the solver
    _transient_solver->clear_assembly();
    
    _discipline       = NULL;
    _transient_solver = NULL;
    _system           = NULL;
}



void
MAST::TransientAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libMesh::NonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
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
    std::auto_ptr<MAST::ElementBase> physics_elem;
    

    // stores the localized solution, velocity, acceleration, etc. vectors.
    // These pointers will have to be deleted
    std::vector<libMesh::NumericVector<Real>*>
    local_qtys;
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init_for_system_and_solution(*_system, X);

    // ask the solver to localize the relevant solutions
    _transient_solver->build_local_quantities(X, local_qtys);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        _transient_solver->_set_element_data(dof_indices,
                                             local_qtys,
                                             *physics_elem);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);

        // perform the element level calculations
        _transient_solver->_elem_calculations(*physics_elem,
                                              dof_indices,
                                              J!=NULL?true:false,
                                              vec, mat);
        
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
    
    // delete pointers to the local solutions
    for (unsigned int i=0; i<local_qtys.size(); i++)
        delete local_qtys[i];
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();

    if (R) R->close();
    if (J) J->close();
}



bool
MAST::TransientAssembly::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {

    libMesh::NonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, vel;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    // localize the solution and velocity for element assembly
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution(_build_localized_vector(transient_sys,
                                               _transient_solver->solution(0)).release()),
    localized_velocity(_build_localized_vector(transient_sys,
                                               _transient_solver->velocity(0)).release());
    
    
    // ask the solver to provide the velocity estimate
    const libMesh::NumericVector<Real>
    &solution = *localized_solution,
    &velocity = *localized_velocity;

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init_for_system_and_solution(*_system, solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vel.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol(i) = solution(dof_indices[i]);
            vel(i) = velocity(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(vel);

        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);

        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution(sol);
        
        // perform the element level calculations
        _transient_solver->_elem_sensitivity_calculations(*physics_elem, dof_indices, vec);
        
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
    }
    
    
    sensitivity_rhs.close();
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();

    return true;
}


