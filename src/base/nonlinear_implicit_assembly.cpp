/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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



MAST::NonlinearImplicitAssembly::
NonlinearImplicitAssembly():
MAST::AssemblyBase() {
    
}



MAST::NonlinearImplicitAssembly::~NonlinearImplicitAssembly() {
    
}




void
MAST::NonlinearImplicitAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(system.system());
    
    nonlin_sys.nonlinear_solver->residual_and_jacobian_object = this;
    nonlin_sys.attach_sensitivity_assemble_object(*this);
}




void
MAST::NonlinearImplicitAssembly::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {

        libMesh::NonlinearImplicitSystem& nonlin_sys =
        dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
        
        nonlin_sys.nonlinear_solver->residual_and_jacobian_object = NULL;
        nonlin_sys.reset_sensitivity_assembly();
    }
    
    _discipline = NULL;
    _system     = NULL;
}



void
MAST::NonlinearImplicitAssembly::
_check_element_numerical_jacobian(MAST::ElementBase& e,
                                  RealVectorX& sol) {
    RealVectorX
    dsol,
    res0,
    dres;

    RealMatrixX
    jac0,
    jac,
    dummy;
    
    unsigned int ndofs = (unsigned int)sol.size();
    res0.setZero(ndofs);
    dres.setZero(ndofs);
    jac0.setZero(ndofs, ndofs);
    jac.setZero(ndofs, ndofs);

    e.set_solution(sol);
    _elem_calculations(e, true, res0, jac0);
    Real delta = 1.0e-5;
    
    for (unsigned int i=0; i<sol.size(); i++) {
        dsol = sol;
        dsol(i) += delta;
        
        e.set_solution(dsol);
        _elem_calculations(e, false, dres, dummy);
        jac.col(i) = (dres-res0)/delta;
    }
    
    // write the numerical and analytical jacobians
    std::cout
    << "Analytical Jacobian: " << std::endl
    << jac0
    << std::endl << std::endl
    << "Numerical Jacobian: " << std::endl
    << jac
    << std::endl << std::endl;
    
}



void
MAST::NonlinearImplicitAssembly::
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
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init_for_system_and_solution(*_system, X);
    
    
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
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        physics_elem->set_solution(sol);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        _elem_calculations(*physics_elem,
                           J!=NULL?true:false,
                           vec, mat);
        
        physics_elem->detach_active_solution_function();

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
            nonlin_sys.get_dof_map().constrain_element_matrix_and_vector(m, v, dof_indices);
        else if (R)
            nonlin_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        else
            nonlin_sys.get_dof_map().constrain_element_matrix(m, dof_indices);
        
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
MAST::NonlinearImplicitAssembly::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     *nonlin_sys.solution).release());
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init_for_system_and_solution(*_system, *nonlin_sys.solution);
    
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
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        physics_elem->sensitivity_param = _discipline->get_parameter(parameters[i]);
        physics_elem->set_solution(sol);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        _elem_sensitivity_calculations(*physics_elem, vec);
        
        // the sensitivity method provides sensitivity of the residual.
        // Hence, this is multiplied with -1 to make it the RHS of the
        // sensitivity equations. 
        vec *= -1.;
        
        physics_elem->detach_active_solution_function();
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        nonlin_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
    }
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();

    return true;
}


