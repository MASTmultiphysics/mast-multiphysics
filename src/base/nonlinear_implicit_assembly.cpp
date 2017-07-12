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
#include "base/nonlinear_implicit_assembly.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::NonlinearImplicitAssembly::
NonlinearImplicitAssembly():
MAST::AssemblyBase(),
_post_assembly(nullptr) {
    
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
    
    _system->system().nonlinear_solver->residual_and_jacobian_object = this;
}



void
MAST::NonlinearImplicitAssembly::reattach_to_system() {

    libmesh_assert(_system);
    
    _system->system().nonlinear_solver->residual_and_jacobian_object = this;
}



void
MAST::NonlinearImplicitAssembly::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {

        _system->system().nonlinear_solver->residual_and_jacobian_object = nullptr;
    }
    
    _discipline    = nullptr;
    _system        = nullptr;
    _post_assembly = nullptr;
}


void
MAST::NonlinearImplicitAssembly::
set_post_assembly_operation(MAST::NonlinearImplicitAssembly::PostAssemblyOperation& post) {
    
    _post_assembly = &post;
}



namespace MAST {

    bool
    is_numerical_zero(const Real v, const Real eps) {
        
        return fabs(v) <= eps;
    }
    
    
    bool
    compare(const Real v1, const Real v2, const Real tol) {
        
        const Real
        eps      = 1.0e-7;
        
        bool rval = false;
        
        // check to see if the values are both small enough
        // to be zero
        if (MAST::is_numerical_zero(v1, eps) &&
            MAST::is_numerical_zero(v2, eps))
            rval = true;
        // check to see if the absolute difference is small enough
        else if (MAST::is_numerical_zero(v1-v2, eps))
            rval = true;
        // check to see if the relative difference is small enough
        else if (fabs(v1) > 0)
            rval = fabs((v1-v2)/v1) <= tol;
        
        return rval;
    }
    
    bool
    compare_matrix(const RealMatrixX& m0, const RealMatrixX& m, const Real tol) {
        
        unsigned int
        m0_rows = (unsigned int) m0.rows(),
        m0_cols = (unsigned int) m0.cols();
        libmesh_assert_equal_to(m0_rows,  m.rows());
        libmesh_assert_equal_to(m0_cols,  m.cols());
        
        
        bool pass = true;
        for (unsigned int i=0; i<m0_rows; i++) {
            for (unsigned int j=0; j<m0_cols; j++)
                if (!MAST::compare(m0(i,j), m(i,j), tol)) {
                    libMesh::out << "Failed comparison at (i,j) = ("
                    << i << ", " << j << ") : "
                    << "expected: " << m0(i,j) << "  , "
                    << "computed: " << m(i,j) << " : "
                    << "diff: " << m0(i,j) - m(i,j) << " , "
                    << "tol: " << tol << std::endl;
                    pass = false;
                }
        }
        
        return pass;
    }
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
    Real delta = 1.0e-8;
    
    for (unsigned int i=0; i<sol.size(); i++) {
        dsol = sol;
        dsol(i) += delta;
        
        e.set_solution(dsol);
        _elem_calculations(e, false, dres, dummy);
        jac.col(i) = (dres-res0)/delta;
    }
    
    // write the numerical and analytical jacobians
    libMesh::out
    << "Analytical Jacobian: " << std::endl
    << jac0
    << std::endl << std::endl
    << "Numerical Jacobian: " << std::endl
    << jac
    << std::endl << std::endl;
    
    MAST::compare_matrix(jac, jac0, 1.0e-5);
    // set the original solution vector for the element
    e.set_solution(sol);
}



void
MAST::NonlinearImplicitAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
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
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
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
                           J!=nullptr?true:false,
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
    
    if (R) R->close();
    if (J) J->close();
}




void
MAST::NonlinearImplicitAssembly::
linearized_jacobian_solution_product (const libMesh::NumericVector<Real>& X,
                                      const libMesh::NumericVector<Real>& dX,
                                      libMesh::NumericVector<Real>& JdX,
                                      libMesh::NonlinearImplicitSystem& S) {
    
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
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_perturbed_solution;
    
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     X).release());
    localized_perturbed_solution.reset(_build_localized_vector(nonlin_sys,
                                                               dX).release());
    
    
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
        dsol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol (i) = (*localized_solution)          (dof_indices[i]);
            dsol(i) = (*localized_perturbed_solution)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_perturbed_solution(dsol);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        _elem_linearized_jacobian_solution_product(*physics_elem,
                                                   vec);
        
        physics_elem->detach_active_solution_function();
        
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
_elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                           RealVectorX& vec) {
    
    // must be implemented in inherited class
    libmesh_error();
}



bool
MAST::NonlinearImplicitAssembly::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
        _sol_function->init( *nonlin_sys.solution);
    
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
        
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution(sol);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        _elem_sensitivity_calculations(*physics_elem, false, vec, mat);
        
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


