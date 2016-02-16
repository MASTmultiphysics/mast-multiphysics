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
#include "base/complex_assembly_base.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"
#include "base/mesh_field_function.h"
#include "solver/complex_solver_base.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::ComplexAssemblyBase::
ComplexAssemblyBase():
MAST::AssemblyBase(),
_if_assemble_real(true),
_complex_solver(NULL),
_base_sol(NULL),
_base_sol_sensitivity(NULL) {
    
}



MAST::ComplexAssemblyBase::~ComplexAssemblyBase() {
    
}




void
MAST::ComplexAssemblyBase::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::ComplexSolverBase& solver,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _if_assemble_real     = true;
    _complex_solver       = &solver;
    _discipline           = &discipline;
    _system               = &system;
    _base_sol             = NULL;
    _base_sol_sensitivity = NULL;
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(system.system());
    
    _complex_solver->set_assembly(*this);
    
    nonlin_sys.nonlinear_solver->residual_and_jacobian_object = this;
    nonlin_sys.attach_sensitivity_assemble_object(*this);
}




void
MAST::ComplexAssemblyBase::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {
        
        libMesh::NonlinearImplicitSystem& nonlin_sys =
        dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
        
        nonlin_sys.nonlinear_solver->residual_and_jacobian_object = NULL;
        nonlin_sys.reset_sensitivity_assembly();
    }
    
    _complex_solver->clear_assembly();
    
    _if_assemble_real     = true;
    _complex_solver       = NULL;
    _discipline           = NULL;
    _system               = NULL;
    _base_sol             = NULL;
    _base_sol_sensitivity = NULL;
}




void
MAST::ComplexAssemblyBase::set_base_solution(libMesh::NumericVector<Real>& sol,
                                             bool if_sens) {
    
    if (!if_sens)
        _base_sol             = &sol;
    else
        _base_sol_sensitivity = &sol;
}



Real
MAST::ComplexAssemblyBase::residual_l2_norm() {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol,
    vec_re;
    RealMatrixX
    mat_re;
    
    ComplexVectorX
    delta_sol,
    vec;
    ComplexMatrixX
    mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    residual_re(nonlin_sys.solution->zero_clone().release()),
    residual_im(nonlin_sys.solution->zero_clone().release()),
    localized_base_solution,
    localized_real_solution(_build_localized_vector
                            (nonlin_sys,
                             _complex_solver->real_solution()).release()),
    localized_imag_solution(_build_localized_vector
                            (nonlin_sys,
                             _complex_solver->imag_solution()).release());
    
    
    if (_base_sol)
        localized_base_solution.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());

    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init_for_system_and_solution(*_system, X);
    
    
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
        
        // set zero velocity for base solution
        physics_elem->set_velocity(sol);
        
        // set the value of the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        physics_elem->set_solution(sol);
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            delta_sol(i) = Complex((*localized_real_solution)(dof_indices[i]),
                                   (*localized_imag_solution)(dof_indices[i]));
        
        physics_elem->set_complex_solution(delta_sol);
        
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        _elem_calculations(*physics_elem,
                           false,
                           vec, mat);
        
        physics_elem->detach_active_solution_function();
        
        // add to the real part of the residual
        vec_re  =  vec.real();
        DenseRealVector v;
        MAST::copy(v, vec_re);
        nonlin_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        residual_re->add_vector(v, dof_indices);
        
        // now add to the imaginary part of the residual
        vec_re  =  vec.imag();
        v.zero();
        MAST::copy(v, vec_re);
        nonlin_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        residual_im->add_vector(v, dof_indices);
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    residual_re->close();
    residual_im->close();
    
    // return the residual
    Real
    l2_real =  residual_re->l2_norm(),
    l2_imag =  residual_im->l2_norm();
    
    return sqrt(pow(l2_real,2) + pow(l2_imag,2));
}





void
MAST::ComplexAssemblyBase::
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
    RealVectorX    sol, vec_re;
    RealMatrixX    mat_re;
    ComplexVectorX delta_sol, vec;
    ComplexMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_base_solution,
    localized_real_solution,
    localized_imag_solution;

    
    // localize the base solution, if it was provided
    if (_base_sol)
        localized_base_solution.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());
    
    
    // use the provided solution as either the real or imaginary solution
    if (_if_assemble_real) {
        
        // localize sol to real vector
        localized_real_solution.reset(_build_localized_vector(nonlin_sys,
                                                              X).release());
        // localize sol to imag vector
        localized_imag_solution.reset(_build_localized_vector
                                      (nonlin_sys,
                                       _complex_solver->imag_solution()).release());
    }
    else {

        // localize sol to real vector
        localized_real_solution.reset(_build_localized_vector
                                      (nonlin_sys,
                                       _complex_solver->real_solution()).release());
        // localize sol to imag vector
        localized_imag_solution.reset(_build_localized_vector(nonlin_sys,
                                                              X).release());
    }
    
    
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
            delta_sol(i) = Complex((*localized_real_solution)(dof_indices[i]),
                                   (*localized_imag_solution)(dof_indices[i]));
        
        physics_elem->set_complex_solution(delta_sol);
        
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);

        
        // perform the element level calculations
        _elem_calculations(*physics_elem,
                           J!=NULL?true:false,
                           vec, mat);
        
        physics_elem->detach_active_solution_function();
        
        // extract the real or the imaginary part of the matrix/vector
        //  The complex system of equations
        //     (J_R + i J_I) (x_R + i x_I) - (f_R + i f_I) = 0
        //  is rewritten as
        //     [ J_R   -J_I] {x_R}  -  {f_R}  = {0}
        //     [ J_I    J_R] {x_I}  -  {f_I}  = {0}
        //
        if (_if_assemble_real) {
            
            if (R) vec_re  =  vec.real();
            if (J) mat_re  =  mat.real();
        }
        else {
            
            // the imaginary residual also uses the real part of the Jacobian
            if (R) vec_re  =  vec.imag();
            if (J) mat_re  =  mat.real();
        }
        
        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        DenseRealMatrix m;
        if (R) MAST::copy(v, vec_re);
        if (J) MAST::copy(m, mat_re);
        
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
MAST::ComplexAssemblyBase::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol, vec_re;
    ComplexVectorX vec, delta_sol;
    ComplexMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_base_solution,
    localized_base_solution_sens,
    localized_real_solution,
    localized_imag_solution;
    
    
    // localize the base solution, if it was provided
    if (_base_sol)
        localized_base_solution.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());
    if (_base_sol_sensitivity)
        localized_base_solution.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol_sensitivity).release());
    
    
    // localize sol to real vector
    localized_real_solution.reset(_build_localized_vector
                                  (nonlin_sys,
                                   _complex_solver->real_solution()).release());
    // localize sol to imag vector
    localized_imag_solution.reset(_build_localized_vector
                                  (nonlin_sys,
                                   _complex_solver->imag_solution()).release());

    
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
        delta_sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        // set the value of the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        physics_elem->set_solution(sol);
        
        // set value of the base solution sensitivity
        if (_base_sol_sensitivity)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution_sens)(dof_indices[i]);
        physics_elem->set_solution(sol, true);
        
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            delta_sol(i) = Complex((*localized_real_solution)(dof_indices[i]),
                                   (*localized_imag_solution)(dof_indices[i]));
        
        physics_elem->set_complex_solution(delta_sol);
        
        physics_elem->sensitivity_param = _discipline->get_parameter(&(parameters[i].get()));
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        _elem_sensitivity_calculations(*physics_elem, false, vec, mat);
        
        // the sensitivity method provides sensitivity of the residual.
        // Hence, this is multiplied with -1 to make it the RHS of the
        // sensitivity equations.
        vec *= -1.;
        
        physics_elem->detach_active_solution_function();
        
        // extract the real or the imaginary part of the matrix/vector
        if (_if_assemble_real)
            vec_re  =  vec.real();
        else
            vec_re  =  vec.imag();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec_re);
        
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
