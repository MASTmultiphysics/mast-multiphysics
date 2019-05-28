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
#include "base/complex_assembly_base.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/parameter.h"
#include "base/nonlinear_system.h"
#include "base/complex_assembly_elem_operations.h"
#include "mesh/geom_elem.h"
#include "numerics/utility.h"
#include "solver/complex_solver_base.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"



MAST::ComplexAssemblyBase::ComplexAssemblyBase():
MAST::AssemblyBase(),
_base_sol                 (nullptr),
_base_sol_sensitivity     (nullptr) {
    
}



MAST::ComplexAssemblyBase::~ComplexAssemblyBase() {
    
}





void
MAST::ComplexAssemblyBase::set_base_solution(const libMesh::NumericVector<Real>& sol,
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
MAST::ComplexAssemblyBase::clear_base_solution(bool if_sens) {
    
    if (!if_sens)
        _base_sol             = nullptr;
    else
        _base_sol_sensitivity = nullptr;
}




const libMesh::NumericVector<Real>&
MAST::ComplexAssemblyBase::base_sol(bool if_sens) const {
    
    if (!if_sens)
        return *_base_sol;
    else
        return *_base_sol_sensitivity;
}





Real
MAST::ComplexAssemblyBase::residual_l2_norm(const libMesh::NumericVector<Real>& real,
                                            const libMesh::NumericVector<Real>& imag) {
    
    START_LOG("complex_solve()", "Residual-L2");

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    residual_re(nonlin_sys.solution->zero_clone().release()),
    residual_im(nonlin_sys.solution->zero_clone().release()),
    localized_base_solution,
    localized_real_solution(build_localized_vector(nonlin_sys, real).release()),
    localized_imag_solution(build_localized_vector(nonlin_sys, imag).release());
    
    
    if (_base_sol)
        localized_base_solution.reset(build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());

    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::ComplexAssemblyElemOperations&
    ops = dynamic_cast<MAST::ComplexAssemblyElemOperations&>(*_elem_ops);
    
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
        delta_sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        // set zero velocity for base solution
        ops.set_elem_velocity(sol);
        
        // set the value of the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        ops.set_elem_solution(sol);
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            delta_sol(i) = Complex((*localized_real_solution)(dof_indices[i]),
                                   (*localized_imag_solution)(dof_indices[i]));
        
        ops.set_elem_complex_solution(delta_sol);
        
        
//        if (_sol_function)
//            ops.attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        ops.elem_calculations(false,
                                             vec, mat);
        ops.clear_elem();
        
//        ops.detach_active_solution_function();
        
        // add to the real part of the residual
        vec_re  =  vec.real();
        DenseRealVector v;
        MAST::copy(v, vec_re);
        dof_map.constrain_element_vector(v, dof_indices);
        residual_re->add_vector(v, dof_indices);
        
        // now add to the imaginary part of the residual
        vec_re  =  vec.imag();
        v.zero();
        MAST::copy(v, vec_re);
        dof_map.constrain_element_vector(v, dof_indices);
        residual_im->add_vector(v, dof_indices);
        dof_indices.clear();
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
    
    STOP_LOG("complex_solve()", "Residual-L2");

    return sqrt(pow(l2_real,2) + pow(l2_imag,2));
}



void
MAST::ComplexAssemblyBase::
residual_and_jacobian_field_split (const libMesh::NumericVector<Real>& X_R,
                                   const libMesh::NumericVector<Real>& X_I,
                                   libMesh::NumericVector<Real>& R_R,
                                   libMesh::NumericVector<Real>& R_I,
                                   libMesh::SparseMatrix<Real>&  J_R,
                                   libMesh::SparseMatrix<Real>&  J_I) {

    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    R_R.zero();
    R_I.zero();
    J_R.zero();
    J_I.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol, vec_re;
    RealMatrixX    mat_re;
    ComplexVectorX delta_sol, vec;
    ComplexMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_base_solution,
    localized_real_solution,
    localized_imag_solution;
    
    
    // localize the base solution, if it was provided
    if (_base_sol)
        localized_base_solution.reset(build_localized_vector(nonlin_sys,
                                                              *_base_sol).release());
    
    
    // localize sol to real vector
    localized_real_solution.reset(build_localized_vector(nonlin_sys,
                                                              X_R).release());
    // localize sol to imag vector
    localized_imag_solution.reset(build_localized_vector(nonlin_sys,
                                                              X_I).release());
    
    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::ComplexAssemblyElemOperations&
    ops = dynamic_cast<MAST::ComplexAssemblyElemOperations&>(*_elem_ops);

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
        delta_sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        // first set the velocity to be zero
        ops.set_elem_velocity(sol);
        
        // next, set the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++)
            delta_sol(i) = Complex((*localized_real_solution)(dof_indices[i]),
                                   (*localized_imag_solution)(dof_indices[i]));
        
        ops.set_elem_complex_solution(delta_sol);
        
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        ops.elem_calculations(true,
                                             vec,
                                             mat);
        ops.clear_elem();

        vec *= -1.;
        
        //physics_elem->detach_active_solution_function();
        
        // extract the real or the imaginary part of the matrix/vector
        //  The complex system of equations
        //     (J_R + i J_I) (x_R + i x_I) + (r_R + i r_I) = 0
        //  is rewritten as
        //     [ J_R   -J_I] {x_R}  +  {r_R}  = {0}
        //     [ J_I    J_R] {x_I}  +  {r_I}  = {0}
        //
        DenseRealVector v;
        DenseRealMatrix m;

        // copy the real part of the residual and Jacobian
        MAST::copy(v, vec.real());
        MAST::copy(m, mat.real());
        
        dof_map.constrain_element_matrix_and_vector(m, v, dof_indices);
        R_R.add_vector(v, dof_indices);
        J_R.add_matrix(m, dof_indices);

        
        // copy the imag part of the residual and Jacobian
        v.zero();
        m.zero();
        MAST::copy(v, vec.imag());
        MAST::copy(m, mat.imag());
        
        dof_map.constrain_element_matrix_and_vector(m, v, dof_indices);
        R_I.add_vector(v, dof_indices);
        J_I.add_matrix(m, dof_indices);
        dof_indices.clear();
    }
    
    
    // if a solution function is attached, clear it
    //if (_sol_function)
    //    _sol_function->clear();
    
    R_R.close();
    R_I.close();
    J_R.close();
    J_I.close();
    
    libMesh::out << "R_R: " << R_R.l2_norm() << "   R_I: " << R_I.l2_norm() << std::endl;
}







void
MAST::ComplexAssemblyBase::
residual_and_jacobian_blocked (const libMesh::NumericVector<Real>& X,
                               libMesh::NumericVector<Real>& R,
                               libMesh::SparseMatrix<Real>&  J,
                               MAST::Parameter* p) {

    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    START_LOG("residual_and_jacobian()", "ComplexSolve");
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    R.zero();
    J.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol;
    ComplexVectorX
    delta_sol,
    vec;
    ComplexMatrixX
    mat,
    dummy;

    // get the petsc vector and matrix objects
    Mat
    jac_bmat = dynamic_cast<libMesh::PetscMatrix<Real>&>(J).mat();
    
    PetscInt ierr;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    const std::vector<libMesh::dof_id_type>&
    send_list = nonlin_sys.get_dof_map().get_send_list();
    
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_base_solution,
    localized_complex_sol(libMesh::NumericVector<Real>::build(nonlin_sys.comm()).release());
    
    // prepare a send list for localization of the complex solution
    std::vector<libMesh::dof_id_type>
    complex_send_list(2*send_list.size());
    
    for (unsigned int i=0; i<send_list.size(); i++) {
        complex_send_list[2*i  ] = 2*send_list[i];
        complex_send_list[2*i+1] = 2*send_list[i]+1;
    }

    localized_complex_sol->init(2*nonlin_sys.n_dofs(),
                                2*nonlin_sys.n_local_dofs(),
                                complex_send_list,
                                false,
                                libMesh::GHOSTED);
    X.localize(*localized_complex_sol, complex_send_list);
    
    // localize the base solution, if it was provided
    if (_base_sol)
        localized_base_solution.reset(build_localized_vector(nonlin_sys,
                                                             *_base_sol).release());
    
    

    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::ComplexAssemblyElemOperations&
    ops = dynamic_cast<MAST::ComplexAssemblyElemOperations&>(*_elem_ops);

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
        delta_sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        // first set the velocity to be zero
        ops.set_elem_velocity(sol);
        
        // next, set the base solution, if provided
        if (_base_sol)
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_base_solution)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
        // set the value of the small-disturbance solution
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            
            // get the complex block for this dof
            delta_sol(i) = Complex((*localized_complex_sol)(2*dof_indices[i]),
                                   (*localized_complex_sol)(2*dof_indices[i]+1));
        }
        
        ops.set_elem_complex_solution(delta_sol);
        
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        ops.elem_calculations(true, vec, mat);
        
        // if sensitivity was requested, then ask the element for sensitivity
        // of the residual
        if (p) {
            
            // set the sensitivity of complex sol to zero
            delta_sol.setZero();
            ops.set_elem_complex_solution_sensitivity(delta_sol);
            vec.setZero();
            ops.elem_sensitivity_calculations(*p, vec);
        }
        
        ops.clear_elem();

        //physics_elem->detach_active_solution_function();
        
        // extract the real or the imaginary part of the matrix/vector
        //  The complex system of equations
        //     (J_R + i J_I) (x_R + i x_I) + (r_R + i r_I) = 0
        //  is rewritten as
        //     [ J_R   -J_I] {x_R}  +  {r_R}  = {0}
        //     [ J_I    J_R] {x_I}  +  {r_I}  = {0}
        //
        DenseRealVector v_R, v_I;
        DenseRealMatrix m_R, m_I1, m_I2;
        std::vector<Real> vals(4);
        
        // copy the real part of the residual and Jacobian
        MAST::copy( m_R, mat.real());
        MAST::copy(m_I1, mat.imag()); m_I1 *= -1.;   // this is the -J_I component
        MAST::copy(m_I2, mat.imag());                // this is the J_I component
        MAST::copy( v_R, vec.real());
        MAST::copy( v_I, vec.imag());
        dof_map.constrain_element_matrix(m_R,  dof_indices);
        dof_map.constrain_element_matrix(m_I1, dof_indices);
        dof_map.constrain_element_matrix(m_I2, dof_indices);
        dof_map.constrain_element_vector(v_R,  dof_indices);
        dof_map.constrain_element_vector(v_I,  dof_indices);
        
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            
            R.add(2*dof_indices[i],     v_R(i));
            R.add(2*dof_indices[i]+1,   v_I(i));
            
            for (unsigned int j=0; j<dof_indices.size(); j++) {
                vals[0] = m_R (i,j);
                vals[1] = m_I1(i,j);
                vals[2] = m_I2(i,j);
                vals[3] = m_R (i,j);
                ierr = MatSetValuesBlocked(jac_bmat,
                                           1, (PetscInt*)&dof_indices[i],
                                           1, (PetscInt*)&dof_indices[j],
                                           &vals[0],
                                           ADD_VALUES);
            }
        }
    }
    
    
    // if a solution function is attached, clear it
    //if (_sol_function)
    //    _sol_function->clear();
    
    R.close();
    J.close();
    
    libMesh::out << "R: " << R.l2_norm() << std::endl;
    STOP_LOG("residual_and_jacobian()", "ComplexSolve");
}






bool
MAST::ComplexAssemblyBase::
sensitivity_assemble (const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {

    libmesh_error(); // not implemented. Call the blocked assembly instead.
    
    return false;
}
