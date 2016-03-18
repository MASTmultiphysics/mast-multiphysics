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
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "numerics/utility.h"
#include "base/real_output_function.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/parameter_vector.h"



MAST::StructuralFluidInteractionAssembly::
StructuralFluidInteractionAssembly():
MAST::NonlinearImplicitAssembly(),
_base_sol(NULL),
_base_sol_sensitivity(NULL) {
    
}



MAST::StructuralFluidInteractionAssembly::
~StructuralFluidInteractionAssembly() {
    
}


void
MAST::StructuralFluidInteractionAssembly::
clear_discipline_and_system( ) {

    _base_sol             = NULL;
    _base_sol_sensitivity = NULL;
    
    MAST::NonlinearImplicitAssembly::clear_discipline_and_system();
}




void
MAST::StructuralFluidInteractionAssembly::set_base_solution(const libMesh::NumericVector<Real>& sol,
                                                            bool if_sens) {
    
    if (!if_sens) {
        
        // make sure that the pointer has been cleared
        libmesh_assert(!_base_sol);
        _base_sol             = &sol;
    }
    else {

        libmesh_assert(!_base_sol_sensitivity);
        _base_sol_sensitivity = &sol;
    }
}




void
MAST::StructuralFluidInteractionAssembly::clear_base_solution(bool if_sens) {
    
    if (!if_sens)
        _base_sol             = NULL;
    else
        _base_sol_sensitivity = NULL;
}




void
MAST::StructuralFluidInteractionAssembly::
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
        _elem_calculations(*physics_elem, true, vec, mat);
        
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
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J) J->close();
}



void
MAST::StructuralFluidInteractionAssembly::
assemble_reduced_order_quantity
(std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat, basis_mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    if (_base_sol)
        localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                         *_base_sol).release());
    
    // also create localized solution vectos for the bassis vectors
    const unsigned int
    n_basis = (unsigned int)basis.size();
    
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = _build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init_for_system_and_solution(*_system, *_base_sol);
    
    
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
        basis_mat.setZero(ndofs, n_basis);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            
            if (_base_sol)
                sol(i) = (*localized_solution)(dof_indices[i]);
            
            
            for (unsigned int j=0; j<n_basis; j++)
                basis_mat(i,j) = (*localized_basis[j])(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(vec);     // set to zero value
        physics_elem->set_acceleration(vec); // set to zero value
        
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // now iterative over all qty types in the map and assemble them
        std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
        it   = mat_qty_map.begin(),
        end  = mat_qty_map.end();
        
        for ( ; it != end; it++) {
            
            _qty_type = it->first;
            _elem_calculations(*physics_elem, true, vec, mat);
            
            DenseRealMatrix m;
            MAST::copy(m, mat);
            dof_map.constrain_element_matrix(m, dof_indices);
            MAST::copy(mat, m);
            
            // now add to the reduced order matrix
            (*it->second) += basis_mat.transpose() * mat * basis_mat;
        }
        
        physics_elem->detach_active_solution_function();
        
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];

    // sum the matrix and provide it to each processor
    std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
    it  = mat_qty_map.begin(),
    end = mat_qty_map.end();
    
    
    for ( ; it != end; it++)
        MAST::parallel_sum(_system->system().comm(), *(it->second));
}





void
MAST::StructuralFluidInteractionAssembly::
assemble_reduced_order_quantity_sensitivity
(const libMesh::ParameterVector& parameters,
 const unsigned int i,
 std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {
    
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, dsol;
    RealMatrixX mat, basis_mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    
    if (_base_sol) {
        
        // make sure that the solution sensitivity is provided
        libmesh_assert(_base_sol_sensitivity);
        
        localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                         *_base_sol).release());
        localized_solution_sens.reset(_build_localized_vector(nonlin_sys,
                                                              *_base_sol_sensitivity).release());
    }
    
    // also create localized solution vectos for the bassis vectors
    const unsigned int
    n_basis = (unsigned int)basis.size();
    
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = _build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init_for_system_and_solution(*_system, *_base_sol);
    
    
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
        basis_mat.setZero(ndofs, n_basis);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol(i)  = (*localized_solution)(dof_indices[i]);
            dsol(i) = (*localized_solution_sens)(dof_indices[i]);
            for (unsigned int j=0; j<n_basis; j++)
                basis_mat(i,j) = (*localized_basis[j])(dof_indices[i]);
        }
        
        physics_elem->sensitivity_param  = _discipline->get_parameter(&(parameters[i].get()));
        physics_elem->set_solution(sol);
        physics_elem->set_solution(dsol, true);
        physics_elem->set_velocity(vec);     // set to zero value
        physics_elem->set_acceleration(vec); // set to zero value
        
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // now iterative over all qty types in the map and assemble them
        std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
        it   = mat_qty_map.begin(),
        end  = mat_qty_map.end();
        
        for ( ; it != end; it++) {
            
            _qty_type = it->first;
            _elem_sensitivity_calculations(*physics_elem, true, vec, mat);

            DenseRealMatrix m;
            MAST::copy(m, mat);
            dof_map.constrain_element_matrix(m, dof_indices);
            MAST::copy(mat, m);
            
            // now add to the reduced order matrix
            (*it->second) += basis_mat.transpose() * mat * basis_mat;
        }
        
        physics_elem->detach_active_solution_function();
        
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];
    
    // sum the matrix and provide it to each processor
    std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
    it  = mat_qty_map.begin(),
    end = mat_qty_map.end();
    
    
    for ( ; it != end; it++)
        MAST::parallel_sum(_system->system().comm(), *(it->second));
}





std::auto_ptr<MAST::ElementBase>
MAST::StructuralFluidInteractionAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}




void
MAST::StructuralFluidInteractionAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {

    libmesh_assert(if_jac);
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealMatrixX
    dummy = RealMatrixX::Zero(mat.rows(), mat.cols());
    mat.setZero();
    vec.setZero();
    
    switch (_qty_type) {
        case MAST::MASS: {
            
            e.inertial_residual(true, vec, mat, dummy, dummy);
        }
            break;

            
        case MAST::DAMPING: {
            
            e.inertial_residual(true, vec, dummy, mat, dummy);
            e.side_external_residual(true,
                                     vec,
                                     mat,
                                     dummy,
                                     _discipline->side_loads());
            e.volume_external_residual(true,
                                       vec,
                                       mat,
                                       dummy,
                                       _discipline->volume_loads());
        }
            break;

        case MAST::STIFFNESS: {

            // create
            
            e.internal_residual(true, vec, mat, false);
            e.inertial_residual(true, vec, dummy, dummy, mat);
            e.side_external_residual(true,
                                     vec,
                                     dummy,
                                     mat,
                                     _discipline->side_loads());
            e.volume_external_residual(true,
                                       vec,
                                       dummy,
                                       mat,
                                       _discipline->volume_loads());
        }
            break;
            
        default:
            libmesh_error(); // should not get here
    }

    
}






void
MAST::StructuralFluidInteractionAssembly::
_elem_aerodynamic_force_calculations(MAST::ElementBase& elem,
                                     ComplexVectorX& vec) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    ComplexMatrixX
    dummy = ComplexMatrixX::Zero(vec.size(), vec.size());
    vec.setZero();
    
    e.side_frequency_domain_external_residual(false,
                                              vec,
                                              dummy,
                                              _discipline->side_loads());
    e.volume_frequency_domain_external_residual(false,
                                                vec,
                                                dummy,
                                                _discipline->volume_loads());
}




void
MAST::StructuralFluidInteractionAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               bool if_jac,
                               RealVectorX& vec,
                               RealMatrixX& mat) {
    
    libmesh_assert(if_jac);

    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealMatrixX
    dummy = RealMatrixX::Zero(vec.size(), vec.size());
    mat.setZero();
    vec.setZero();
    
    switch (_qty_type) {
        case MAST::MASS: {

            e.inertial_residual_sensitivity(true, vec, mat, dummy, dummy);
        }
            break;
            
        case MAST::DAMPING: {
            
            e.inertial_residual_sensitivity(true, vec, dummy, mat, dummy);
            e.side_external_residual_sensitivity(true,
                                                 vec,
                                                 mat,
                                                 dummy,
                                                 _discipline->side_loads());
            e.volume_external_residual_sensitivity(true,
                                                   vec,
                                                   mat,
                                                   dummy,
                                                   _discipline->volume_loads());
        }
            break;
            
            
        case MAST::STIFFNESS: {
            
            e.internal_residual_sensitivity(true, vec, mat, false);
            e.inertial_residual_sensitivity(true, vec, dummy, dummy, mat);
            e.side_external_residual_sensitivity(true,
                                                 vec,
                                                 dummy,
                                                 mat,
                                                 _discipline->side_loads());
            e.volume_external_residual_sensitivity(true,
                                                   vec,
                                                   dummy,
                                                   mat,
                                                   _discipline->volume_loads());
        }
            break;
        
        default:
            libmesh_error(); // should not get here
    }

    
}



