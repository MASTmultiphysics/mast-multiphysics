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
_base_sol(nullptr),
_base_sol_sensitivity(nullptr) {
    
}



MAST::StructuralFluidInteractionAssembly::
~StructuralFluidInteractionAssembly() {
    
}




void
MAST::StructuralFluidInteractionAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;

    // this method does not provide residual assembly routines. Hence,
    // nothing is done to attach itself to the system.
}



void
MAST::StructuralFluidInteractionAssembly::
clear_discipline_and_system( ) {

    _base_sol             = nullptr;
    _base_sol_sensitivity = nullptr;
    
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
        _base_sol             = nullptr;
    else
        _base_sol_sensitivity = nullptr;
}




void
MAST::StructuralFluidInteractionAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    // this method should not be called for this class.
    libmesh_error();
}



void
MAST::StructuralFluidInteractionAssembly::
assemble_reduced_order_quantity
(std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    unsigned int
    n_basis = (unsigned int)basis.size();

    // initialize the quantities to zero matrices
    std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
    it  = mat_qty_map.begin(),
    end = mat_qty_map.end();
    
    for ( ; it != end; it++)
        *it->second = RealMatrixX::Zero(n_basis, n_basis);

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
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = _build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init( *_base_sol);
    
    
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
    
    unsigned int
    n_basis = (unsigned int)basis.size();
    
    
    // initialize the quantities to zero matrices
    std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
    it  = mat_qty_map.begin(),
    end = mat_qty_map.end();
    
    for ( ; it != end; it++)
        *it->second = RealMatrixX::Zero(n_basis, n_basis);

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
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = _build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init( *_base_sol);
    
    
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
            
            if (_base_sol) {
                
                sol(i)  = (*localized_solution)(dof_indices[i]);
                dsol(i) = (*localized_solution_sens)(dof_indices[i]);
            }
            
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
            
            e.internal_residual(true, vec, mat);
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
    
    e.linearized_frequency_domain_side_external_residual(false,
                                                         vec,
                                                         dummy,
                                                         _discipline->side_loads());
    e.linearized_frequency_domain_volume_external_residual(false,
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
            
            e.internal_residual_sensitivity(true, vec, mat);
            
            // if the linearization is about a base state, then the sensitivity of
            // the base state will influence the sensitivity of the Jacobian
            if (_base_sol)
                e.internal_residual_jac_dot_state_sensitivity(mat);

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



