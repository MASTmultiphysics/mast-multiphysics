/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include "elasticity/structural_assembly.h"
#include "elasticity/fluid_structure_assembly_elem_operations.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/assembly_elem_operation.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"



MAST::StructuralFluidInteractionAssembly::
StructuralFluidInteractionAssembly():
MAST::AssemblyBase(),
_base_sol             (nullptr),
_base_sol_sensitivity (nullptr) {
    
}



MAST::StructuralFluidInteractionAssembly::
~StructuralFluidInteractionAssembly() {
    
}





void
MAST::StructuralFluidInteractionAssembly::
clear_discipline_and_system( ) {
    
    _base_sol             = nullptr;
    _base_sol_sensitivity = nullptr;
    
    MAST::AssemblyBase::clear_discipline_and_system();
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
assemble_reduced_order_quantity
(std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {
    
    libmesh_assert(_elem_ops);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    if (_base_sol)
        localized_solution.reset(build_localized_vector(nonlin_sys,
                                                        *_base_sol).release());
    
    // also create localized solution vectos for the bassis vectors
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init( *_base_sol, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::FluidStructureAssemblyElemOperations
    &ops = dynamic_cast<MAST::FluidStructureAssemblyElemOperations&>(*_elem_ops);
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
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
        
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);
        
        MAST::GeomElem geom_elem;
        _elem_ops->set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        _elem_ops->init(geom_elem);
        _elem_ops->set_elem_solution(sol);
        _elem_ops->set_elem_velocity(vec);     // set to zero value
        _elem_ops->set_elem_acceleration(vec); // set to zero value
        
        
        // now iterative over all qty types in the map and assemble them
        it   = mat_qty_map.begin();
        end  = mat_qty_map.end();
        
        for ( ; it != end; it++) {
            
            ops.set_qty_to_evaluate(it->first);
            ops.elem_calculations(true, vec, mat);
            
            DenseRealMatrix m;
            MAST::copy(m, mat);
            dof_map.constrain_element_matrix(m, dof_indices);
            MAST::copy(mat, m);
            
            // now add to the reduced order matrix
            (*it->second) += basis_mat.transpose() * mat * basis_mat;
        }
        
        _elem_ops->clear_elem();
        //        physics_elem->detach_active_solution_function();
        
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];
    
    // sum the matrix and provide it to each processor
    it  = mat_qty_map.begin();
    end = mat_qty_map.end();
    
    
    for ( ; it != end; it++)
        MAST::parallel_sum(_system->system().comm(), *(it->second));
}





void
MAST::StructuralFluidInteractionAssembly::
assemble_reduced_order_quantity_sensitivity
(const MAST::FunctionBase& f,
 std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {
    
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    
    if (_base_sol) {
        
        // make sure that the solution sensitivity is provided
        libmesh_assert(_base_sol_sensitivity);
        
        localized_solution.reset(build_localized_vector(nonlin_sys,
                                                        *_base_sol).release());
        localized_solution_sens.reset(build_localized_vector(nonlin_sys,
                                                             *_base_sol_sensitivity).release());
    }
    
    // also create localized solution vectos for the bassis vectors
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = build_localized_vector(nonlin_sys, *basis[i]).release();
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init( *_base_sol, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::FluidStructureAssemblyElemOperations
    &ops = dynamic_cast<MAST::FluidStructureAssemblyElemOperations&>(*_elem_ops);
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dsol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        basis_mat.setZero(ndofs, n_basis);
        
        MAST::GeomElem geom_elem;
        _elem_ops->set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        _elem_ops->init(geom_elem);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            
            if (_base_sol) {
                
                ops.use_base_sol_for_sensitivity(true);
                sol(i)  = (*localized_solution)(dof_indices[i]);
                dsol(i) = (*localized_solution_sens)(dof_indices[i]);
            }
            
            for (unsigned int j=0; j<n_basis; j++)
                basis_mat(i,j) = (*localized_basis[j])(dof_indices[i]);
        }
        
        _elem_ops->set_elem_solution(sol);
        _elem_ops->set_elem_solution_sensitivity(dsol);
        _elem_ops->set_elem_velocity(vec);     // set to zero value
        _elem_ops->set_elem_acceleration(vec); // set to zero value
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);
        
        // now iterative over all qty types in the map and assemble them
        it   = mat_qty_map.begin();
        end  = mat_qty_map.end();
        
        for ( ; it != end; it++) {
            
            ops.set_qty_to_evaluate(it->first);
            ops.elem_sensitivity_calculations(f, true, vec, mat);
            
            DenseRealMatrix m;
            MAST::copy(m, mat);
            dof_map.constrain_element_matrix(m, dof_indices);
            MAST::copy(mat, m);
            
            // now add to the reduced order matrix
            (*it->second) += basis_mat.transpose() * mat * basis_mat;
        }
        
        _elem_ops->clear_elem();
        //        physics_elem->detach_active_solution_function();
        
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];
    
    // sum the matrix and provide it to each processor
    it  = mat_qty_map.begin();
    end = mat_qty_map.end();
    
    
    for ( ; it != end; it++)
        MAST::parallel_sum(_system->system().comm(), *(it->second));
}

