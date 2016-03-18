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
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "solver/complex_solver_base.h"
#include "boundary_condition/flexible_surface_motion.h"
#include "fluid/small_disturbance_pressure_function.h"
#include "base/complex_assembly_base.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "property_cards/element_property_card_base.h"
#include "numerics/utility.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/parameter_vector.h"





MAST::FSIGeneralizedAeroForceAssembly::FSIGeneralizedAeroForceAssembly():
MAST::StructuralFluidInteractionAssembly(),
_freq(NULL),
_structural_comm(NULL),
_fluid_comm(NULL),
_fluid_complex_solver(NULL),
_pressure_function(NULL),
_motion(NULL)
{ }






MAST::FSIGeneralizedAeroForceAssembly::~FSIGeneralizedAeroForceAssembly() {
    
}





void
MAST::FSIGeneralizedAeroForceAssembly::
init(MAST::FrequencyFunction&                freq,
     libMesh::Parallel::Communicator&        structural_comm,
     libMesh::Parallel::Communicator&        fluid_comm,
     MAST::ComplexSolverBase*                complex_solver,
     MAST::SmallDisturbancePressureFunction* pressure_func,
     MAST::FlexibleSurfaceMotion*            motion_func) {
    
    
    // make sure that this has not been already set
    libmesh_assert(!_freq);
    
    _freq                    = &freq;
    _structural_comm         = &structural_comm;
    _fluid_comm              = &fluid_comm;
    
    
    // presently it is assumed that the structural comm is a subset of
    // the fluid comm.
    MPI_Comm
    s_comm      = _structural_comm->get(),
    f_comm      = _fluid_comm->get();
    
    if (s_comm != MPI_COMM_NULL) {
        
        MPI_Group
        s_group,
        f_group,
        union_group,
        diff_group;
        
        int
        diff_size = 0;
        
        // the difference between the union of s_comm and f_comm, and f_comm
        // should be null
        MPI_Comm_group      (s_comm, &s_group);
        MPI_Comm_group      (f_comm, &f_group);
        MPI_Group_union     (s_group, f_group, &union_group);
        MPI_Group_difference(union_group, f_group, &diff_group);
        MPI_Group_size      (diff_group, &diff_size);
        
        MPI_Group_free(&s_group);
        MPI_Group_free(&f_group);
        MPI_Group_free(&union_group);
        MPI_Group_free(&diff_group);
        
        libmesh_assert(!diff_size);

        // make sure the pointer is provided
        libmesh_assert(motion_func);
        _motion                  = motion_func;
    }
    else {
        
        // make sure the pointer is null
        libmesh_assert(!motion_func);
    }
    
    if (f_comm != MPI_COMM_NULL) {
        
        // make sure the pointers are provided
        libmesh_assert(complex_solver);
        libmesh_assert(pressure_func);
        
        _fluid_complex_solver    = complex_solver;
        _pressure_function       = pressure_func;
    }
    else {
        
        // make sure the pointers are null
        libmesh_assert(!complex_solver);
        libmesh_assert(!pressure_func);
    }
}




void
MAST::FSIGeneralizedAeroForceAssembly::clear_discipline_and_system() {
    
    _freq                  = NULL;
    _structural_comm       = NULL;
    _fluid_comm            = NULL;
    _motion                = NULL;
    _pressure_function     = NULL;
    _fluid_complex_solver  = NULL;
}



void
MAST::FSIGeneralizedAeroForceAssembly::
assemble_generalized_aerodynamic_force_matrix(std::vector<libMesh::NumericVector<Real>*>& basis,
                                              ComplexMatrixX& mat) {
    
    // make sure the data provided is sane
    libmesh_assert(_structural_comm);
    libmesh_assert(_fluid_comm);
    
    const bool
    if_structural =  _structural_comm->get() != MPI_COMM_NULL,
    if_fluid      =  _fluid_comm->get()      != MPI_COMM_NULL;
    
    // also create localized solution vectos for the bassis vectors
    unsigned int
    n_basis = (unsigned int)basis.size();
    

    
    if (if_structural) {
        
        libmesh_assert(_motion);
        libmesh_assert(_discipline);
        libmesh_assert(_system);
        libmesh_assert(n_basis);
        
        // make sure tha basis size is the same on this communicator
        _structural_comm->verify(n_basis);
    }
    else {
        
        // structural basis size should be zero on processors where the
        // structural comm is invalid
        libmesh_assert(!n_basis);
    }

    // send the basis size to the other communicator. We can do this
    // since it is assumed that the structural comm is a subset of the
    // fluid comm
    _fluid_comm->broadcast(n_basis);
    
    if (if_fluid) {
        
        libmesh_assert(_pressure_function);
        libmesh_assert(_fluid_complex_solver);
    }
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol;
    ComplexVectorX vec;
    RealMatrixX    basis_mat;

    mat.setZero(n_basis, n_basis);

    std::vector<libMesh::dof_id_type> dof_indices;
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);

    if (if_structural) {
        
        if (_base_sol)
            localized_solution.reset(_build_localized_vector(_system->system(),
                                                             *_base_sol).release());
        
        for (unsigned int i=0; i<n_basis; i++)
            localized_basis[i] = _build_localized_vector(_system->system(), *basis[i]).release();
        
        
        // if a solution function is attached, initialize it
        if (_sol_function && _base_sol)
            _sol_function->init_for_system_and_solution(*_system, *_base_sol);
    }
    
    
    // iterate over each structural mode to calculate the
    // fluid small-disturbance solution
    for (unsigned int i=0; i<n_basis; i++) {
        
        // set up the fluid flexible-surface boundary condition for this mode
        if (if_structural)
            _motion->init(*_freq, *localized_basis[i]);
        
        
        // solve the complex smamll-disturbance fluid-equations
        _fluid_complex_solver->solve_block_matrix();
        
        
        // use this solution to initialize the structural boundary conditions
        _pressure_function->init(_fluid_complex_solver->get_assembly().base_sol(),
                                 _fluid_complex_solver->real_solution(),
                                 _fluid_complex_solver->imag_solution());
        
        
        if (if_structural) {
            
            const libMesh::DofMap& dof_map = _system->system().get_dof_map();
            
            // assemble the complex small-disturbance force vector force vector
            libMesh::MeshBase::const_element_iterator       el     =
            _system->system().get_mesh().active_local_elements_begin();
            const libMesh::MeshBase::const_element_iterator end_el =
            _system->system().get_mesh().active_local_elements_end();
            
            
            for ( ; el != end_el; ++el) {
                
                const libMesh::Elem* elem = *el;
                
                dof_map.dof_indices (elem, dof_indices);
                
                physics_elem.reset(_build_elem(*elem).release());
                
                // get the solution
                unsigned int ndofs = (unsigned int)dof_indices.size();
                sol.setZero(ndofs);
                vec.setZero(ndofs);
                basis_mat.setZero(ndofs, n_basis);
                
                for (unsigned int i=0; i<dof_indices.size(); i++) {
                    
                    if (_base_sol)
                        sol(i) = (*localized_solution)(dof_indices[i]);
                    
                    for (unsigned int j=0; j<n_basis; j++)
                        basis_mat(i,j) = (*localized_basis[j])(dof_indices[i]);
                }
                
                
                physics_elem->set_solution(sol);
                sol.setZero();
                physics_elem->set_velocity(sol);     // set to zero value
                physics_elem->set_acceleration(sol); // set to zero value
                
                
                if (_sol_function)
                    physics_elem->attach_active_solution_function(*_sol_function);
                
                _elem_aerodynamic_force_calculations(*physics_elem, vec);
                
                DenseRealVector v1;
                RealVectorX     v2;
                
                // constrain and set the real component
                MAST::copy(v1, vec.real());
                dof_map.constrain_element_vector(v1, dof_indices);
                MAST::copy(v2, v1);
                vec.real() =  v2;
                
                // constrain and set the imag component
                MAST::copy(v1, vec.imag());
                dof_map.constrain_element_vector(v1, dof_indices);
                MAST::copy(v2, v1);
                vec.imag() =  v2;
                
                // project the force vector on all the structural modes for the
                // i^th column of the generalized aerodynamic force matrix
                mat.col(i) += basis_mat.transpose() * vec;
                
                physics_elem->detach_active_solution_function();
            }
        }
    }
    
    
    
    if (if_structural) {
        
        // if a solution function is attached, clear it
        if (_sol_function)
            _sol_function->clear();
        
        
        // delete the localized basis vectors
        for (unsigned int i=0; i<basis.size(); i++)
            delete localized_basis[i];
    }
    
    // sum the matrix and provide it to each processor
    // this assumes that the structural comm is a subset of fluid comm
    MAST::parallel_sum(*_fluid_comm, mat);
}




void
MAST::FSIGeneralizedAeroForceAssembly::
assemble_reduced_order_quantity
(std::vector<libMesh::NumericVector<Real>*>& basis,
 std::map<MAST::StructuralQuantityType, RealMatrixX*>& mat_qty_map) {

    // make sure the data provided is sane
    libmesh_assert(_structural_comm);
    libmesh_assert(_fluid_comm);
    
    const bool
    if_structural =  _structural_comm->get() != MPI_COMM_NULL;

    unsigned int
    n_basis = (unsigned int)basis.size();

    
    // the structural comm should have the same basis
    if (if_structural)
        _structural_comm->verify(n_basis);
    else {
        
        // basis size should be zero on ranks where structural comm is not valid
        libmesh_assert(!n_basis);
    }

    // tell all processors about the basis size
    _fluid_comm->broadcast(n_basis);
    
    // initialize the quantities to zero matrices
    std::map<MAST::StructuralQuantityType, RealMatrixX*>::iterator
    it  = mat_qty_map.begin(),
    end = mat_qty_map.end();
    
    for ( ; it != end; it++)
        *it->second = RealMatrixX::Zero(n_basis, n_basis);

    
    
    if (if_structural) {
        
        libmesh_assert(_motion);
        libmesh_assert(_discipline);
        libmesh_assert(_system);
        
        
        // make sure that a valid basis is provided on the structural comm
        libmesh_assert(n_basis);

        // iterate over each element, initialize it and get the relevant
        // analysis quantities
        RealVectorX vec, sol;
        RealMatrixX mat, basis_mat;
        
        std::vector<libMesh::dof_id_type> dof_indices;
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        std::auto_ptr<MAST::ElementBase> physics_elem;
        
        std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
        if (_base_sol)
            localized_solution.reset(_build_localized_vector(_system->system(),
                                                             *_base_sol).release());
        
        // also create localized solution vectos for the bassis vectors
        
        
        std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
        for (unsigned int i=0; i<n_basis; i++)
            localized_basis[i] = _build_localized_vector(_system->system(),
                                                         *basis[i]).release();
        
        
        // if a solution function is attached, initialize it
        if (_sol_function && _base_sol)
            _sol_function->init_for_system_and_solution(*_system, *_base_sol);
        
        
        libMesh::MeshBase::const_element_iterator       el     =
        _system->system().get_mesh().active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el =
        _system->system().get_mesh().active_local_elements_end();
        
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
    }
    
    // sum the matrix and provide it to each processor
    it  = mat_qty_map.begin();
    end = mat_qty_map.end();
    
    for ( ; it != end; it++)
        // this assumes that the structural comm is a subset of fluid comm
        MAST::parallel_sum(*_fluid_comm, *(it->second));
}

