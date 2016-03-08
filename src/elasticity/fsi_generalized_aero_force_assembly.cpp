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
_fluid_complex_solver(NULL),
_pressure_function(NULL),
_motion(NULL)
{ }






MAST::FSIGeneralizedAeroForceAssembly::~FSIGeneralizedAeroForceAssembly() {
    
}





void
MAST::FSIGeneralizedAeroForceAssembly::
init(MAST::FrequencyFunction& freq,
     MAST::ComplexSolverBase& complex_solver,
     MAST::SmallDisturbancePressureFunction& pressure_func,
     MAST::FlexibleSurfaceMotion& motion_func) {
    
    
    // make sure that this has not been already set
    libmesh_assert(!_freq);
    
    _freq                    = &freq;
    _fluid_complex_solver    = &complex_solver;
    _pressure_function       = &pressure_func;
    _motion                  = &motion_func;
}





void
MAST::FSIGeneralizedAeroForceAssembly::
assemble_generalized_aerodynamic_force_matrix(const libMesh::NumericVector<Real>& X,
                                              std::vector<libMesh::NumericVector<Real>*>& basis,
                                              ComplexMatrixX& mat) {
    
    libMesh::NonlinearImplicitSystem& nonlin_sys =
    dynamic_cast<libMesh::NonlinearImplicitSystem&>(_system->system());
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol;
    ComplexVectorX vec;
    RealMatrixX    basis_mat;

    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    std::auto_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(nonlin_sys,
                                                     X).release());
    
    // also create localized solution vectos for the bassis vectors
    const unsigned int
    n_basis = (unsigned int)basis.size();
    
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = _build_localized_vector(nonlin_sys, *basis[i]).release();
    

    mat.setZero(n_basis, n_basis);
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init_for_system_and_solution(*_system, X);
    
    
    
    // iterate over each structural mode to calculate the
    // fluid small-disturbance solution
    for (unsigned int i=0; i<n_basis; i++) {
        
        // set up the fluid flexible-surface boundary condition for this mode
        _motion->init(*_freq, *localized_basis[i]);
        
        
        // solve the complex smamll-disturbance fluid-equations
        _fluid_complex_solver->solve_block_matrix();
        
        
        // use this solution to initialize the structural boundary conditions
        _pressure_function->init(_fluid_complex_solver->get_assembly().base_sol(),
                                 _fluid_complex_solver->real_solution(),
                                 _fluid_complex_solver->imag_solution());
        
        
        // assemble the complex small-disturbance force vector force vector
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
            basis_mat.setZero(ndofs, n_basis);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
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
            nonlin_sys.get_dof_map().constrain_element_vector(v1, dof_indices);
            MAST::copy(v2, v1);
            vec.real() =  v2;
            
            // constrain and set the imag component
            MAST::copy(v1, vec.imag());
            nonlin_sys.get_dof_map().constrain_element_vector(v1, dof_indices);
            MAST::copy(v2, v1);
            vec.imag() =  v2;
            
            // project the force vector on all the structural modes for the
            // i^th column of the generalized aerodynamic force matrix
            mat.col(i) += basis_mat.transpose() * vec;
            
            physics_elem->detach_active_solution_function();
        }
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];
}


