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
#include "elasticity/fsi_generalized_aero_force_assembly.h"
#include "elasticity/structural_element_base.h"
#include "elasticity/structural_assembly.h"
#include "elasticity/fluid_structure_assembly_elem_operations.h"
#include "base/complex_mesh_field_function.h"
#include "base/complex_assembly_base.h"
#include "base/complex_assembly_elem_operations.h"
#include "base/physics_discipline_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "fluid/pressure_function.h"
#include "fluid/frequency_domain_pressure_function.h"
#include "property_cards/element_property_card_base.h"
#include "solver/complex_solver_base.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"





MAST::FSIGeneralizedAeroForceAssembly::FSIGeneralizedAeroForceAssembly():
MAST::StructuralFluidInteractionAssembly(),
_fluid_complex_solver           (nullptr),
_fluid_complex_assembly         (nullptr),
_pressure_function              (nullptr),
_freq_domain_pressure_function  (nullptr),
_complex_displ                  (nullptr)
{ }






MAST::FSIGeneralizedAeroForceAssembly::~FSIGeneralizedAeroForceAssembly() {
    
}





void
MAST::FSIGeneralizedAeroForceAssembly::
init(MAST::FluidStructureAssemblyElemOperations&  fsi_elem_ops,
     MAST::ComplexSolverBase&                     complex_solver,
     MAST::ComplexAssemblyBase&                   complex_assembly,
     MAST::ComplexAssemblyElemOperations&         fluid_elem_ops,
     MAST::PressureFunction&                      pressure_func,
     MAST::FrequencyDomainPressureFunction&       freq_pressure_func,
     MAST::ComplexMeshFieldFunction&              displ_func) {
    
    libmesh_assert(!_fluid_complex_solver);
    
    _complex_displ                  = &displ_func;
    _fluid_complex_solver           = &complex_solver;
    _fluid_complex_assembly         = &complex_assembly;
    _pressure_function              = &pressure_func;
    _freq_domain_pressure_function  = &freq_pressure_func;

    this->set_elem_operation_object(fsi_elem_ops);
    complex_solver.set_assembly(complex_assembly);
    complex_assembly.set_elem_operation_object(fluid_elem_ops);
}




void
MAST::FSIGeneralizedAeroForceAssembly::clear_discipline_and_system() {

    _fluid_complex_assembly->clear_elem_operation_object();
    _fluid_complex_solver->clear_assembly();
    
    _elem_ops                          = nullptr;
    _complex_displ                     = nullptr;
    _pressure_function                 = nullptr;
    _freq_domain_pressure_function     = nullptr;
    _fluid_complex_solver              = nullptr;
    _fluid_complex_assembly            = nullptr;
    
    this->clear_elem_operation_object();
    MAST::StructuralFluidInteractionAssembly::clear_discipline_and_system();
}



void
MAST::FSIGeneralizedAeroForceAssembly::
assemble_generalized_aerodynamic_force_matrix
(std::vector<libMesh::NumericVector<Real>*>& basis,
 ComplexMatrixX& mat,
 MAST::Parameter* p) {
    
    // make sure the data provided is sane
    libmesh_assert(_complex_displ);
    
    
    // also create localized solution vectos for the bassis vectors
    unsigned int
    n_basis = (unsigned int)basis.size();
    

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX    sol;
    ComplexVectorX vec;
    RealMatrixX    basis_mat;

    mat.setZero(n_basis, n_basis);

    std::vector<libMesh::dof_id_type> dof_indices;
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_zero;
    std::vector<libMesh::NumericVector<Real>*> localized_basis(n_basis);

    if (_base_sol)
        localized_solution.reset(build_localized_vector(_system->system(),
                                                         *_base_sol).release());
    
    for (unsigned int i=0; i<n_basis; i++)
        localized_basis[i] = build_localized_vector(_system->system(), *basis[i]).release();
    
    //create a zero-clone copy for the imaginary component of the solution
    localized_zero.reset(localized_basis[0]->zero_clone().release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function && _base_sol)
        _sol_function->init( *_base_sol, false);
    
    MAST::FluidStructureAssemblyElemOperations&
    ops = dynamic_cast<MAST::FluidStructureAssemblyElemOperations&>(*_elem_ops);
    
    // iterate over each structural mode to calculate the
    // fluid small-disturbance solution
    for (unsigned int i=0; i<n_basis; i++) {
        
        // set up the fluid flexible-surface boundary condition for this mode
        _complex_displ->clear();
        _complex_displ->init(*localized_basis[i], *localized_zero);
        
        
        // solve the complex smamll-disturbance fluid-equations
        _fluid_complex_solver->solve_block_matrix(p);
        
        // use this solution to initialize the structural boundary conditions
        _pressure_function->init(_fluid_complex_assembly->base_sol());
        
        // use this solution to initialize the structural boundary conditions
        _freq_domain_pressure_function->init
        (_fluid_complex_assembly->base_sol(),
         _fluid_complex_solver->real_solution(p != nullptr),
         _fluid_complex_solver->imag_solution(p != nullptr));

        
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        
        // assemble the complex small-disturbance force vector force vector
        libMesh::MeshBase::const_element_iterator       el     =
        _system->system().get_mesh().active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el =
        _system->system().get_mesh().active_local_elements_end();
        
        
        
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
            vec.setZero(ndofs);
            basis_mat.setZero(ndofs, n_basis);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                
                if (_base_sol)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                
                for (unsigned int j=0; j<n_basis; j++)
                    basis_mat(i,j) = (*localized_basis[j])(dof_indices[i]);
            }
            
            
            ops.set_elem_solution(sol);
            sol.setZero();
            ops.set_elem_velocity(sol);     // set to zero value
            ops.set_elem_acceleration(sol); // set to zero value
            
            
//            if (_sol_function)
//                physics_elem->attach_active_solution_function(*_sol_function);
            
            ops.elem_aerodynamic_force_calculations(vec);
            ops.clear_elem();
            
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
            
//            physics_elem->detach_active_solution_function();
        }
    }
    
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    
    // delete the localized basis vectors
    for (unsigned int i=0; i<basis.size(); i++)
        delete localized_basis[i];
    
    // sum the matrix and provide it to each processor
    // this assumes that the structural comm is a subset of fluid comm
    MAST::parallel_sum(_system->system().comm(), mat);
}

