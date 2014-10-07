
// MAST includes
#include "base/transient_assembly.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "solver/transient_solver_base.h"
#include "numerics/utility.h"


// libMesh includes
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::TransientAssembly::
TransientAssembly(MAST::PhysicsDisciplineBase& discipline,
                  MAST::TransientSolverBase& solver,
                  MAST::SystemInitialization& sys):
MAST::AssemblyBase(discipline, sys),
_transient_solver(solver) {
    libMesh::TransientNonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::TransientNonlinearImplicitSystem&>(_system.system());
    
    transient_sys.nonlinear_solver->residual_and_jacobian_object = this;
    transient_sys.attach_sensitivity_assemble_object(*this);
}



MAST::TransientAssembly::~TransientAssembly() {
    libMesh::TransientNonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::TransientNonlinearImplicitSystem&>(_system.system());
    
    transient_sys.nonlinear_solver->residual_and_jacobian_object = NULL;
    transient_sys.reset_sensitivity_assembly();
}



void
MAST::TransientAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libMesh::TransientNonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::TransientNonlinearImplicitSystem&>(_system.system());
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &transient_sys);
    
    if (R) R->zero();
    if (J) J->zero();
    
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, vel;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    
    // ask the solver to localize the relevant solutions
    _transient_solver._localize_solutions();
    
    const libMesh::NumericVector<Real>
    &solution = _transient_solver.solution(0),
    &velocity = _transient_solver.velocity(0);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vel.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol(i) = solution(dof_indices[i]);
            vel(i) = velocity(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(vel);
        
        // perform the element level calculations
        _transient_solver._elem_calculations(*physics_elem,
                                             dof_indices,
                                             J!=NULL?true:false,
                                             vec, mat);
        
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
            transient_sys.get_dof_map().constrain_element_matrix_and_vector(m, v, dof_indices);
        else if (R)
            transient_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        else
            transient_sys.get_dof_map().constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        if (R) R->add_vector(v, dof_indices);
        if (J) J->add_matrix(m, dof_indices);
    }
    
    
    if (R) R->close();
    if (J) J->close();
    
}



bool
MAST::TransientAssembly::
sensitivity_assemble (const libMesh::ParameterVector& parameters,
                      const unsigned int i,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libMesh::TransientNonlinearImplicitSystem& transient_sys =
    dynamic_cast<libMesh::TransientNonlinearImplicitSystem&>(_system.system());
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, vel;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    // ask the solver to localize the relevant solutions
    _transient_solver._localize_solutions();
    
    const libMesh::NumericVector<Real>
    &solution = _transient_solver.solution(0),
    &velocity = _transient_solver.velocity(0);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vel.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol(i) = solution(dof_indices[i]);
            vel(i) = velocity(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(vel);
        
        physics_elem->sensitivity_param = _discipline.get_parameter(parameters[i]);
        physics_elem->set_solution(sol);
        
        // perform the element level calculations
        _transient_solver._elem_sensitivity_calculations(*physics_elem, dof_indices, vec);
        
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        transient_sys.get_dof_map().constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
    }
    
    
    sensitivity_rhs.close();
    
    return true;
}


