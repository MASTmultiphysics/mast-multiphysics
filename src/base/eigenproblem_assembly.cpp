
// MAST includes
#include "base/eigenproblem_assembly.h"
#include "base/system_initialization.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"



MAST::EigenproblemAssembly::EigenproblemAssembly():
MAST::AssemblyBase(),
_if_exchange_A_and_B_matrices(false) {
    
}



MAST::EigenproblemAssembly::~EigenproblemAssembly() {
    
}




void
MAST::EigenproblemAssembly::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
    
    // now attach this to the system
    libMesh::EigenSystem& eigen_sys =
    dynamic_cast<libMesh::EigenSystem&>(system.system());
    eigen_sys.attach_assemble_object(*this);
    eigen_sys.attach_eigenproblem_sensitivity_assemble_object(*this);
}




void
MAST::EigenproblemAssembly::
clear_discipline_and_system( ) {
    
    if (_system && _discipline) {

        libMesh::EigenSystem& eigen_sys =
        dynamic_cast<libMesh::EigenSystem&>(_system->system());
        
        eigen_sys.reset_assembly();
        eigen_sys.reset_eigenproblem_sensitivity_assembly();
    }
    
    _discipline = NULL;
    _system     = NULL;
}




void
MAST::EigenproblemAssembly::assemble() {
    
    libMesh::EigenSystem& eigen_sys =
    dynamic_cast<libMesh::EigenSystem&>(_system->system());
    
    // zero the solution since it is not needed for eigenproblem
    eigen_sys.solution->zero();
    
    libMesh::SparseMatrix<Real>&  matrix_A =
    *(dynamic_cast<libMesh::EigenSystem&>(_system->system()).matrix_A);
    libMesh::SparseMatrix<Real>&  matrix_B =
    *(dynamic_cast<libMesh::EigenSystem&>(_system->system()).matrix_B);
    
    matrix_A.zero();
    matrix_B.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    libMesh::MeshBase::const_element_iterator       el     =
    eigen_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    eigen_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol.get()) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*_base_sol)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        
        _elem_calculations(*physics_elem, mat_A, mat_B);

        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);

        // constrain the element matrices.
        eigen_sys.get_dof_map().constrain_element_matrix(A, dof_indices);
        eigen_sys.get_dof_map().constrain_element_matrix(B, dof_indices);
        
        // add to the global matrices
        if (_if_exchange_A_and_B_matrices)
        {
            matrix_A.add_matrix (B, dof_indices); // load dependent
            matrix_B.add_matrix (A, dof_indices); // load independent
        }
        else
        {
            matrix_A.add_matrix (A, dof_indices); // load independent
            matrix_B.add_matrix (B, dof_indices); // load dependent
        }
    }
    
    
}





bool
MAST::EigenproblemAssembly::
sensitivity_assemble(const libMesh::ParameterVector& parameters,
                     const unsigned int i,
                     libMesh::SparseMatrix<Real>* sensitivity_A,
                     libMesh::SparseMatrix<Real>* sensitivity_B) {
    
    libMesh::EigenSystem& eigen_sys =
    dynamic_cast<libMesh::EigenSystem&>(_system->system());
    
    // zero the solution since it is not needed for eigenproblem
    eigen_sys.solution->zero();
    
    libMesh::SparseMatrix<Real>&  matrix_A =
    *(dynamic_cast<libMesh::EigenSystem*>(_system)->matrix_A);
    libMesh::SparseMatrix<Real>&  matrix_B =
    *(dynamic_cast<libMesh::EigenSystem*>(_system)->matrix_B);
    
    matrix_A.zero();
    matrix_B.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    RealMatrixX mat_A, mat_B;
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = eigen_sys.get_dof_map();
    std::auto_ptr<MAST::ElementBase> physics_elem;
    
    libMesh::MeshBase::const_element_iterator       el     =
    eigen_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    eigen_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, tell the element about it
        if (_base_sol.get()) {
            // make sure that the sensitivity is also provided
            libmesh_assert(_base_sol_sensitivity.get());
            
            // set the element's base solution
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*_base_sol)(dof_indices[i]);
        }
        
        physics_elem->set_base_solution(sol);
        
        // set the element's base solution sensitivity
        if (_base_sol_sensitivity.get()) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*_base_sol_sensitivity)(dof_indices[i]);
        }
        
        physics_elem->set_base_solution(sol, true);
        
        // tell the element about the sensitivity parameter
        physics_elem->sensitivity_param = _discipline->get_parameter(parameters[i]);
        
        _elem_sensitivity_calculations(*physics_elem, mat_A, mat_B);

        // copy to the libMesh matrix for further processing
        DenseRealMatrix A, B;
        MAST::copy(A, mat_A);
        MAST::copy(B, mat_B);
        
        // constrain the element matrices.
        eigen_sys.get_dof_map().constrain_element_matrix(A, dof_indices);
        eigen_sys.get_dof_map().constrain_element_matrix(B, dof_indices);
        
        // add to the global matrices
        if (_if_exchange_A_and_B_matrices)
        {
            matrix_A.add_matrix (B, dof_indices);
            matrix_B.add_matrix (A, dof_indices);
        }
        else
        {
            matrix_A.add_matrix (A, dof_indices);
            matrix_B.add_matrix (B, dof_indices);
        }
    }
    
    return true;
}

