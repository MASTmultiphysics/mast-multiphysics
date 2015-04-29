/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/physics_discipline_base.h"
#include "numerics/utility.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"


MAST::StructuralModalEigenproblemAssembly::
StructuralModalEigenproblemAssembly():
MAST::EigenproblemAssembly() { }




MAST::StructuralModalEigenproblemAssembly::
~StructuralModalEigenproblemAssembly()
{ }



void
MAST::StructuralModalEigenproblemAssembly::assemble() {
    
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
    RealVectorX sol, dummy;
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
        dummy.setZero(ndofs);
        mat_A.setZero(ndofs, ndofs);
        mat_B.setZero(ndofs, ndofs);
        
        // if the base solution is provided, then tell the element about it
        if (_base_sol.get()) {
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*_base_sol)(dof_indices[i]);
        }
        
        physics_elem->set_solution(sol);
        physics_elem->set_velocity(dummy);
        physics_elem->set_acceleration(dummy);
        
        
        // set the incompatible mode solution if required by the
        // element
        MAST::StructuralElementBase& p_elem =
        dynamic_cast<MAST::StructuralElementBase&>(*physics_elem);
        if (p_elem.if_incompatible_modes()) {
            // check if the vector exists in the map
            if (!_incompatible_sol.count(elem))
                _incompatible_sol[elem] = RealVectorX::Zero(p_elem.incompatible_mode_size());
            p_elem.set_incompatible_mode_solution(_incompatible_sol[elem]);
        }

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



std::auto_ptr<MAST::ElementBase>
MAST::StructuralModalEigenproblemAssembly::_build_elem(const libMesh::Elem& elem) {
    
    
    const MAST::ElementPropertyCardBase& p =
    dynamic_cast<const MAST::ElementPropertyCardBase&>(_discipline->get_property_card(elem));
    
    MAST::ElementBase* rval =
    MAST::build_structural_element(*_system, elem, p, false).release();
    
    return std::auto_ptr<MAST::ElementBase>(rval);
}



void
MAST::StructuralModalEigenproblemAssembly::
_elem_calculations(MAST::ElementBase& elem,
                   RealMatrixX& mat_A,
                   RealMatrixX& mat_B) {
    
    MAST::StructuralElementBase& e =
    dynamic_cast<MAST::StructuralElementBase&>(elem);
    
    RealVectorX vec = RealVectorX::Zero(mat_A.rows()); // dummy vector
    RealMatrixX mat = RealMatrixX::Zero(mat_A.rows(), mat_A.cols()); // dummy matrix
    mat_A.setZero();
    mat_B.setZero();
    
    // calculate the Jacobian components
    e.internal_residual(true, vec, mat_A, true);
//    e.side_external_residual<Real>(true, vec, mat_A, _discipline->side_loads());
//    e.volume_external_residual<Real>(true, vec, mat_A, _discipline->volume_loads());
    
    // calculate the mass matrix components
    e.inertial_residual(true, vec, mat_B, mat, mat_A);
}
        



void
MAST::StructuralModalEigenproblemAssembly::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               RealMatrixX& mat_A,
                               RealMatrixX& mat_B) {
    // to be implelmented
    libmesh_error();
}


