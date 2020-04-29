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
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_transient_assembly.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "solver/transient_solver_base.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::LevelSetReinitializationTransientAssembly::
LevelSetReinitializationTransientAssembly():
MAST::TransientAssembly(),
_ref_sol   (nullptr) {
    
}




MAST::LevelSetReinitializationTransientAssembly::
~LevelSetReinitializationTransientAssembly() {
    
}



void
MAST::LevelSetReinitializationTransientAssembly::
set_reference_solution(const libMesh::NumericVector<Real>& sol) {

    libmesh_assert(!_ref_sol);
    _ref_sol = &sol;
}



void
MAST::LevelSetReinitializationTransientAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::TransientSolverBase
    &solver = dynamic_cast<MAST::TransientSolverBase&>(*_elem_ops);
    MAST::NonlinearSystem
    &transient_sys = _system->system();
    
    MAST::LevelSetTransientAssemblyElemOperations
    &level_set_elem_ops = dynamic_cast<MAST::LevelSetTransientAssemblyElemOperations&>
    (solver.get_elem_operation_object());
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &transient_sys);
    
    if (R) R->zero();
    if (J) J->zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, ref_sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = transient_sys.get_dof_map();
    
    
    
    // stores the localized solution, velocity, acceleration, etc. vectors.
    // These pointers will have to be deleted
    std::vector<libMesh::NumericVector<Real>*>
    local_qtys;

    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(transient_sys,
                                                    *_ref_sol).release());

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X, false);
    
    // ask the solver to localize the relevant solutions
    solver.build_local_quantities(X, local_qtys);
    
    libMesh::MeshBase::const_element_iterator       el     =
    transient_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    transient_sys.get_mesh().active_local_elements_end();
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        solver.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        solver.init(geom_elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        vec.setZero(ndofs);
        ref_sol.setZero(ndofs);
        mat.setZero(ndofs, ndofs);

        for (unsigned int i=0; i<dof_indices.size(); i++)
            ref_sol(i) = (*localized_solution)(dof_indices[i]);

        solver.set_element_data(dof_indices, local_qtys);
        level_set_elem_ops.set_elem_reference_solution(ref_sol);
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        solver.elem_calculations(J!=nullptr?true:false, vec, mat);
        solver.clear_elem();
        
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
        dof_indices.clear();
    }
    
    // delete pointers to the local solutions
    for (unsigned int i=0; i<local_qtys.size(); i++)
        delete local_qtys[i];
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J && close_matrix) J->close();
}
