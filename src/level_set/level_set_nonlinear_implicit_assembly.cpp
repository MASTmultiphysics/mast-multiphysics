/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2018  Manav Bhatia
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
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_intersection.h"
#include "level_set/sub_cell_fe.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "base/elem_base.h"
#include "numerics/utility.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::LevelSetNonlinearImplicitAssembly::
LevelSetNonlinearImplicitAssembly():
MAST::NonlinearImplicitAssembly(),
_level_set     (nullptr),
_intersection  (nullptr) {
    
}



MAST::LevelSetNonlinearImplicitAssembly::~LevelSetNonlinearImplicitAssembly() {
 
    if (_intersection)
        delete _intersection;
}




void
MAST::LevelSetNonlinearImplicitAssembly::
attach_discipline_and_system(MAST::NonlinearImplicitAssemblyElemOperations& elem_ops,
                             MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system,
                             MAST::FieldFunction<Real>& level_set) {
    
    MAST::NonlinearImplicitAssembly::attach_discipline_and_system(elem_ops,
                                                                  discipline,
                                                                  system);
    _system->system().attach_constraint_object(*this);
    _level_set = &level_set;
    _intersection = new MAST::LevelSetIntersection;
}



void
MAST::LevelSetNonlinearImplicitAssembly::clear_discipline_and_system() {
    
    MAST::NonlinearImplicitAssembly::clear_discipline_and_system();
    _level_set = nullptr;
    
    if (_intersection) {
        delete _intersection;
        _intersection = nullptr;
    }
}



void
MAST::LevelSetNonlinearImplicitAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                    X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);
        
        
        if (!_intersection->if_intersection() &&
            !_intersection->if_elem_on_positive_phi()) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            // copy to the libMesh matrix for further processing
            DenseRealVector v;
            DenseRealMatrix m;
            if (R)
                v.resize(ndofs);
            if (J)
                m.resize(ndofs, ndofs);
            
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
        else {
            
            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                dof_map.dof_indices (elem, dof_indices);
                
                _implicit_elem_ops->init(*sub_elem);
                
                // get the solution
                unsigned int ndofs = (unsigned int)dof_indices.size();
                sol.setZero(ndofs);
                vec.setZero(ndofs);
                mat.setZero(ndofs, ndofs);
                
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                
                _implicit_elem_ops->set_elem_solution(sol);
                
//                if (_sol_function)
//                    physics_elem->attach_active_solution_function(*_sol_function);
                
                // perform the element level calculations
                _implicit_elem_ops->elem_calculations(J!=nullptr?true:false,
                                                      vec, mat);
                
//                physics_elem->detach_active_solution_function();
                
                _implicit_elem_ops->clear_elem();
                
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
        }
        
        _intersection->clear();
    }
    
    // call the post assembly object, if provided by user
    if (_post_assembly)
        _post_assembly->post_assembly(X, R, J, S);
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J) J->close();
}



void
MAST::LevelSetNonlinearImplicitAssembly::constrain() {

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    

    std::vector<libMesh::dof_id_type>
    dof_indices,
    negative_and_intersected_dof_indices,
    intersected_dof_indices;
    
    negative_and_intersected_dof_indices.reserve(dof_map.n_local_dofs());
    intersected_dof_indices.reserve(dof_map.n_local_dofs());
    
    // our intent is to constrain only those dofs that belong to the
    // unintersected elements on the negative phi side of level set, AND
    // if they do not belong to any elements intersected by the level set.
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        MAST::LevelSetIntersection intersect;
        intersect.init(*_level_set, *elem, nonlin_sys.time);
        
        // if the element is entirely on the negative side of the level set,
        // we will constrain all dofs of the element to zero
        
        // this identifies if the element is NOT entirely on the positive side
        if (!intersect.if_elem_on_positive_phi()) {
            
            dof_indices.clear();
            dof_map.dof_indices(elem, dof_indices);
            negative_and_intersected_dof_indices.insert(negative_and_intersected_dof_indices.end(),
                                                        dof_indices.begin(),
                                                        dof_indices.end());
        }
        
        // this identifies if there is an intersection on the element
        if (intersect.if_intersection()) {
            
            dof_indices.clear();
            dof_map.dof_indices(elem, dof_indices);
            intersected_dof_indices.insert(intersected_dof_indices.end(),
                                           dof_indices.begin(),
                                           dof_indices.end());
        }
    }
    
    // now that we have the dofs of the negative and intersected regions,
    // we use only those dofs that have no connection to the unintersected
    // negative phi dofs.
    std::set<libMesh::dof_id_type>
    constrained_dofs(negative_and_intersected_dof_indices.begin(),
                     negative_and_intersected_dof_indices.end()),
    intersected_dofs(intersected_dof_indices.begin(),
                     intersected_dof_indices.end());
    
    // clear the unwanted space
    negative_and_intersected_dof_indices.clear();
    intersected_dof_indices.clear();
    
    dof_indices.clear();
    dof_indices.reserve(constrained_dofs.size());
    
    std::set_difference(constrained_dofs.begin(),
                        constrained_dofs.end(),
                        intersected_dofs.begin(),
                        intersected_dofs.end(),
                        std::inserter(dof_indices, dof_indices.begin()));

    constrained_dofs.clear();
    intersected_dofs.clear();
    
    // now that we have the final set of dofs to be constrained, clear
    // temprary storage and operate on these
    std::vector<libMesh::dof_id_type>::const_iterator
    dof_it  = dof_indices.begin(),
    dof_end = dof_indices.end();
    
    for ( ; dof_it != dof_end; dof_it++) {
        
        libMesh::DofConstraintRow c_row;
        dof_map.add_constraint_row(*dof_it, c_row, true);
    }
}



std::unique_ptr<MAST::FEBase>
MAST::LevelSetNonlinearImplicitAssembly::build_fe(const libMesh::Elem& elem) {
    
    std::unique_ptr<MAST::FEBase> fe;
    
    if (_elem_ops->if_use_local_elem() &&
        elem.dim() < 3) {
        
        MAST::SubCellFE*
        local_fe = new MAST::SubCellFE(*_system,
                                       *_intersection);
        // FIXME: we would ideally like to send this to the elem ops object for
        // setting of any local data. But the code has not been setup to do that
        // for SubCellFE. 
        //_elem_ops->set_local_fe_data(*local_fe);

        fe.reset(local_fe);
    }
    else {
        
        fe.reset(new MAST::SubCellFE(*_system,
                                     *_intersection));
    }
    
    return fe;
}

