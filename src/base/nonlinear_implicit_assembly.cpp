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
#include "base/nonlinear_implicit_assembly.h"
#include "base/system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "boundary_condition/point_load_condition.h"
#include "numerics/utility.h"
#include "mesh/geom_elem.h"

// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::NonlinearImplicitAssembly::
NonlinearImplicitAssembly():MAST::AssemblyBase(),
_post_assembly           (nullptr),
_res_l2_norm             (0.),
_first_iter_res_l2_norm  (-1.) {
    
}



MAST::NonlinearImplicitAssembly::~NonlinearImplicitAssembly() {
    
}



void
MAST::NonlinearImplicitAssembly::
set_post_assembly_operation(MAST::NonlinearImplicitAssembly::PostAssemblyOperation& post) {
    
    _post_assembly = &post;
}



void
MAST::NonlinearImplicitAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

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
        _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id())) {

            if (J) {
                
                unsigned int ndofs = (unsigned int)dof_indices.size();
                mat.setIdentity(ndofs, ndofs);
                mat *= 1.e-24;
                DenseRealMatrix m;
                MAST::copy(m, mat);
                dof_map.constrain_element_matrix(m, dof_indices);
                J->add_matrix(m, dof_indices);
                dof_indices.clear();
            }
        }
        else {
                        
            MAST::GeomElem geom_elem;
            ops.set_elem_data(elem->dim(), *elem, geom_elem);
            geom_elem.init(*elem, *_system);
            
            ops.init(geom_elem);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            vec.setZero(ndofs);
            mat.setZero(ndofs, ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
            
            ops.set_elem_solution(sol);
            
            
            //        if (_sol_function)
            //            physics_elem->attach_active_solution_function(*_sol_function);
            
            //_check_element_numerical_jacobian(*physics_elem, sol);
            
            // perform the element level calculations
            ops.elem_calculations(J!=nullptr?true:false,
                                  vec, mat);
            
            //        physics_elem->detach_active_solution_function();
            
            ops.clear_elem();
            
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
    }

    
    // add the point loads if any in the discipline
    if (R && _discipline->point_loads().size()) {
        
        const MAST::PointLoadSetType&
        loads = _discipline->point_loads();
        
        vec = RealVectorX::Zero(_system->n_vars());
        
        MAST::PointLoadSetType::const_iterator
        it    = loads.begin(),
        end   = loads.end();
        
        const libMesh::dof_id_type
        first_dof  = dof_map.first_dof(nonlin_sys.comm().rank()),
        end_dof   = dof_map.end_dof(nonlin_sys.comm().rank());

        for ( ; it != end; it++) {
            
            // get the point load function
            const MAST::FieldFunction<RealVectorX>
            &func = (*it)->get<MAST::FieldFunction<RealVectorX>>("load");
            
            // get the nodes on which this object defines the load
            const std::set<const libMesh::Node*>
            nodes = (*it)->get_nodes();
            
            std::set<const libMesh::Node*>::const_iterator
            n_it    = nodes.begin(),
            n_end   = nodes.end();
            
            for (; n_it != n_end; n_it++) {
                
                // load at the node
                vec.setZero();
                func(**n_it, nonlin_sys.time, vec);
                // multiply with -1 to be consistent with res(X) = 0, which
                // requires taking the force vector on RHS to the LHS
                vec *= -1.;
                
                dof_map.dof_indices(*n_it, dof_indices);

                libmesh_assert_equal_to(dof_indices.size(), vec.rows());

                // zero the components of the vector if they do not
                // belong to this processor
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    if (dof_indices[i] <   first_dof  ||
                        dof_indices[i] >=  end_dof)
                        vec(i) = 0.;

                DenseRealVector v;
                MAST::copy(v, vec);
                
                dof_map.constrain_element_vector(v, dof_indices);
                R->add_vector(v, dof_indices);
                dof_indices.clear();
            }
        }
    }
    
    // call the post assembly object, if provided by user
    if (_post_assembly)
        _post_assembly->post_assembly(X, R, J, S);
    

    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) {
        
        R->close();
        _res_l2_norm = R->l2_norm();
        if (_first_iter_res_l2_norm < 0.)
            _first_iter_res_l2_norm = _res_l2_norm;
    }
    if (J && close_matrix) J->close();
}




void
MAST::NonlinearImplicitAssembly::
linearized_jacobian_solution_product (const libMesh::NumericVector<Real>& X,
                                      const libMesh::NumericVector<Real>& dX,
                                      libMesh::NumericVector<Real>& JdX,
                                      libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    // zero the solution vector
    JdX.zero();
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, dsol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_perturbed_solution;
    
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                     X).release());
    localized_perturbed_solution.reset(build_localized_vector(nonlin_sys,
                                                               dX).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;
            
        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dsol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol (i) = (*localized_solution)          (dof_indices[i]);
            dsol(i) = (*localized_perturbed_solution)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        ops.set_elem_perturbed_solution(dsol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        ops.elem_linearized_jacobian_solution_product(vec);
        
        //physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        JdX.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    JdX.close();
}



void
MAST::NonlinearImplicitAssembly::
second_derivative_dot_solution_assembly (const libMesh::NumericVector<Real>& X,
                                         bool if_localize_sol,
                                         const libMesh::NumericVector<Real>& dX,
                                         bool if_localize_sol_sens,
                                         libMesh::SparseMatrix<Real>& d_JdX_dX,
                                         libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    // zero the matrix
    d_JdX_dX.zero();
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol, dsol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    const libMesh::NumericVector<Real>
    *sol_vec  = nullptr,
    *dsol_vec = nullptr;
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_perturbed_solution;
    
    if (if_localize_sol) {
        localized_solution.reset(build_localized_vector(nonlin_sys, X).release());
        sol_vec = localized_solution.get();
    }
    else
        sol_vec = &X;
    
    if (if_localize_sol_sens) {
        localized_perturbed_solution.reset(build_localized_vector(nonlin_sys, dX).release());
        dsol_vec = localized_perturbed_solution.get();
    }
    else
        dsol_vec = &dX;
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;

        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        dsol.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++) {
            sol (i) = (*sol_vec)          (dof_indices[i]);
            dsol(i) = (*dsol_vec)(dof_indices[i]);
        }
        
        ops.set_elem_solution(sol);
        ops.set_elem_solution_sensitivity(dsol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        ops.elem_second_derivative_dot_solution_assembly(mat);
        
//        physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealMatrix m;
        MAST::copy(m, mat);
        
        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_matrix(m, dof_indices);
        
        // add to the global matrices
        d_JdX_dX.add_matrix(m, dof_indices);
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    d_JdX_dX.close();
}





bool
MAST::NonlinearImplicitAssembly::
sensitivity_assemble (const libMesh::NumericVector<Real>& X,
                      bool if_localize_sol,
                      const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs,
                      bool close_vector) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, vec1, sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    const libMesh::NumericVector<Real>
    *sol_vec = nullptr;
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    
    if (if_localize_sol) {
        localized_solution.reset(build_localized_vector(nonlin_sys, X).release());
        sol_vec = localized_solution.get();
    }
    else
        sol_vec = &X;
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( *nonlin_sys.solution, false);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
    
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;

        // no sensitivity computation assembly is neeed in these cases
        if (_param_dependence &&
            // if object is specified and elem does not depend on it
            !_param_dependence->if_elem_depends_on_parameter(*elem, f))
            continue;

        dof_map.dof_indices (elem, dof_indices);
        
        MAST::GeomElem geom_elem;
        ops.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        ops.init(geom_elem);

        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        vec1.setZero(ndofs);

        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*sol_vec)(dof_indices[i]);
        
        ops.set_elem_solution(sol);
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        ops.elem_sensitivity_calculations(f, vec);
        if (f.is_topology_parameter()) {
            ops.elem_topology_sensitivity_calculations(f, vec1);
            vec += vec1;
        }
        
        
//        physics_elem->detach_active_solution_function();
        ops.clear_elem();

        // copy to the libMesh matrix for further processing
        DenseRealVector v;
        MAST::copy(v, vec);

        // constrain the quantities to account for hanging dofs,
        // Dirichlet constraints, etc.
        dof_map.constrain_element_vector(v, dof_indices);
        
        // add to the global matrices
        sensitivity_rhs.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    // add the point loads if any in the discipline
    if (_discipline->point_loads().size()) {
        
        const MAST::PointLoadSetType&
        loads = _discipline->point_loads();
        
        vec = RealVectorX::Zero(_system->n_vars());
        
        MAST::PointLoadSetType::const_iterator
        it    = loads.begin(),
        end   = loads.end();
        
        const libMesh::dof_id_type
        first_dof  = dof_map.first_dof(nonlin_sys.comm().rank()),
        end_dof   = dof_map.end_dof(nonlin_sys.comm().rank());
        
        for ( ; it != end; it++) {
            
            // get the point load function
            const MAST::FieldFunction<RealVectorX>
            &func = (*it)->get<MAST::FieldFunction<RealVectorX>>("load");
            
            // get the nodes on which this object defines the load
            const std::set<const libMesh::Node*>
            nodes = (*it)->get_nodes();
            
            std::set<const libMesh::Node*>::const_iterator
            n_it    = nodes.begin(),
            n_end   = nodes.end();
            
            for (; n_it != n_end; n_it++) {
                
                // load at the node
                vec.setZero();
                func.derivative(f, **n_it, nonlin_sys.time, vec);
                vec *= -1.;
                
                dof_map.dof_indices(*n_it, dof_indices);
                
                libmesh_assert_equal_to(dof_indices.size(), vec.rows());

                // zero the components of the vector if they do not
                // belong to this processor
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    if (dof_indices[i] <   first_dof  ||
                        dof_indices[i] >=  end_dof)
                        vec(i) = 0.;

                DenseRealVector v;
                MAST::copy(v, vec);
                
                dof_map.constrain_element_vector(v, dof_indices);
                sensitivity_rhs.add_vector(v, dof_indices);
                dof_indices.clear();
            }
        }
    }

    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    if (close_vector)
        sensitivity_rhs.close();
    
    return true;
}


