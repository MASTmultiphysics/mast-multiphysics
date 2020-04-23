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
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_intersection.h"
#include "level_set/interface_dof_handler.h"
#include "level_set/level_set_void_solution.h"
#include "level_set/level_set_intersected_elem.h"
#include "level_set/sub_cell_fe.h"
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_implicit_assembly_elem_operations.h"
#include "base/output_assembly_elem_operations.h"
#include "base/elem_base.h"
#include "numerics/utility.h"


// libMesh includes
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"



MAST::LevelSetNonlinearImplicitAssembly::
LevelSetNonlinearImplicitAssembly(bool enable_dof_handler):
MAST::NonlinearImplicitAssembly(),
_enable_dof_handler              (enable_dof_handler),
_evaluate_output_on_negative_phi (false),
_level_set                       (nullptr),
_indicator                       (nullptr),
_intersection                    (nullptr),
_dof_handler                     (nullptr),
_void_solution_monitor           (nullptr),
_velocity                        (nullptr),
_filter                          (nullptr) {
    
}



MAST::LevelSetNonlinearImplicitAssembly::~LevelSetNonlinearImplicitAssembly() {
 
    if (_intersection)          delete _intersection;
    if (_dof_handler)           delete _dof_handler;
    if (_void_solution_monitor) delete _void_solution_monitor;
}


MAST::LevelSetIntersection&
MAST::LevelSetNonlinearImplicitAssembly::get_intersection() {
    
    libmesh_assert(_level_set);
    return *_intersection;
}



void
MAST::LevelSetNonlinearImplicitAssembly::set_evaluate_output_on_negative_phi(bool f) {
 
    _evaluate_output_on_negative_phi = f;
}


bool
MAST::LevelSetNonlinearImplicitAssembly::if_use_dof_handler() const {
    return _enable_dof_handler;
}


MAST::LevelSetInterfaceDofHandler&
MAST::LevelSetNonlinearImplicitAssembly::get_dof_handler() {
    
    libmesh_assert(_level_set);
    libmesh_assert(_dof_handler);
    return *_dof_handler;
}


void
MAST::LevelSetNonlinearImplicitAssembly::
set_level_set_function(MAST::FieldFunction<Real>& level_set,
                       const MAST::FilterBase& filter) {

    libmesh_assert(!_level_set);
    libmesh_assert(!_intersection);
    libmesh_assert(_system);
    
    _level_set    = &level_set;
    _filter       = &filter;
    _intersection = new MAST::LevelSetIntersection();
    if (_enable_dof_handler) {
        _dof_handler  = new MAST::LevelSetInterfaceDofHandler();
        _dof_handler->init(*_system, *_intersection, *_level_set);
        _void_solution_monitor = new MAST::LevelSetVoidSolution();
        _void_solution_monitor->init(*this, *_intersection, *_dof_handler);
    }
}



void
MAST::LevelSetNonlinearImplicitAssembly::
set_indicator_function(MAST::FieldFunction<RealVectorX>& indicator) {
    
    libmesh_assert(!_indicator);
    libmesh_assert(_system);
    
    _indicator    = &indicator;
}



void
MAST::LevelSetNonlinearImplicitAssembly::clear_level_set_function() {
    
    _level_set = nullptr;
    _filter    = nullptr;
    
    if (_intersection) {
        delete _intersection;
        delete _dof_handler;
        delete _void_solution_monitor;
        
        _intersection          = nullptr;
        _dof_handler           = nullptr;
        _void_solution_monitor = nullptr;
    }
}




void
MAST::LevelSetNonlinearImplicitAssembly::
set_level_set_velocity_function(MAST::FieldFunction<RealVectorX>& velocity) {
    
    libmesh_assert(_level_set);
    libmesh_assert(_intersection);
    libmesh_assert(_system);
    libmesh_assert(!_velocity);
    
    _velocity = &velocity;
}



void
MAST::LevelSetNonlinearImplicitAssembly::clear_level_set_velocity_function() {

    _velocity = nullptr;
}



#if MAST_ENABLE_PLPLOT == 1

#include <numeric>
#include "plplot/plplot.h"

void plot_elem(const libMesh::Elem& e) {

    unsigned int
    n  = e.n_nodes();

    RealVectorX
    x  = RealVectorX::Zero(n+1),
    y  = RealVectorX::Zero(n+1);

    for (unsigned int i=0; i<n+1; i++) {
        x(i) = e.point(i%n)(0);
        y(i) = e.point(i%n)(1);
    }

    plline(n+1, x.data(), y.data());
    plflush();
}


void plot_points(const std::vector<libMesh::Point>& pts) {

    unsigned int
    n  = pts.size();

    RealVectorX
    x  = RealVectorX::Zero(n),
    y  = RealVectorX::Zero(n);

    for (unsigned int i=0; i<n; i++) {
        x(i) = pts[i](0);
        y(i) = pts[i](1);
    }

    plpoin(n, x.data(), y.data(), -1);
    plflush();
}


void
MAST::LevelSetNonlinearImplicitAssembly::plot_sub_elems(bool plot_reference_elem,
                                                        bool plot_low_phi_elem,
                                                        bool plot_high_phi_elem) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    plsdev("xwin");
    plinit();
    plenv(0,.3,0,.3,0,0);
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        {
            MAST::FEBase fe(*_system);
            fe.init(*elem);
            const std::vector<Real>& JxW = fe.get_JxW();
            std::cout << "==== original JxW: " << std::accumulate(JxW.begin(),
                                                                  JxW.end(), 0.) << std::endl;
        }


        if (plot_reference_elem) {
            plcol0(1); // red color for elements
            plot_elem(*elem);
        }
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);
        
        // now get elements on either side and plot
        const std::vector<const libMesh::Elem *> &
        elems_low = _intersection->get_sub_elems_negative_phi(),
        elems_hi = _intersection->get_sub_elems_positive_phi();
        
        
        if (plot_low_phi_elem) {
            
            for (unsigned int i = 0; i < elems_low.size(); i++) {
                
                plcol0(3); // green color for sub elements
                plot_elem(*elems_low[i]);
        
                // create FE
                std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_system, *_intersection));
                fe->init(*elems_low[i]);
                const std::vector<Real>& JxW = fe->get_JxW();
                const std::vector<libMesh::Point>& xyz = fe->get_xyz();
                std::cout << "low: JxW: " << std::accumulate(JxW.begin(),
                                                             JxW.end(), 0.) << std::endl;
                plot_points(xyz);
            }
        }
        
        
        if (plot_high_phi_elem) {
            
            for (unsigned int i=0; i<elems_hi.size(); i++) {
                
                plcol0(15); // white color for sub elements
                plot_elem(*elems_hi[i]);
                
                // create FE
                std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_system, *_intersection));
                fe->init(*elems_hi[i]);
                const std::vector<Real>& JxW = fe->get_JxW();
                const std::vector<libMesh::Point>& xyz = fe->get_xyz();
                std::cout << "hi: JxW: " << std::accumulate(JxW.begin(),
                                                            JxW.end(), 0.) << std::endl;
                plot_points(xyz);
            }
        }
        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}

#endif // HAVE_PLPLOT




void
MAST::LevelSetNonlinearImplicitAssembly::
residual_and_jacobian (const libMesh::NumericVector<Real>& X,
                       libMesh::NumericVector<Real>* R,
                       libMesh::SparseMatrix<Real>*  J,
                       libMesh::NonlinearImplicitSystem& S) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);
    libmesh_assert(_level_set);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    // make sure that the system for which this object was created,
    // and the system passed through the function call are the same
    libmesh_assert_equal_to(&S, &(nonlin_sys));
    
    if (R) R->zero();
    if (J) J->zero();
    
    const Real
    tol   = 1.e-10;
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    vec,
    sol,
    sub_elem_vec,
    res_factored_u,
    nd_indicator = RealVectorX::Ones(1),
    indicator    = RealVectorX::Zero(1);
    RealMatrixX
    mat,
    sub_elem_mat,
    jac_factored_uu;
    
    std::vector<libMesh::dof_id_type>
    dof_indices,
    material_rows;
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
    
    MAST::NonlinearImplicitAssemblyElemOperations
    &ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time,
                            nonlin_sys.get_mesh().max_elem_id(),
                            nonlin_sys.get_mesh().max_node_id());
        dof_map.dof_indices (elem, dof_indices);
        
        // use the indicator if it was provided
        if (_indicator) {
            
            nd_indicator.setZero(elem->n_nodes());
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                (*_indicator)(elem->node_ref(i), nonlin_sys.time, indicator);
                nd_indicator(i) = indicator(0);
            }
        }

        unsigned int ndofs = (unsigned int)dof_indices.size();

        // Petsc needs that every diagonal term be provided some contribution,
        // even if zero. Otherwise, it complains about lack of diagonal entry.
        // So, if the element is NOT completely on the positive side, we still
        // add a zero matrix to get around this issue.
        if ((_intersection->if_elem_on_negative_phi() ||
             nd_indicator.maxCoeff() < tol) && J) {
            
            DenseRealMatrix m(ndofs, ndofs);
            //dof_map.constrain_element_matrix(m, dof_indices);
            for (unsigned int i=0; i<ndofs; i++)
                m(i,i) = 1.e-14;
            dof_map.constrain_element_matrix(m, dof_indices);
            J->add_matrix(m, dof_indices);
        }
        
        
        if (nd_indicator.maxCoeff() > tol &&
            _intersection->if_elem_has_positive_phi_region()) {
            
            // get the solution
            sol.setZero(ndofs);
            vec.setZero(ndofs);
            sub_elem_vec.setZero(ndofs);
            mat.setZero(ndofs, ndofs);
            sub_elem_mat.setZero(ndofs, ndofs);

            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
                _dof_handler->solution_of_factored_element(*elem, sol);
            
            // if the element has been marked for factorization,
            // get the factorized jacobian and residual contributions
            if (_dof_handler && _dof_handler->if_factor_element(*elem)) {
                
                // the Jacobian is based on the homogenizaton method to maintain
                // a well conditioned global Jacobian.
                MAST::GeomElem geom_elem;
                ops.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*elem, *_system);
                
                ops.init(geom_elem);
                ops.set_elem_solution(sol);
                ops.elem_calculations(true, vec, mat);
                ops.clear_elem();
                mat *= _intersection->get_positive_phi_volume_fraction();
                // the residual based on homogenization is used for factorized
                // elements since the factorization depends on exact Jacobian
                // and the homogenization based Jacobian is not exact linearization
                // of residual based on sub cells.
                vec *= _intersection->get_positive_phi_volume_fraction();

                _dof_handler->element_factored_residual_and_jacobian(*elem,
                                                                     mat,
                                                                     vec,
                                                                     material_rows,
                                                                     jac_factored_uu,
                                                                     res_factored_u);
                
                // zero the element matrix and update with the
                // factored matrix
                mat.setIdentity(); mat *= 1.e-6; // set a small value on the diagonal to avoid nans
                vec.setZero();
                
                for (unsigned int i=0; i<material_rows.size(); i++) {
                    
                    vec(material_rows[i])   = res_factored_u(i);
                    
                    for (unsigned int j=0; j<material_rows.size(); j++)
                        mat(material_rows[i], material_rows[j]) = jac_factored_uu(i,j);
                }
            }
            else {
                
                // get the residual from the sub elements
                mat.setZero();
                vec.setZero();

                const std::vector<const libMesh::Elem *> &
                elems_hi = _intersection->get_sub_elems_positive_phi();
                
                std::vector<const libMesh::Elem*>::const_iterator
                hi_sub_elem_it  = elems_hi.begin(),
                hi_sub_elem_end = elems_hi.end();
                
                for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                    
                    const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                    
                    MAST::LevelSetIntersectedElem geom_elem;
                    ops.set_elem_data(elem->dim(), *elem, geom_elem);
                    geom_elem.init(*sub_elem, *_system, *_intersection);
                    
                    ops.init(geom_elem);
                    ops.set_elem_solution(sol);
                    
                    ops.elem_calculations(J!=nullptr?true:false, sub_elem_vec, sub_elem_mat);
                    
                    mat += sub_elem_mat;
                    vec += sub_elem_vec;
                    
                    ops.clear_elem();
                }
            }
            
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
        
        _intersection->clear();
    }
    
    // call the post assembly object, if provided by user
    if (_post_assembly)
        _post_assembly->post_assembly(X, R, J, S);
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    if (R) R->close();
    if (J && close_matrix) J->close();
}


bool
MAST::LevelSetNonlinearImplicitAssembly::
sensitivity_assemble (const libMesh::NumericVector<Real>& X,
                      bool if_localize_sol,
                      const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);
    // we need the velocity for topology parameter
    if (f.is_topology_parameter()) libmesh_assert(_velocity);

    // we need the velocity for topology parameter
    const MAST::LevelSetParameter
    *p_ls = nullptr;
    if (f.is_topology_parameter()) {
        libmesh_assert(_velocity);
        p_ls = dynamic_cast<const MAST::LevelSetParameter*>(&f);
    }

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    const Real
    tol = 1.e-10;

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    vec1,
    vec2,
    vec_total,
    sol,
    res_factored_u,
    nd_indicator = RealVectorX::Ones(1),
    indicator    = RealVectorX::Zero(1);
    RealMatrixX
    mat,
    jac_factored_uu;

    std::vector<libMesh::dof_id_type>
    dof_indices,
    material_rows;
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
        _sol_function->init( *nonlin_sys.solution);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        // no sensitivity computation assembly is neeed in these cases
        if (_param_dependence &&
            // if object is specified and elem does not depend on it
            !_param_dependence->if_elem_depends_on_parameter(*elem, f))
            continue;

        _intersection->init(*_level_set, *elem, nonlin_sys.time,
                            nonlin_sys.get_mesh().max_elem_id(),
                            nonlin_sys.get_mesh().max_node_id());

        dof_map.dof_indices (elem, dof_indices);
        
        // use the indicator if it was provided
        if (_indicator) {
            
            nd_indicator.setZero(elem->n_nodes());
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                (*_indicator)(elem->node_ref(i), nonlin_sys.time, indicator);
                nd_indicator(i) = indicator(0);
            }
        }

        if (nd_indicator.maxCoeff() > tol &&
            _intersection->if_elem_has_positive_phi_region()) {

            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            vec_total.setZero(ndofs);
            vec1.setZero(ndofs);
            vec2.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*sol_vec)(dof_indices[i]);

            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
                _dof_handler->solution_of_factored_element(*elem, sol);

            const std::vector<const libMesh::Elem *> &
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                MAST::LevelSetIntersectedElem geom_elem;
                ops.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*sub_elem, *_system, *_intersection);
                
                ops.init(geom_elem);
                ops.set_elem_solution(sol);
                
                //        if (_sol_function)
                //            physics_elem->attach_active_solution_function(*_sol_function);
                
                // perform the element level calculations
                ops.elem_sensitivity_calculations(f, vec1);
                
                // if the quantity is also defined as a topology parameter,
                // then calculate sensitivity from boundary movement.
                if (f.is_topology_parameter()) {
                    
                    ops.elem_topology_sensitivity_calculations(f,
                                                               *_velocity,
                                                               vec2);
                    vec1 += vec2;
                }
                ops.clear_elem();
                
                // if the element has been identified for factoring then
                // we will factor the residual and replace it in the
                // residual vector. For this we need to compute the
                // element Jacobian matrix
                if (_dof_handler && _dof_handler->if_factor_element(*elem)) {

                    vec2.setZero(ndofs);
                    mat.setZero(ndofs, ndofs);
                    
                    MAST::GeomElem geom_elem;
                    ops.set_elem_data(elem->dim(), *elem, geom_elem);
                    geom_elem.init(*elem, *_system);
                    
                    ops.init(geom_elem);
                    ops.set_elem_solution(sol);
                    ops.elem_calculations(true, vec2, mat);
                    ops.clear_elem();
                    mat *= _intersection->get_positive_phi_volume_fraction();

                    _dof_handler->element_factored_residual_and_jacobian(*elem,
                                                                         mat,
                                                                         vec1,
                                                                         material_rows,
                                                                         jac_factored_uu,
                                                                         res_factored_u);

                    vec1.setZero();
                    
                    for (unsigned int i=0; i<material_rows.size(); i++)
                        vec1(material_rows[i])   = res_factored_u(i);
                }

                vec_total += vec1;
                

                //        physics_elem->detach_active_solution_function();
                ops.clear_elem();
            }

            // copy to the libMesh matrix for further processing
            DenseRealVector v;
            MAST::copy(v, vec_total);
            
            // constrain the quantities to account for hanging dofs,
            // Dirichlet constraints, etc.
            dof_map.constrain_element_vector(v, dof_indices);
            
            // add to the global matrices
            sensitivity_rhs.add_vector(v, dof_indices);
            dof_indices.clear();
        }

        _intersection->clear();
    }
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->clear();
    
    sensitivity_rhs.close();
    
    return true;
}




void
MAST::LevelSetNonlinearImplicitAssembly::
calculate_output(const libMesh::NumericVector<Real>& X,
                 bool if_localize_sol,
                 MAST::OutputAssemblyElemOperations& output) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);
    
    output.zero_for_analysis();

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.set_assembly(*this);
    
    const Real
    tol   = 1.e-10;

    RealVectorX
    sol,
    nd_indicator = RealVectorX::Ones(1),
    indicator    = RealVectorX::Zero(1);
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
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
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        // use the indicator if it was provided
        if (_indicator) {
            
            nd_indicator.setZero(elem->n_nodes());
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                (*_indicator)(elem->node_ref(i), nonlin_sys.time, indicator);
                nd_indicator(i) = indicator(0);
            }
        }
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time,
                            nonlin_sys.get_mesh().max_elem_id(),
                            nonlin_sys.get_mesh().max_node_id());
        
        if (_evaluate_output_on_negative_phi &&
            _intersection->get_sub_elems_negative_phi().size() ) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*sol_vec)(dof_indices[i]);
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
            _dof_handler->solution_of_factored_element(*elem, sol);
            
            const std::vector<const libMesh::Elem *> &
            elems_low = _intersection->get_sub_elems_negative_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            low_sub_elem_it  = elems_low.begin(),
            low_sub_elem_end = elems_low.end();
            
            for (; low_sub_elem_it != low_sub_elem_end; low_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *low_sub_elem_it;
                
                //                if (_sol_function)
                //                    physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_elem;
                output.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_elem);
                output.set_elem_solution(sol);
                output.evaluate();
                output.clear_elem();
            }
        }

        if (nd_indicator.maxCoeff() > tol &&
            _intersection->if_elem_has_positive_phi_region()) {

            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*sol_vec)(dof_indices[i]);

            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
                _dof_handler->solution_of_factored_element(*elem, sol);

            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                //                if (_sol_function)
                //                    physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_elem;
                output.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_elem);
                output.set_elem_solution(sol);
                output.evaluate();
                output.clear_elem();
            }
        }
        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    output.clear_assembly();
}




void
MAST::LevelSetNonlinearImplicitAssembly::
calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                            bool if_localize_sol,
                            MAST::OutputAssemblyElemOperations& output,
                            libMesh::NumericVector<Real>& dq_dX) {
    
    libmesh_assert(_discipline);
    libmesh_assert(_system);
    
    output.zero_for_sensitivity();

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.set_assembly(*this);
    
    dq_dX.zero();

    const Real
    tol   = 1.e-10;

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    vec1,
    vec2,
    vec_total,
    sol,
    res_factored_u,
    nd_indicator = RealVectorX::Ones(1),
    indicator    = RealVectorX::Zero(1);
    RealMatrixX
    mat,
    jac_factored_uu;

    std::vector<libMesh::dof_id_type>
    dof_indices,
    material_rows;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    MAST::NonlinearImplicitAssemblyElemOperations&
    ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);

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
        _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;

        _intersection->init(*_level_set, *elem, nonlin_sys.time,
                            nonlin_sys.get_mesh().max_elem_id(),
                            nonlin_sys.get_mesh().max_node_id());

        // use the indicator if it was provided
        if (_indicator) {
            
            nd_indicator.setZero(elem->n_nodes());
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                (*_indicator)(elem->node_ref(i), nonlin_sys.time, indicator);
                nd_indicator(i) = indicator(0);
            }
        }


        if (_evaluate_output_on_negative_phi &&
            _intersection->get_sub_elems_negative_phi().size()) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            vec1.setZero(ndofs);
            vec2.setZero(ndofs);
            vec_total.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*sol_vec)(dof_indices[i]);
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
            _dof_handler->solution_of_factored_element(*elem, sol);
            
            const std::vector<const libMesh::Elem *> &
            elems_low = _intersection->get_sub_elems_negative_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            low_sub_elem_it  = elems_low.begin(),
            low_sub_elem_end = elems_low.end();
            
            for (; low_sub_elem_it != low_sub_elem_end; low_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *low_sub_elem_it;
                
                //        if (_sol_function)
                //            physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_sub_elem;
                output.set_elem_data(elem->dim(), *elem, geom_sub_elem);
                geom_sub_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_sub_elem);
                output.set_elem_solution(sol);
                output.output_derivative_for_elem(vec1);
                output.clear_elem();
                
                vec_total += vec1;
            }
            
            DenseRealVector v;
            MAST::copy(v, vec_total);
            dof_map.constrain_element_vector(v, dof_indices);
            dq_dX.add_vector(v, dof_indices);
            dof_indices.clear();
        }


        if (nd_indicator.maxCoeff() > tol &&
            _intersection->if_elem_has_positive_phi_region()) {

            dof_map.dof_indices (elem, dof_indices);
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            vec1.setZero(ndofs);
            vec2.setZero(ndofs);
            vec_total.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*sol_vec)(dof_indices[i]);
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
                _dof_handler->solution_of_factored_element(*elem, sol);

            const std::vector<const libMesh::Elem *> &
            elems_hi = _intersection->get_sub_elems_positive_phi();

            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                //        if (_sol_function)
                //            physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_sub_elem;
                output.set_elem_data(elem->dim(), *elem, geom_sub_elem);
                geom_sub_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_sub_elem);
                output.set_elem_solution(sol);
                output.output_derivative_for_elem(vec1);
                output.clear_elem();

                if (_dof_handler && _dof_handler->if_factor_element(*elem)) {
                    
                    vec2.setZero(ndofs);
                    mat.setZero(ndofs, ndofs);
                    
                    MAST::GeomElem geom_elem;
                    ops.set_elem_data(elem->dim(), *elem, geom_elem);
                    geom_elem.init(*elem, *_system);
                    
                    ops.init(geom_elem);
                    ops.set_elem_solution(sol);
                    ops.elem_calculations(true, vec2, mat);
                    ops.clear_elem();
                    mat *= _intersection->get_positive_phi_volume_fraction();
                    
                    _dof_handler->element_factored_residual_and_jacobian(*elem,
                                                                         mat.transpose(),
                                                                         vec1,
                                                                         material_rows,
                                                                         jac_factored_uu,
                                                                         res_factored_u);
                    
                    vec1.setZero();
                    
                    for (unsigned int i=0; i<material_rows.size(); i++)
                        vec1(material_rows[i])   = res_factored_u(i);
                }

                vec_total += vec1;
            }

            DenseRealVector v;
            MAST::copy(v, vec_total);
            dof_map.constrain_element_vector(v, dof_indices);
            dq_dX.add_vector(v, dof_indices);
            dof_indices.clear();
        }
        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    dq_dX.close();
    output.clear_assembly();
}




void
MAST::LevelSetNonlinearImplicitAssembly::
calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                    bool if_localize_sol,
                                    const libMesh::NumericVector<Real>* dXdp,
                                    bool if_localize_sol_sens,
                                    const MAST::FunctionBase& p,
                                    MAST::OutputAssemblyElemOperations& output) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);
    
    // we need the velocity for topology parameter
    const MAST::LevelSetParameter
    *p_ls = nullptr;
    if (p.is_topology_parameter()) {
        libmesh_assert(_velocity);
        p_ls = dynamic_cast<const MAST::LevelSetParameter*>(&p);
    }

    output.zero_for_sensitivity();

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.set_assembly(*this);
    
    const Real
    tol   = 1.e-10;

    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol,
    dsol,
    nd_indicator = RealVectorX::Ones(1),
    indicator    = RealVectorX::Zero(1);
    
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
        localized_perturbed_solution.reset(build_localized_vector(nonlin_sys, *dXdp).release());
        dsol_vec = localized_perturbed_solution.get();
    }
    else
        dsol_vec = dXdp;

    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        // no sensitivity computation assembly is neeed in these cases
        if (_param_dependence &&
            // if object is specified and elem does not depend on it
            !_param_dependence->if_elem_depends_on_parameter(*elem, p) &&
            // and if no sol_sens is given
            (!dXdp ||
             // or if it can be ignored for elem
             (dXdp && _param_dependence->override_flag)))
            continue;
        
        // use the indicator if it was provided
        if (_indicator) {
            
            nd_indicator.setZero(elem->n_nodes());
            for (unsigned int i=0; i<elem->n_nodes(); i++) {
                (*_indicator)(elem->node_ref(i), nonlin_sys.time, indicator);
                nd_indicator(i) = indicator(0);
            }
        }

        _intersection->init(*_level_set, *elem, nonlin_sys.time,
                            nonlin_sys.get_mesh().max_elem_id(),
                            nonlin_sys.get_mesh().max_node_id());
        
        if (_evaluate_output_on_negative_phi &&
            _intersection->get_sub_elems_negative_phi().size()) {
            
            dof_map.dof_indices (elem, dof_indices);
            
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            dsol.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                sol(i)  = (*sol_vec)(dof_indices[i]);
                if (dXdp)
                    dsol(i) = (*dsol_vec)(dof_indices[i]);
            }
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
            _dof_handler->solution_of_factored_element(*elem, sol);
            
            const std::vector<const libMesh::Elem *> &
            elems_low = _intersection->get_sub_elems_negative_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            low_sub_elem_it  = elems_low.begin(),
            low_sub_elem_end = elems_low.end();
            
            for (; low_sub_elem_it != low_sub_elem_end; low_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *low_sub_elem_it;
                
                // if (_sol_function)
                //   physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_elem;
                output.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_elem);
                output.set_elem_solution(sol);
                output.set_elem_solution_sensitivity(dsol);
                output.evaluate_sensitivity(p);
                if (p.is_topology_parameter())
                output.evaluate_topology_sensitivity(p, *_velocity);
                
                output.clear_elem();
            }
        }


        if (nd_indicator.maxCoeff() > tol &&
            _intersection->if_elem_has_positive_phi_region()) {

            dof_map.dof_indices (elem, dof_indices);
            
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            dsol.setZero(ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++) {
                sol(i)  = (*sol_vec)(dof_indices[i]);
                if (dXdp)
                    dsol(i) = (*dsol_vec)(dof_indices[i]);
            }
            
            // if the element has been marked for factorization then
            // get the void solution from the storage
            if (_dof_handler && _dof_handler->if_factor_element(*elem))
                _dof_handler->solution_of_factored_element(*elem, sol);

            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                // if (_sol_function)
                //   physics_elem->attach_active_solution_function(*_sol_function);
                MAST::LevelSetIntersectedElem geom_elem;
                output.set_elem_data(elem->dim(), *elem, geom_elem);
                geom_elem.init(*sub_elem, *_system, *_intersection);
                
                output.init(geom_elem);
                output.set_elem_solution(sol);
                output.set_elem_solution_sensitivity(dsol);
                output.evaluate_sensitivity(p);
                if (p.is_topology_parameter())
                    output.evaluate_topology_sensitivity(p, *_velocity);
                
                output.clear_elem();
            }
        }
        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    output.clear_assembly();
}

