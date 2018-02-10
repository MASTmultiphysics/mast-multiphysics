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
#include "base/output_assembly_elem_operations.h"
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


MAST::LevelSetIntersection&
MAST::LevelSetNonlinearImplicitAssembly::get_intersection() {
    
    libmesh_assert(_level_set);
    return *_intersection;
}


void
MAST::LevelSetNonlinearImplicitAssembly::
set_level_set_function(MAST::FieldFunction<Real>& level_set) {

    libmesh_assert(!_level_set);
    libmesh_assert(!_intersection);
    
    _level_set    = &level_set;
    _intersection = new MAST::LevelSetIntersection;
}



void
MAST::LevelSetNonlinearImplicitAssembly::clear_level_set_function() {
    
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
    
    MAST::NonlinearImplicitAssemblyElemOperations
    &ops = dynamic_cast<MAST::NonlinearImplicitAssemblyElemOperations&>(*_elem_ops);
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);

        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();

        // Petsc needs that every diagonal term be provided some contribution,
        // even if zero. Otherwise, it complains about lack of diagonal entry.
        // So, if the element is NOT completely on the positive side, we still
        // add a zero matrix to get around this issue.
        if (!_intersection->if_elem_on_positive_phi() && J) {
            
            DenseRealMatrix m(ndofs, ndofs);
            J->add_matrix(m, dof_indices);
        }

        
        if (_intersection->if_elem_has_positive_phi_region()) {
            
            sol.setZero(ndofs);
            vec.setZero(ndofs);
            mat.setZero(ndofs, ndofs);
            
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);

            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                ops.init(*sub_elem);
                ops.set_elem_solution(sol);
                
//                if (_sol_function)
//                    physics_elem->attach_active_solution_function(*_sol_function);
                
                // perform the element level calculations
                ops.elem_calculations(J!=nullptr?true:false,
                                      vec, mat);
                
//                physics_elem->detach_active_solution_function();
                
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

    char s[] = ".";
    //plstring(n, x.data(), y.data(), s);
    plpoin(n, x.data(), y.data(), -1);
}


void
MAST::LevelSetNonlinearImplicitAssembly::plot_sub_elems() {
    
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

        
        plcol0(1); // red color for elements
        plot_elem(*elem);
        plflush();
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);
        
        // now get elements on either side and plot
        const std::vector<const libMesh::Elem *> &
        elems_low = _intersection->get_sub_elems_negative_phi(),
        elems_hi = _intersection->get_sub_elems_positive_phi();
        
        plcol0(15); // white color for sub elements
        
        for (unsigned int i = 0; i < elems_low.size(); i++) {
            plot_elem(*elems_low[i]);
            plflush();
            // create FE
            std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_system, *_intersection));
            fe->init(*elems_low[i]);
            const std::vector<Real>& JxW = fe->get_JxW();
            const std::vector<libMesh::Point>& xyz = fe->get_xyz();
            std::cout << "low: JxW: " << std::accumulate(JxW.begin(),
                                                         JxW.end(), 0.) << std::endl;
            plot_points(xyz);
            plflush();
        }
                
        for (unsigned int i=0; i<elems_hi.size(); i++) {
            plot_elem(*elems_hi[i]);
            plflush();
            
            // create FE
            std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_system, *_intersection));
            fe->init(*elems_hi[i]);
            const std::vector<Real>& JxW = fe->get_JxW();
            const std::vector<libMesh::Point>& xyz = fe->get_xyz();
            std::cout << "hi: JxW: " << std::accumulate(JxW.begin(),
                                                        JxW.end(), 0.) << std::endl;
            plot_points(xyz);
            plflush();
        }

        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    this->clear_elem_operation_object();
}

#endif // HAVE_PLPLOT



bool
MAST::LevelSetNonlinearImplicitAssembly::
sensitivity_assemble (const MAST::FunctionBase& f,
                      libMesh::NumericVector<Real>& sensitivity_rhs) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_elem_ops);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    sensitivity_rhs.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = nonlin_sys.get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                    *nonlin_sys.solution).release());
    
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
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);

        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        if (_intersection->if_elem_has_positive_phi_region()) {
            
            const std::vector<const libMesh::Elem *> &
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                ops.init(*sub_elem);
                ops.set_elem_sensitivity_parameter(f);
                ops.set_elem_solution(sol);
                
                //        if (_sol_function)
                //            physics_elem->attach_active_solution_function(*_sol_function);
                
                // perform the element level calculations
                ops.elem_sensitivity_calculations(vec);
                
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
            }
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
                 MAST::OutputAssemblyElemOperations& output) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);
    this->set_elem_operation_object(output);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                    X).release());
    
    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);
        
        
        if (_intersection->if_elem_has_positive_phi_region()) {
            
            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                dof_map.dof_indices (elem, dof_indices);
                
                
                // get the solution
                unsigned int ndofs = (unsigned int)dof_indices.size();
                sol.setZero(ndofs);
                
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                
                //                if (_sol_function)
                //                    physics_elem->attach_active_solution_function(*_sol_function);
                
                output.init(*sub_elem);
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
    this->clear_elem_operation_object();
}




void
MAST::LevelSetNonlinearImplicitAssembly::
calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                            MAST::OutputAssemblyElemOperations& output,
                            libMesh::NumericVector<Real>& dq_dX) {
    
    libmesh_assert(_discipline);
    libmesh_assert(_system);
    
    this->set_elem_operation_object(output);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    dq_dX.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    
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

        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        if (_intersection->if_elem_has_positive_phi_region()) {
            
            const std::vector<const libMesh::Elem *> &
            elems_hi = _intersection->get_sub_elems_positive_phi();

            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                //        if (_sol_function)
                //            physics_elem->attach_active_solution_function(*_sol_function);
                
                output.init(*sub_elem);
                output.set_elem_solution(sol);
                output.output_derivative_for_elem(vec);
                output.clear_elem();
                
                DenseRealVector v;
                MAST::copy(v, vec);
                dof_map.constrain_element_vector(v, dof_indices);
                dq_dX.add_vector(v, dof_indices);
            }
            //        physics_elem->detach_active_solution_function();
        }
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    dq_dX.close();
    this->clear_elem_operation_object();
}




void
MAST::LevelSetNonlinearImplicitAssembly::
calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                    const libMesh::NumericVector<Real>& dXdp,
                                    const MAST::FunctionBase& p,
                                    MAST::OutputAssemblyElemOperations& output) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);

    this->set_elem_operation_object(output);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.zero();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol,
    dsol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    localized_solution.reset(build_localized_vector(nonlin_sys,
                                                    X).release());
    localized_solution_sens.reset(build_localized_vector(nonlin_sys,
                                                         dXdp).release());

    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        _intersection->init(*_level_set, *elem, nonlin_sys.time);
        
        
        if (_intersection->if_elem_has_positive_phi_region()) {
            
            const std::vector<const libMesh::Elem *> &
            //elems_low = intersect.get_sub_elems_negative_phi(),
            elems_hi = _intersection->get_sub_elems_positive_phi();
            
            std::vector<const libMesh::Elem*>::const_iterator
            hi_sub_elem_it  = elems_hi.begin(),
            hi_sub_elem_end = elems_hi.end();
            
            for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
                
                const libMesh::Elem* sub_elem = *hi_sub_elem_it;
                
                dof_map.dof_indices (elem, dof_indices);
                
                
                // get the solution
                unsigned int ndofs = (unsigned int)dof_indices.size();
                sol.setZero(ndofs);
                
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol(i) = (*localized_solution)(dof_indices[i]);
                
                //                if (_sol_function)
                //                    physics_elem->attach_active_solution_function(*_sol_function);
                
                output.init(*sub_elem);
                output.set_elem_solution(sol);
                output.set_elem_solution_sensitivity(dsol);
                output.evaluate_sensitivity(p);
                output.clear_elem();
            }
        }
        
        _intersection->clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    this->clear_elem_operation_object();
}




void
MAST::LevelSetNonlinearImplicitAssembly::constrain() {

    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(_level_set);

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
    libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    

    std::vector<libMesh::dof_id_type>
    group_A,
    group_B,
    dof_indices,
    constrained_dof_indices;
    constrained_dof_indices.reserve(dof_map.n_local_dofs());
    
    // our intent is to constrain only those dofs that belong to the
    // unintersected elements on the negative phi side of level set, AND
    // if they do not belong to any elements intersected by the level set.
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        MAST::LevelSetIntersection intersect;
        intersect.init(*_level_set, *elem, nonlin_sys.time);
        
        // if the element is entirely on the negative side of the level set,
        // we will constrain all dofs of the element to zero
        
        dof_indices.clear();

        if (intersect.if_elem_on_negative_phi()) {
            // This identifies if the element is entirely on the negative side
            // All of these dofs are constrained

            dof_map.dof_indices(elem, dof_indices);
        }
        else if (intersect.get_intersection_mode() == MAST::THROUGH_NODE &&
                 !intersect.if_elem_has_positive_phi_region()) {
            
            // boundary passes through node and the element is on the negative
            // side. Then, constrain everything but the boundary node.
            
            group_A.clear();
            group_B.clear();
            
            // all dofs on the element
            dof_map.dof_indices(elem, group_A);
            dof_indices.reserve(group_A.size());
            
            // all dofs on the node should not be constrained
            const libMesh::Node
            *node = elem->node_ptr(intersect.node_on_boundary());
            dof_map.dof_indices(node, group_B);
            
            std::set_difference(group_A.begin(),
                                group_A.end(),
                                group_B.begin(),
                                group_B.end(),
                                std::inserter(dof_indices, dof_indices.begin()));
        }
        else if (intersect.get_intersection_mode() == MAST::THROUGH_NODE &&
                 !intersect.if_elem_has_positive_phi_region()) {
            
            // boundary passes through edge and the element is on the negative
            // side. Then, constrain everything but the boundary node.

            group_A.clear();
            group_B.clear();
            
            // all dofs on the element
            dof_map.dof_indices(elem, group_A);
            dof_indices.reserve(group_A.size());
            
            // all dofs on the node, which should not be constraint
            std::unique_ptr<const libMesh::Elem>
            side(elem->side_ptr(intersect.edge_on_boundary()).release());
            dof_map.dof_indices(side.get(), group_B);
            
            std::set_difference(group_A.begin(),
                                group_A.end(),
                                group_B.begin(),
                                group_B.end(),
                                std::inserter(dof_indices, dof_indices.begin()));
        }
     

        constrained_dof_indices.insert(constrained_dof_indices.end(),
                                       dof_indices.begin(),
                                       dof_indices.end());
    }
    
    // create a set so that we only deal with unique set of ids.
    std::set<libMesh::dof_id_type>
    constrained_dofs(constrained_dof_indices.begin(),
                     constrained_dof_indices.end());
    
    // now, constrain everythign in the set
    std::set<libMesh::dof_id_type>::const_iterator
    dof_it  = constrained_dofs.begin(),
    dof_end = constrained_dofs.end();

    for ( ; dof_it != dof_end; dof_it++) {
        
        // if the dof is already Dirichlet constrained, then we do not
        // add another constraint on it
        if (!dof_map.is_constrained_dof(*dof_it)) {
            
            libMesh::DofConstraintRow c_row;
            dof_map.add_constraint_row(*dof_it, c_row, true);
        }
    }
}



std::unique_ptr<MAST::FEBase>
MAST::LevelSetNonlinearImplicitAssembly::build_fe(const libMesh::Elem& elem) {
    
    libmesh_assert(_elem_ops);
    libmesh_assert(_system);
    libmesh_assert(_intersection);
    
    std::unique_ptr<MAST::FEBase> fe;
    
    if (_elem_ops->if_use_local_elem() &&
        elem.dim() < 3) {
        
        MAST::SubCellFE*
        local_fe = new MAST::SubCellFE(*_system, *_intersection);
        // FIXME: we would ideally like to send this to the elem ops object for
        // setting of any local data. But the code has not been setup to do that
        // for SubCellFE. 
        //_elem_ops->set_local_fe_data(*local_fe);

        fe.reset(local_fe);
    }
    else {
        
        fe.reset(new MAST::SubCellFE(*_system, *_intersection));
    }
    
    return fe;
}

