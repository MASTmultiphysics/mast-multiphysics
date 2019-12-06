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
#include "base/assembly_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_elem_operation.h"
#include "base/output_assembly_elem_operations.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"
#include "numerics/utility.h"
#include "boundary_condition/point_load_condition.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"


MAST::AssemblyBase::AssemblyBase():
close_matrix      (true),
_elem_ops         (nullptr),
_discipline       (nullptr),
_system           (nullptr),
_sol_function     (nullptr),
_solver_monitor   (nullptr),
_param_dependence (nullptr) {
    
}




MAST::AssemblyBase::~AssemblyBase() {
    
    this->clear_discipline_and_system();
    this->clear_elem_operation_object();
}



const MAST::PhysicsDisciplineBase&
MAST::AssemblyBase::discipline() const {
    
    libmesh_assert_msg(_discipline,
                       "Error: Discipline not yet attached to Assembly.");
    return *_discipline;
}



MAST::PhysicsDisciplineBase&
MAST::AssemblyBase::discipline() {
    
    libmesh_assert_msg(_discipline,
                       "Error: Discipline not yet attached to Assembly.");
    return *_discipline;
}




const MAST::NonlinearSystem&
MAST::AssemblyBase::system() const {
    
    libmesh_assert_msg(_discipline,
                       "Error: System not yet attached to Assembly.");
    return _system->system();
}


MAST::NonlinearSystem&
MAST::AssemblyBase::system() {
    
    libmesh_assert_msg(_discipline,
                       "Error: System not yet attached to Assembly.");
    return _system->system();
}


MAST::AssemblyElemOperations&
MAST::AssemblyBase::get_elem_ops() {
    
    libmesh_assert_msg(_elem_ops,
                       "Error: Not yet initialized.");
    return *_elem_ops;
}


MAST::SystemInitialization&
MAST::AssemblyBase::system_init()  {
    
    libmesh_assert_msg(_system,
                       "Error: System not yet attached to Assembly.");
    return *_system;
}


void
MAST::AssemblyBase::set_solver_monitor(MAST::AssemblyBase::SolverMonitor& monitor) {
    
    libmesh_assert(!_solver_monitor);
    _solver_monitor = &monitor;
}



void
MAST::AssemblyBase::attach_elem_parameter_dependence_object
(MAST::AssemblyBase::ElemParameterDependence& dep) {

    libmesh_assert(!_param_dependence);
    
    _param_dependence = &dep;
}



void
MAST::AssemblyBase::clear_elem_parameter_dependence_object() {
    
    _param_dependence = nullptr;
}



MAST::AssemblyBase::SolverMonitor*
MAST::AssemblyBase::get_solver_monitor() {
    
    return _solver_monitor;
}


void
MAST::AssemblyBase::clear_solver_monitor() {
    
    _solver_monitor = nullptr;
}


void
MAST::AssemblyBase::
set_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                          MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
}



void
MAST::AssemblyBase::
clear_discipline_and_system() {
    
    close_matrix      = true;
    _discipline       = nullptr;
    _system           = nullptr;
    _param_dependence = nullptr;
}



void
MAST::AssemblyBase::set_elem_operation_object(MAST::AssemblyElemOperations& elem_ops) {
    
    libmesh_assert_msg(!_elem_ops,
                       "Error: Assembly should be cleared before attaching Elem.");
    
    _elem_ops   = &elem_ops;
    _elem_ops->set_assembly(*this);
}



void
MAST::AssemblyBase::clear_elem_operation_object() {

    if (_elem_ops) {
        _elem_ops->clear_assembly();
        _elem_ops    = nullptr;
    }
}



std::unique_ptr<libMesh::NumericVector<Real> >
MAST::AssemblyBase::
build_localized_vector(const libMesh::System& sys,
                       const libMesh::NumericVector<Real>& global) const {
    
    libMesh::NumericVector<Real>* local =
    libMesh::NumericVector<Real>::build(sys.comm()).release();
    
    const std::vector<libMesh::dof_id_type>& send_list =
    sys.get_dof_map().get_send_list();
    
    local->init(sys.n_dofs(),
                sys.n_local_dofs(),
                send_list,
                false,
                libMesh::GHOSTED);
    global.localize(*local, send_list);
    
    return std::unique_ptr<libMesh::NumericVector<Real> >(local);
}




void
MAST::AssemblyBase::attach_solution_function(MAST::MeshFieldFunction& f){
    
    // make sure that no prior association is specified
    libmesh_assert(!_sol_function);
    
    _sol_function = &f;
}




void
MAST::AssemblyBase::detach_solution_function() {
    _sol_function = nullptr;
}




void
MAST::AssemblyBase::calculate_output(const libMesh::NumericVector<Real>& X,
                                     bool if_localize_sol,
                                     MAST::OutputAssemblyElemOperations& output) {
    
    libmesh_assert(_discipline);
    libmesh_assert(_system);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.set_assembly(*this);
    
    output.zero_for_analysis();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    const libMesh::NumericVector<Real>*
    sol_vec = nullptr;
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    if (if_localize_sol) {
        localized_solution.reset(build_localized_vector(nonlin_sys, X).release());
        sol_vec = localized_solution.get();
    }
    else
        sol_vec = &X;
    
    
    // if a solution function is attached, initialize it
    //if (_sol_function)
    //    _sol_function->init( X, false);
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    

    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;

        //if (_sol_function)
        //    physics_elem->attach_active_solution_function(*_sol_function);
        
        MAST::GeomElem geom_elem;
        output.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);
        
        if (!output.if_evaluate_for_element(geom_elem)) continue;
        
        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*sol_vec)(dof_indices[i]);
        
        output.init(geom_elem);
        output.set_elem_solution(sol);
        output.evaluate();
        output.clear_elem();

        //physics_elem->detach_active_solution_function();
    }
    
    /**
     * Now we need to iterate over node for node based operations such as 
     * the compliance due to point loads, a weighted sum of displacements,
     * etc.
     */
    
    // If point loads exist, assumble a vector of them
    std::unique_ptr<libMesh::NumericVector<Real>> Fp = 
        nonlin_sys.solution->zero_clone();
    calculate_point_load_vector(*Fp);
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_Fp;
    localized_Fp.reset(build_localized_vector(nonlin_sys, *Fp).release());
    
    /**
     *  Now that we have point loads, we can now pass this information along
     *  with the displacements to the output object.
     */
    libMesh::MeshBase::node_iterator nd = 
        nonlin_sys.get_mesh().local_nodes_begin();
        
    const libMesh::MeshBase::node_iterator end_nd = 
        nonlin_sys.get_mesh().local_nodes_end();
    
    for ( ; nd != end_nd; ++nd)
    {
        libMesh::Node* node = *nd;
        
        // Get the global DOFs belonging to this node
        dof_map.dof_indices(node, dof_indices);
        
        // get the solution for this node
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        RealVectorX fpt = RealVectorX::Zero(ndofs);
        for (unsigned int i=0; i<dof_indices.size(); i++)
        {
            sol(i) = (*localized_solution)(dof_indices[i]);
            fpt(i) = (*localized_Fp)(dof_indices[i]);
        }
        
        output.set_node_data(node);
        output.evaluate_for_node(sol, fpt);
        output.clear_node();
        dof_indices.clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    output.clear_assembly();
}




void
MAST::AssemblyBase::
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
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    const libMesh::NumericVector<Real>*
    sol_vec = nullptr;
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    
    if (if_localize_sol) {
        localized_solution.reset(build_localized_vector(nonlin_sys, X).release());
        sol_vec = localized_solution.get();
    }
    else
        sol_vec = &X;
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;

        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;
        
        MAST::GeomElem geom_elem;
        output.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);

        if (!output.if_evaluate_for_element(geom_elem)) continue;

        dof_map.dof_indices (elem, dof_indices);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*sol_vec)(dof_indices[i]);
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);

        output.init(geom_elem);
        output.set_elem_solution(sol);
        output.output_derivative_for_elem(vec);
        output.clear_elem();
        
        DenseRealVector v;
        MAST::copy(v, vec);
        dof_map.constrain_element_vector(v, dof_indices);
        dq_dX.add_vector(v, dof_indices);
        dof_indices.clear();
    }
    
    /**
     * Now we need to iterate over node for node based operations such as 
     * the compliance due to point loads, a weighted sum of displacements,
     * etc.
     */
    
    // If point loads exist, assumble a vector of them
    std::unique_ptr<libMesh::NumericVector<Real>> Fp = 
        nonlin_sys.solution->zero_clone();
    calculate_point_load_vector(*Fp);
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_Fp;
    localized_Fp.reset(build_localized_vector(nonlin_sys, *Fp).release());
    
    /**
     *  Now that we have point loads, we can now pass this information along
     *  with the displacements to the output object.
     */
    libMesh::MeshBase::node_iterator nd = 
        nonlin_sys.get_mesh().local_nodes_begin();
        
    const libMesh::MeshBase::node_iterator end_nd = 
        nonlin_sys.get_mesh().local_nodes_end();
    
    for ( ; nd != end_nd; ++nd)
    {
        libMesh::Node* node = *nd;
        
        // Get the global DOFs belonging to this node
        dof_map.dof_indices(node, dof_indices);
        
        // get the solution for this node
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        RealVectorX fpt = RealVectorX::Zero(ndofs);
        vec.setZero(ndofs);
        for (unsigned int i=0; i<dof_indices.size(); i++)
        {
            sol(i) = (*localized_solution)(dof_indices[i]);
            fpt(i) = (*localized_Fp)(dof_indices[i]);
        }
        
        output.set_node_data(node);
        output.output_derivative_for_node(sol, fpt, vec);
        output.clear_node();
        
        DenseRealVector v;
        MAST::copy(v, vec);
        dof_map.constrain_element_vector(v, dof_indices);
        dq_dX.add_vector(v, dof_indices);
        dof_indices.clear();
    }
        
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    dq_dX.close();
    
    output.clear_assembly();
}


void
MAST::AssemblyBase::
calculate_output_direct_sensitivity(const libMesh::NumericVector<Real>& X,
                                    bool if_localize_sol,
                                    const libMesh::NumericVector<Real>* dXdp,
                                    bool if_localize_sol_sens,
                                    const MAST::FunctionBase& p,
                                    MAST::OutputAssemblyElemOperations& output) {

    libmesh_assert(_discipline);
    libmesh_assert(_system);

    output.zero_for_sensitivity();

    MAST::NonlinearSystem& nonlin_sys = _system->system();
    output.set_assembly(*this);
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX
    sol,
    dsol;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = _system->system().get_dof_map();
    
    const libMesh::NumericVector<Real>
    *sol_vec  = nullptr,
    *dsol_vec = nullptr;
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sens;
    
    if (if_localize_sol) {
        localized_solution.reset(build_localized_vector(nonlin_sys, X).release());
        sol_vec = localized_solution.get();
    }
    else
        sol_vec = &X;
    
    if (dXdp) {
        if (if_localize_sol_sens) {
            localized_solution_sens.reset(build_localized_vector(nonlin_sys, *dXdp).release());
            dsol_vec = localized_solution_sens.get();
        }
        else
            dsol_vec = dXdp;
    }

    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X, false);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    nonlin_sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    nonlin_sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        if (diagonal_elem_subdomain_id.count(elem->subdomain_id()))
            continue;

        // no sensitivity computation assembly is neeed in these cases
        if (_param_dependence &&
            // if object is specified and elem does not depend on it
            !_param_dependence->if_elem_depends_on_parameter(*elem, p) &&
            // and if no sol_sens is given
            (!dXdp ||
             // or if it can be ignored for elem
             (dXdp && _param_dependence->override_flag)))
            continue;
            
        MAST::GeomElem geom_elem;
        output.set_elem_data(elem->dim(), *elem, geom_elem);
        geom_elem.init(*elem, *_system);

        if (!output.if_evaluate_for_element(geom_elem)) continue;
        
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
        
        //        if (_sol_function)
        //            physics_elem->attach_active_solution_function(*_sol_function);
        

        output.init(geom_elem);
        output.set_elem_solution(sol);
        output.set_elem_solution_sensitivity(dsol);
        output.evaluate_sensitivity(p);
        output.clear_elem();
        
        //        physics_elem->detach_active_solution_function();
    }
    
    /**
     * Now we need to iterate over node for node based operations such as 
     * the compliance due to point loads, a weighted sum of displacements,
     * etc.
     */
    
    // If point loads exist, assumble a vector of them
    std::unique_ptr<libMesh::NumericVector<Real>> Fp = 
        nonlin_sys.solution->zero_clone();
    calculate_point_load_vector(*Fp);
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_Fp;
    localized_Fp.reset(build_localized_vector(nonlin_sys, *Fp).release());
    
    // If the point loads exist, assemble of vector of their derivatives
    std::unique_ptr<libMesh::NumericVector<Real>> dFp = 
        nonlin_sys.solution->zero_clone();
    calculate_point_load_derivative_vector(p, *dFp);
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_dFp;
    localized_dFp.reset(build_localized_vector(nonlin_sys, *dFp).release());
        
    /**
     *  Now that we have point loads, we can now pass this information along
     *  with the displacements to the output object.
     */
    libMesh::MeshBase::node_iterator nd = 
        nonlin_sys.get_mesh().local_nodes_begin();
        
    const libMesh::MeshBase::node_iterator end_nd = 
        nonlin_sys.get_mesh().local_nodes_end();
    
    for ( ; nd != end_nd; ++nd)
    {
        libMesh::Node* node = *nd;
        
        // Get the global DOFs belonging to this node
        dof_map.dof_indices(node, dof_indices);
        
        // get the solution for this node
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        RealVectorX fpt = RealVectorX::Zero(ndofs);
        RealVectorX dfpt = RealVectorX::Zero(ndofs);
        for (unsigned int i=0; i<dof_indices.size(); i++)
        {
            sol(i) =  (*localized_solution)(dof_indices[i]);
            fpt(i) =  (*localized_Fp)(dof_indices[i]);
            dfpt(i) = (*localized_dFp)(dof_indices[i]);
        }
        
        output.set_node_data(node);
        output.evaluate_sensitivity_for_node(p, sol, fpt, dfpt);
        output.clear_node();
        
        dof_indices.clear();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
    
    output.clear_assembly();
}




Real
MAST::AssemblyBase::
calculate_output_adjoint_sensitivity(const libMesh::NumericVector<Real>& X,
                                     bool if_localize_sol,
                                     const libMesh::NumericVector<Real>& adj_sol,
                                     const MAST::FunctionBase& p,
                                     MAST::AssemblyElemOperations& elem_ops,
                                     MAST::OutputAssemblyElemOperations& output,
                                     const bool include_partial_sens) {

    libmesh_assert(_discipline);
    libmesh_assert(_system);
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();

    libMesh::NumericVector<Real>
    &dres_dp = nonlin_sys.add_sensitivity_rhs();
    
    this->set_elem_operation_object(elem_ops);
    this->sensitivity_assemble(X, if_localize_sol, p, dres_dp);
    this->clear_elem_operation_object();

    Real
    dq_dp = adj_sol.dot(dres_dp);

    if (include_partial_sens) {

        // calculate the partial sensitivity of the output, which is done
        // with zero solution vector
        this->calculate_output_direct_sensitivity(X, if_localize_sol, nullptr, false, p, output);

        dq_dp += output.output_sensitivity_total(p);
    }

    return dq_dp;
}


void
MAST::AssemblyBase::calculate_output_adjoint_sensitivity_multiple_parameters_no_direct
(const libMesh::NumericVector<Real>&           X,
 bool                                          if_localize_sol,
 const libMesh::NumericVector<Real>&           adj_sol,
 const std::vector<const MAST::FunctionBase*>& p_vec,
 MAST::AssemblyElemOperations&                 elem_ops,
 MAST::OutputAssemblyElemOperations&           output,
 std::vector<Real>&                            sens) {
    
    libmesh_assert(_discipline);
    libmesh_assert(_system);
    libmesh_assert_equal_to(sens.size(), p_vec.size());
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();

    // zero the sensitivity data first
    std::fill(sens.begin(), sens.end(), 0.);

     // add vectors before computing sensitivity
    for (unsigned int i=0; i<p_vec.size(); i++) nonlin_sys.add_sensitivity_rhs(i);

    // first compute all the residual vectors without closing them. later we will close them
    // and then compute the sensitivity
    for (unsigned int i=0; i<p_vec.size(); i++) {
        
        const MAST::FunctionBase& p = *p_vec[i];
        
        libMesh::NumericVector<Real>
        &dres_dp = nonlin_sys.get_sensitivity_rhs(i);
        
        this->set_elem_operation_object(elem_ops);
        this->sensitivity_assemble(X, if_localize_sol, p, dres_dp, false);
        this->clear_elem_operation_object();
    }
    
    for (unsigned int i=0; i<p_vec.size(); i++) {
        
        const MAST::FunctionBase& p = *p_vec[i];
        
        libMesh::NumericVector<Real>
        &dres_dp = nonlin_sys.add_sensitivity_rhs(i);
        dres_dp.close();
        sens[i] = adj_sol.dot(dres_dp);
    }
}


void MAST::AssemblyBase::calculate_point_load_vector(
    libMesh::NumericVector<Real>& Fp)
{
    Fp.zero();
    if (_discipline->point_loads().size())
    {
        const MAST::PointLoadSetType& loads = _discipline->point_loads();
        
        RealVectorX vec = RealVectorX::Zero(_system->n_vars());
        
        MAST::PointLoadSetType::const_iterator it  = loads.begin();
        MAST::PointLoadSetType::const_iterator end = loads.end();
        
        std::vector<libMesh::dof_id_type> dof_indices;
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        
        const libMesh::dof_id_type
        first_dof  = dof_map.first_dof(_system->system().comm().rank()),
        end_dof	   = dof_map.end_dof(_system->system().comm().rank());

        for ( ; it != end; it++) 
        {
            // get the point load function
            const MAST::FieldFunction<RealVectorX> &func = 
                (*it)->get<MAST::FieldFunction<RealVectorX>>("load");
            
            // get the nodes on which this object defines the load
            const std::set<const libMesh::Node*> nodes = (*it)->get_nodes();
            
            std::set<const libMesh::Node*>::const_iterator
            n_it    = nodes.begin(),
            n_end   = nodes.end();
            
            for (; n_it != n_end; n_it++) 
            {
                
                // load at the node
                vec.setZero();
                func(**n_it, _system->system().time, vec);
                
                dof_map.dof_indices(*n_it, dof_indices);

                libmesh_assert_equal_to(dof_indices.size(), vec.rows());

                // zero the components of the vector if they do not
                // belong to this processor
                for (unsigned int i=0; i<dof_indices.size(); i++)
                {
                    if (dof_indices[i] <   first_dof  || dof_indices[i] >=  end_dof)
                        vec(i) = 0.;
                }

                DenseRealVector v;
                MAST::copy(v, vec);
                
                dof_map.constrain_element_vector(v, dof_indices);
                Fp.add_vector(v, dof_indices);
                dof_indices.clear();
            }
        }
    }
    Fp.close();
}


void MAST::AssemblyBase::calculate_point_load_derivative_vector(
     const MAST::FunctionBase& p, libMesh::NumericVector<Real>& dpFp_dpparam)
{
    dpFp_dpparam.zero();
    if (_discipline->point_loads().size())
    {
        const MAST::PointLoadSetType& loads = _discipline->point_loads();
        
        RealVectorX vec = RealVectorX::Zero(_system->n_vars());
        
        MAST::PointLoadSetType::const_iterator it  = loads.begin();
        MAST::PointLoadSetType::const_iterator end = loads.end();
        
        std::vector<libMesh::dof_id_type> dof_indices;
        const libMesh::DofMap& dof_map = _system->system().get_dof_map();
        
        const libMesh::dof_id_type
        first_dof  = dof_map.first_dof(_system->system().comm().rank()),
        end_dof	   = dof_map.end_dof(_system->system().comm().rank());

        for ( ; it != end; it++) 
        {
            // get the point load function
            const MAST::FieldFunction<RealVectorX> &func = 
                (*it)->get<MAST::FieldFunction<RealVectorX>>("load");
            
            // get the nodes on which this object defines the load
            const std::set<const libMesh::Node*> nodes = (*it)->get_nodes();
            
            std::set<const libMesh::Node*>::const_iterator
            n_it    = nodes.begin(),
            n_end   = nodes.end();
            
            for (; n_it != n_end; n_it++) 
            {
                
                // load at the node
                vec.setZero();
                func.derivative(p, **n_it, _system->system().time, vec);
                
                dof_map.dof_indices(*n_it, dof_indices);

                libmesh_assert_equal_to(dof_indices.size(), vec.rows());

                // zero the components of the vector if they do not
                // belong to this processor
                for (unsigned int i=0; i<dof_indices.size(); i++)
                {
                    if (dof_indices[i] <   first_dof  || dof_indices[i] >=  end_dof)
                        vec(i) = 0.;
                }

                DenseRealVector v;
                MAST::copy(v, vec);
                
                dof_map.constrain_element_vector(v, dof_indices);
                dpFp_dpparam.add_vector(v, dof_indices);
                dof_indices.clear();
            }
        }
    }
    dpFp_dpparam.close();
}
