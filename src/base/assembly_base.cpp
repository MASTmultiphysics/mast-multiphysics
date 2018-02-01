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
#include "base/assembly_base.h"
#include "base/system_initialization.h"
#include "base/mesh_field_function.h"
#include "base/elem_base.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "mesh/local_elem_fe.h"
#include "base/assembly_elem_operation.h"
#include "base/output_assembly_elem_operations.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"


MAST::AssemblyBase::AssemblyBase():
_elem_ops         (nullptr),
_discipline       (nullptr),
_system           (nullptr),
_sol_function     (nullptr),
_solver_monitor   (nullptr) {
    
}




MAST::AssemblyBase::~AssemblyBase() {
    
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
attach_discipline_and_system(MAST::AssemblyElemOperations& elem_ops,
                             MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system && !_elem_ops,
                       "Error: Assembly should be cleared before attaching System.");
    
    _elem_ops   = &elem_ops;
    _discipline = &discipline;
    _system     = &system;
    
    _elem_ops->set_assembly(*this);
}



void
MAST::AssemblyBase::
clear_discipline_and_system( ) {
    
    _elem_ops->clear_assembly();
    
    _elem_ops      = nullptr;
    _discipline    = nullptr;
    _system        = nullptr;
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



std::unique_ptr<MAST::FEBase>
MAST::AssemblyBase::build_fe(const libMesh::Elem& elem) {

    std::unique_ptr<MAST::FEBase> fe;

    if (_elem_ops->if_use_local_elem() &&
        elem.dim() < 3) {
        
        MAST::LocalElemFE*
        local_fe = new MAST::LocalElemFE(*_system);
        _elem_ops->set_local_fe_data(*local_fe, elem);
        
        fe.reset(local_fe);
    }
    else {
        
        fe.reset(new MAST::FEBase(*_system));
    }
    
    return fe;
}




void
MAST::AssemblyBase::calculate_output(const libMesh::NumericVector<Real>& X,
                                     MAST::OutputAssemblyElemOperations& output) {
    
    
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
        
        dof_map.dof_indices (elem, dof_indices);
        
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        
        //if (_sol_function)
        //    physics_elem->attach_active_solution_function(*_sol_function);
        
        
        // perform the element level calculations
        output.evaluate();
        
        //physics_elem->detach_active_solution_function();
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}




void
MAST::AssemblyBase::
calculate_output_derivative(const libMesh::NumericVector<Real>& X,
                            MAST::OutputAssemblyElemOperations& output,
                            libMesh::NumericVector<Real>& dq_dX) {
    
    
    MAST::NonlinearSystem& nonlin_sys = _system->system();
    
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
        
        dof_map.dof_indices (elem, dof_indices);
        
        output.init(*elem);
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        output.set_elem_solution(sol);
        
        
//        if (_sol_function)
//            physics_elem->attach_active_solution_function(*_sol_function);
        
        //_check_element_numerical_jacobian(*physics_elem, sol);
        
        // perform the element level calculations
        output.evaluate_derivative();
        
//        physics_elem->detach_active_solution_function();
        
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}





