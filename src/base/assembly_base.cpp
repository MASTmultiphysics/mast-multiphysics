/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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
#include "mesh/fe_base.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/parameter_vector.h"


MAST::AssemblyBase::AssemblyBase():
_discipline(nullptr),
_system(nullptr),
_sol_function(nullptr) {
    
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




std::unique_ptr<libMesh::NumericVector<Real> >
MAST::AssemblyBase::_build_localized_vector(const libMesh::System& sys,
                                            const libMesh::NumericVector<Real>& global) {
    
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
MAST::AssemblyBase::calculate_outputs(const libMesh::NumericVector<Real>& X) {
    
    
    MAST::NonlinearSystem& sys = _system->system();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    std::unique_ptr<MAST::ElementBase> physics_elem;
    
    std::unique_ptr<libMesh::NumericVector<Real> > localized_solution;
    localized_solution.reset(_build_localized_vector(sys, X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init(X);
    
    
    libMesh::MeshBase::const_element_iterator       el     =
    sys.get_mesh().active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
    sys.get_mesh().active_local_elements_end();
    
    
    for ( ; el != end_el; ++el) {
        
        const libMesh::Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        physics_elem.reset(_build_elem(*elem).release());
        
        // get the solution
        unsigned int ndofs = (unsigned int)dof_indices.size();
        sol.setZero(ndofs);
        vec.setZero(ndofs);
        mat.setZero(ndofs, ndofs);
        
        for (unsigned int i=0; i<dof_indices.size(); i++)
            sol(i) = (*localized_solution)(dof_indices[i]);
        
        physics_elem->set_solution(sol);
        
        if (_sol_function)
            physics_elem->attach_active_solution_function(*_sol_function);
        
        // perform the element level calculations
        _elem_outputs(*physics_elem,
                      _discipline->volume_output(),
                      _discipline->side_output());
        
        physics_elem->detach_active_solution_function();
    }
    
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}





void
MAST::AssemblyBase::
calculate_output_sensitivity(libMesh::ParameterVector &params,
                             const bool if_total_sensitivity,
                             const libMesh::NumericVector<Real> &X) {
    
    
    MAST::NonlinearSystem& sys = _system->system();
    
    // iterate over each element, initialize it and get the relevant
    // analysis quantities
    RealVectorX vec, sol, sol_sens;
    RealMatrixX mat;
    
    std::vector<libMesh::dof_id_type> dof_indices;
    const libMesh::DofMap& dof_map = sys.get_dof_map();
    std::unique_ptr<MAST::ElementBase> physics_elem;
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    localized_solution,
    localized_solution_sensitivity;
    
    localized_solution.reset(_build_localized_vector(sys, X).release());
    
    
    // if a solution function is attached, initialize it
    if (_sol_function)
        _sol_function->init( X);
    
    // iterate over the parameters
    for ( unsigned int i=0; i<params.size(); i++) {
        
        if (if_total_sensitivity)
            localized_solution_sensitivity.reset
            (_build_localized_vector(sys,
                                     sys.get_sensitivity_solution(i)).release());
    
        libMesh::MeshBase::const_element_iterator       el     =
        sys.get_mesh().active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el =
        sys.get_mesh().active_local_elements_end();
        
        
        for ( ; el != end_el; ++el) {
            
            const libMesh::Elem* elem = *el;
            
            dof_map.dof_indices (elem, dof_indices);
            
            physics_elem.reset(_build_elem(*elem).release());
            
            // get the solution
            unsigned int ndofs = (unsigned int)dof_indices.size();
            sol.setZero(ndofs);
            sol_sens.setZero(ndofs);
            vec.setZero(ndofs);
            mat.setZero(ndofs, ndofs);

            // tell the element about the sensitivity paramete
            physics_elem->sensitivity_param =
            _discipline->get_parameter(&(params[i].get()));
            
            // get the solution
            for (unsigned int i=0; i<dof_indices.size(); i++)
                sol(i) = (*localized_solution)(dof_indices[i]);

            // tell the element about the solution
            physics_elem->set_solution(sol);

            // get the solution sensitivity
            if (if_total_sensitivity)
                for (unsigned int i=0; i<dof_indices.size(); i++)
                    sol_sens(i) = (*localized_solution_sensitivity)(dof_indices[i]);
        
            // tell the solution about the sensitivity
            physics_elem->set_solution(sol_sens, true);

            if (_sol_function)
                physics_elem->attach_active_solution_function(*_sol_function);
            
            // perform the element level calculations
            _elem_output_sensitivity(*physics_elem,
                                     _discipline->volume_output(),
                                     _discipline->side_output());
            
            physics_elem->detach_active_solution_function();
        }
        
    }
    
    // if a solution function is attached, clear it
    if (_sol_function)
        _sol_function->clear();
}



std::unique_ptr<MAST::FEBase>
MAST::AssemblyBase::build_fe(const libMesh::Elem& elem) {
    
    std::unique_ptr<MAST::FEBase>
    fe(new MAST::FEBase(*_system));
    
    return fe;
}



void
MAST::AssemblyBase::
_elem_outputs(MAST::ElementBase &elem,
              std::multimap<libMesh::subdomain_id_type,MAST::OutputFunctionBase *> &vol_output,
              std::multimap<libMesh::boundary_id_type,MAST::OutputFunctionBase *> &side_output) {
    
    
    // ask the element to provide the outputs
    elem.volume_output_quantity(false, false, vol_output);
    elem.side_output_quantity(false, false,  side_output);
}




void
MAST::AssemblyBase::
_elem_output_sensitivity(MAST::ElementBase &elem,
                         std::multimap<libMesh::subdomain_id_type, MAST::OutputFunctionBase *> &vol_output,
                         std::multimap<libMesh::boundary_id_type,MAST::OutputFunctionBase *> &side_output) {
    
    
    // ask the element to provide the outputs
    elem.volume_output_quantity(false,  // false for adjoints
                                true,   // true for sensitivity
                                vol_output);
    elem.side_output_quantity(false,  // false for adjoints
                              true,   // true for sensitivity
                              side_output);
}



