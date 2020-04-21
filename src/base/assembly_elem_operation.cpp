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
#include "base/assembly_elem_operation.h"
#include "base/elem_base.h"
#include "base/assembly_base.h"
#include "base/physics_discipline_base.h"
#include "base/boundary_condition_base.h"
#include "base/field_function_base.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"






MAST::AssemblyElemOperations::AssemblyElemOperations():
_system           (nullptr),
_discipline       (nullptr),
_assembly         (nullptr),
_physics_elem     (nullptr),
_skip_comm_sum    (false) {
    
}


MAST::AssemblyElemOperations::~AssemblyElemOperations() {
    
    this->clear_elem();
}


void
MAST::AssemblyElemOperations::
set_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                          MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
}



void
MAST::AssemblyElemOperations::
clear_discipline_and_system() {
    
    _discipline    = nullptr;
    _system        = nullptr;
}


void
MAST::AssemblyElemOperations::set_assembly(MAST::AssemblyBase &assembly) {

    libmesh_assert(!_assembly);
    _assembly = &assembly;
}


void
MAST::AssemblyElemOperations::clear_assembly() {
    
    _assembly = nullptr;
}


MAST::SystemInitialization&
MAST::AssemblyElemOperations::get_system_initialization() {

    return *_system;
}


MAST::PhysicsDisciplineBase&
MAST::AssemblyElemOperations::get_discipline() {
    
    return *_discipline;
}


MAST::AssemblyBase&
MAST::AssemblyElemOperations::get_assembly() {
    
    libmesh_assert(_assembly);
    return *_assembly;
}



void
MAST::AssemblyElemOperations::set_elem_solution(const RealVectorX& sol) {

    libmesh_assert(_physics_elem);

    _physics_elem->set_solution(sol);
}


void
MAST::AssemblyElemOperations::set_elem_solution_sensitivity(const RealVectorX& sol) {
    
    libmesh_assert(_physics_elem);

    _physics_elem->set_solution(sol, true);
}



void
MAST::AssemblyElemOperations::set_elem_velocity(const RealVectorX &vel) {
    
    libmesh_assert(_physics_elem);

    _physics_elem->set_velocity(vel);
}


void
MAST::AssemblyElemOperations::set_elem_velocity_sensitivity(const RealVectorX &vel) {
    
    libmesh_assert(_physics_elem);

    _physics_elem->set_velocity(vel, true);
}



void
MAST::AssemblyElemOperations::set_elem_acceleration(const RealVectorX &accel) {
    
    libmesh_assert(_physics_elem);

    _physics_elem->set_acceleration(accel);
}


void
MAST::AssemblyElemOperations::set_elem_acceleration_sensitivity(const RealVectorX &accel) {
    
    libmesh_assert(_physics_elem);
    
    _physics_elem->set_acceleration(accel, true);
}


void
MAST::AssemblyElemOperations::set_elem_perturbed_solution(const RealVectorX &sol) {
    
    libmesh_assert(_physics_elem);
    
    _physics_elem->set_perturbed_solution(sol);
}


void
MAST::AssemblyElemOperations::set_elem_perturbed_velocity(const RealVectorX &vel) {
    
    libmesh_assert(_physics_elem);
    
    _physics_elem->set_perturbed_velocity(vel);
}


void
MAST::AssemblyElemOperations::set_elem_perturbed_acceleration(const RealVectorX &accel) {
    
    libmesh_assert(_physics_elem);
    
    _physics_elem->set_perturbed_acceleration(accel);
}


void
MAST::AssemblyElemOperations::clear_elem() {
    
    if (_physics_elem)
        delete _physics_elem;
        
    _physics_elem = nullptr;
}


std::pair<const MAST::FieldFunction<RealVectorX>*, unsigned int>
MAST::AssemblyElemOperations::
get_elem_boundary_velocity_data() {
    
    libmesh_assert(_physics_elem);

    std::pair<const MAST::FieldFunction<RealVectorX>*, unsigned int>
    val = std::make_pair(nullptr, 0);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>> loads;
    _physics_elem->elem().external_side_loads_for_quadrature_elem(_discipline->side_loads(),
                                                                  loads);
    
    std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>::const_iterator
    it   = loads.begin(),
    end  = loads.end();
    
    for ( ; it != end; it++) {
        
        std::vector<MAST::BoundaryConditionBase*>::const_iterator
        bc_it  = it->second.begin(),
        bc_end = it->second.end();
        
        for ( ; bc_it != bc_end; bc_it++) {
            
            // apply all the types of loading
            switch ((*bc_it)->type()) {
                    
                case MAST::BOUNDARY_VELOCITY: {
                    
                    val.first  = &((*bc_it)->get<MAST::FieldFunction<RealVectorX>>("phi_vel"));
                    val.second = it->first;
                }
                    break;
                    
                    
                default:
                    // nothing to be done here
                    break;
            }
        }
    }
    
    return val;
}

