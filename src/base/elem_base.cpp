/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2019  Manav Bhatia
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
#include "base/elem_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "mesh/fe_base.h"


MAST::ElementBase::ElementBase(MAST::SystemInitialization& sys,
                               MAST::AssemblyBase& assembly,
                               const MAST::GeomElem& elem):
_system                 (sys),
_assembly               (assembly),
_elem                   (elem),
_active_sol_function    (nullptr),
_time                   (_system.system().time) {
    
}


MAST::ElementBase::~ElementBase() {

}


MAST::NonlinearSystem&
MAST::ElementBase::system() {
    return _system.system();
}




const RealVectorX&
MAST::ElementBase::sol(bool if_sens) const {
    if (!if_sens)
        return _sol;
    else
        return _sol_sens;
}


void
MAST::ElementBase::set_solution(const RealVectorX &vec,
                                bool if_sens) {
    
    if (!if_sens)
        _sol = vec;
    else
        _sol_sens = vec;
}



void
MAST::ElementBase::set_perturbed_solution(const RealVectorX &vec,
                                          bool if_sens) {
    
    if (!if_sens)
        _delta_sol = vec;
    else
        _delta_sol_sens = vec;
    
}


void
MAST::ElementBase::set_complex_solution(const ComplexVectorX &vec,
                                        bool if_sens) {
    
    if (!if_sens)
        _complex_sol = vec;
    else
        _complex_sol_sens = vec;
    
}




void
MAST::ElementBase::set_velocity(const RealVectorX &vec,
                                bool if_sens) {
    
    if (!if_sens)
        _vel = vec;
    else
        _vel_sens = vec;
}



void
MAST::ElementBase::set_perturbed_velocity(const RealVectorX &vec,
                                          bool if_sens) {
    
    if (!if_sens)
        _delta_vel = vec;
    else
        _delta_vel_sens = vec;
}



void
MAST::ElementBase::set_acceleration(const RealVectorX &vec,
                                    bool if_sens) {
    
    if (!if_sens)
        _accel = vec;
    else
        _accel_sens = vec;
}



void
MAST::ElementBase::set_perturbed_acceleration(const RealVectorX &vec,
                                              bool if_sens) {
    
    if (!if_sens)
        _delta_accel      = vec;
    else
        _delta_accel_sens = vec;
}



void
MAST::ElementBase::attach_active_solution_function(MAST::FunctionBase &f) {
    
    // make sure that this has not already been set
    libmesh_assert(!_active_sol_function);
    
    _active_sol_function = &f;
}


void
MAST::ElementBase::detach_active_solution_function() {
    
    _active_sol_function = nullptr;
}



