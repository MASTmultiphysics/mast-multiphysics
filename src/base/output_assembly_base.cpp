/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "base/output_assembly_base.h"
#include "base/system_initialization.h"




MAST::OutputAssemblyBase::OutputAssemblyBase():
MAST::AssemblyBase(),
libMesh::System::QOI(),
libMesh::System::QOIDerivative(),
libMesh::System::QOIParameterSensitivity() {
    
}



MAST::OutputAssemblyBase::~OutputAssemblyBase() {
    
}



void
MAST::OutputAssemblyBase::
attach_discipline_and_system(MAST::PhysicsDisciplineBase &discipline,
                             MAST::SystemInitialization &system) {
    
    libmesh_assert_msg(!_discipline && !_system,
                       "Error: Assembly should be cleared before attaching System.");
    
    _discipline = &discipline;
    _system     = &system;
    
    // now attach this to the system
    libMesh::System& sys = system.system();
    sys.attach_QOI_object(*this);
    sys.attach_QOI_derivative_object(*this);
    sys.attach_QOI_parameter_sensitivity_object(*this);
}



void
MAST::OutputAssemblyBase::reattach_to_system() {
    
    libmesh_assert(_system);

    // now attach this to the system
    libMesh::System& sys = _system->system();
    sys.attach_QOI_object(*this);
    sys.attach_QOI_derivative_object(*this);
    sys.attach_QOI_parameter_sensitivity_object(*this);
}



void
MAST::OutputAssemblyBase::
clear_discipline_and_system() {
    
    if (_system && _discipline) {
        libMesh::System& sys = _system->system();
        
        sys.reset_QOI();
        sys.reset_QOI_derivative();
        sys.reset_QOI_parameter_sensitivity();
    }
    
    _discipline = NULL;
    _system     = NULL;
}





