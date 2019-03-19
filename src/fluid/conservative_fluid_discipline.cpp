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
#include "fluid/conservative_fluid_discipline.h"




MAST::ConservativeFluidDiscipline::
ConservativeFluidDiscipline(libMesh::EquationSystems& eq_sys):
MAST::PhysicsDisciplineBase(eq_sys),
_flight_cond(nullptr)
{ }



MAST::ConservativeFluidDiscipline::~ConservativeFluidDiscipline()
{ }


void
MAST::ConservativeFluidDiscipline::set_flight_condition(MAST::FlightCondition &flt) {
    
    // make sure that the pointer is not being overwritten
    libmesh_assert(!_flight_cond);
    
    _flight_cond = &flt;
}


const MAST::FlightCondition&
MAST::ConservativeFluidDiscipline::flight_condition() const {

    // make sure that the pointer has been set
    libmesh_assert(_flight_cond);
    
    return *_flight_cond;
    
}

    
void
MAST::ConservativeFluidDiscipline::clear_flight_condition() {
    
    // clears the pointer
    _flight_cond = nullptr;
 }
