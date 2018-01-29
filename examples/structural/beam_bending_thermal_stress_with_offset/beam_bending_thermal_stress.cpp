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
#include "examples/structural/beam_bending_thermal_stress_with_offset/beam_bending_thermal_stress.h"


MAST::Examples::BeamBendingThermalStress::BeamBendingThermalStress():
MAST::Examples::StructuralExample1D() {
    
}


MAST::Examples::BeamBendingThermalStress::~BeamBendingThermalStress() {
    
}



void
MAST::Examples::BeamBendingThermalStress::_init_section_property() {
    
    this->_init_section_property_with_offset();
}



void
MAST::Examples::BeamBendingThermalStress::_init_loads() {
    
    _init_temperature_load();
}


