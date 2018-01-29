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
#include "examples/structural/plate_thermally_stressed_modal_analysis/plate_thermally_stressed_modal_analysis.h"


MAST::Examples::PlateModalAnalysis::PlateModalAnalysis():
MAST::Examples::StructuralExample2D() {
    
}





MAST::Examples::PlateModalAnalysis::~PlateModalAnalysis() {
    
}



void
MAST::Examples::PlateModalAnalysis::_init_loads() {
    
    this->_init_temperature_load();
    this->_init_pressure_load(false, 0);
}




