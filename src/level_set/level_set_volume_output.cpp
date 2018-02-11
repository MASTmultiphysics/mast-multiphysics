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
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_elem_base.h"
#include "level_set/level_set_intersection.h"
#include "base/assembly_base.h"


MAST::LevelSetVolume::LevelSetVolume(MAST::LevelSetIntersection& intersection):
MAST::OutputAssemblyElemOperations(),
_intersection  (intersection),
_vol           (0.),
_dvol_dp       (0.) {
    
}



MAST::LevelSetVolume::~LevelSetVolume() {

    
}



void
MAST::LevelSetVolume::init(const libMesh::Elem& elem) {
    
    libmesh_assert(!_physics_elem);
    libmesh_assert(_system);
    libmesh_assert(_assembly);
    
    _physics_elem = new MAST::LevelSetElementBase(*_system, *_assembly, elem);
}





void
MAST::LevelSetVolume::zero() {

    _vol     = 0.;
    _dvol_dp = 0.;
}



Real
MAST::LevelSetVolume::output_for_elem() {
    
    libmesh_assert(_physics_elem);
    
    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        MAST::LevelSetElementBase&
        e = dynamic_cast<MAST::LevelSetElementBase&>(*_physics_elem);
        
        return e.volume();
    }
    else
        return 0.;
}



Real
MAST::LevelSetVolume::output_total() {

    return _vol;
}



Real
MAST::LevelSetVolume::output_sensitivity_for_elem(const MAST::FunctionBase& p) {

    libmesh_assert(_physics_elem);

    if (this->if_evaluate_for_element(_physics_elem->elem()) &&
        _intersection.if_elem_has_boundary()) {

        MAST::LevelSetElementBase&
        e = dynamic_cast<MAST::LevelSetElementBase&>(*_physics_elem);
        
        return e.volume_boundary_velocity_on_side(_intersection.get_side_on_interface(_physics_elem->elem()));
    }
    else
        return 0.;
    
}



Real
MAST::LevelSetVolume::output_sensitivity_total(const MAST::FunctionBase& p) {

    return _dvol_dp;
}


void
MAST::LevelSetVolume::evaluate() {

    libmesh_assert(_physics_elem);

    if (this->if_evaluate_for_element(_physics_elem->elem())) {
        
        MAST::LevelSetElementBase&
        e = dynamic_cast<MAST::LevelSetElementBase&>(*_physics_elem);
        
        _vol += e.volume();
    }
}


void
MAST::LevelSetVolume::evaluate_sensitivity(const MAST::FunctionBase& f) {

    libmesh_assert(_physics_elem);

    // sensitivity only exists at the boundary. So, we proceed with calculation
    // only if this element has an intersection in the interior, or with the
    // node or a side. 
    if (this->if_evaluate_for_element(_physics_elem->elem()) &&
        _intersection.if_elem_has_boundary()) {
        
        MAST::LevelSetElementBase&
        e = dynamic_cast<MAST::LevelSetElementBase&>(*_physics_elem);
        
        _dvol_dp += e.volume_boundary_velocity_on_side(_intersection.get_side_on_interface(_physics_elem->elem()));
    }
}


bool
MAST::LevelSetVolume::if_use_local_elem() const {

    return false;
}



void
MAST::LevelSetVolume::set_local_fe_data(MAST::LocalElemFE& fe,
                                        const libMesh::Elem& e) const {

    libmesh_assert(false); // should not get called
}

