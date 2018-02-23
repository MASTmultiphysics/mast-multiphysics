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
#include "mesh/local_elem_fe.h"
#include "mesh/local_1d_elem.h"
#include "mesh/local_2d_elem.h"
#include "mesh/local_3d_elem.h"


MAST::LocalElemFE::LocalElemFE(const MAST::SystemInitialization& sys):
MAST::FEBase  (sys),
_local_elem   (nullptr) {
    
}



MAST::LocalElemFE::~LocalElemFE() {
    
    if (_local_elem) delete _local_elem;
}


void
MAST::LocalElemFE::init(const libMesh::Elem& elem,
                        const std::vector<libMesh::Point>* pts) {

    libmesh_assert(!_initialized);
    
    _create_local_element(elem);
    
    // now that this element has been initialized, use it to initialize
    // the FE object.
    MAST::FEBase::init(_local_elem->local_elem(), pts);
    
    // now initialize the global xyz locations
    const std::vector<libMesh::Point>
    local_xyz    = _fe->get_xyz();
    
    unsigned int
    n = (unsigned int) local_xyz.size();
    _global_xyz.resize(n);
    
    for (unsigned int i=0; i<n; i++)
        _local_elem->global_coordinates_location(local_xyz[i], _global_xyz[i]);
}


void
MAST::LocalElemFE::init_for_side(const libMesh::Elem& elem,
                                 unsigned int s,
                                 bool if_calculate_dphi) {
    
    libmesh_assert(!_initialized);
    
    _create_local_element(elem);
    
    // now that this element has been initialized, use it to initialize
    // the FE object.
    MAST::FEBase::init_for_side(_local_elem->local_elem(), s, if_calculate_dphi);
    
    // now initialize the global xyz locations and normals
    const std::vector<libMesh::Point>
    &local_xyz     = _fe->get_xyz(),
    &local_normals = _fe->get_normals();
    
    unsigned int
    n = (unsigned int) local_xyz.size();
    _global_xyz.resize(n);
    _global_normals.resize(n);
    
    for (unsigned int i=0; i<n; i++) {
        
        _local_elem->global_coordinates_location(local_xyz[i], _global_xyz[i]);
        _local_elem->global_coordinates_normal(local_normals[i], _global_normals[i]);
    }
}




void
MAST::LocalElemFE::set_1d_y_vector(const libMesh::Point& y) {
    
    libmesh_assert_greater(y.norm(), 0.);
    _y_vector_1d_elem = y;
}


const libMesh::Point&
MAST::LocalElemFE::get_1d_y_vector() const {
    
    libmesh_assert_greater(_y_vector_1d_elem.norm(), 0.);
    return _y_vector_1d_elem;
}


void
MAST::LocalElemFE::_create_local_element(const libMesh::Elem& elem) {
    
    libmesh_assert(!_initialized);
    
    switch (elem.dim()) {
        case 1: {
            
            libmesh_assert_greater(_y_vector_1d_elem.norm(), 0.);
            _local_elem = new MAST::Local1DElem(elem, _y_vector_1d_elem);
        }
            break;
            
        case 2:
            _local_elem = new MAST::Local2DElem(elem);
            break;
            
        case 3:
            _local_elem = new MAST::Local3DElem(elem);
            break;
            
        default:
            // should not get here.
            libmesh_error();
            break;
    }
}


