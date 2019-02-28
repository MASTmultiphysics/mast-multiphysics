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
#include "mesh/local_elem_base.h"
#include "base/field_function_base.h"


MAST::LocalElemBase::~LocalElemBase() {
    
}



void
MAST::LocalElemBase::global_coordinates_location(const libMesh::Point& local,
                                                 libMesh::Point& global) const {
    if (!_local_elem) // no local elem is created for the case of 3D elem
        global = local;
    else {
        global = 0.;
        
        // now calculate the global coordinates with respect to the origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                global(j) += _T_mat(j,k)*local(k);
        
        // shift to the global coordinate
        if (_elem.parent())
            global += (*_elem.parent()->node_ptr(0));
        else
            global += (*_elem.node_ptr(0));
    }
}




void
MAST::LocalElemBase::global_coordinates_normal(const libMesh::Point& local,
                                               libMesh::Point& global) const {
    if (!_local_elem)
        for (unsigned int i=0; i<3; i++)
            global(i) = local(i);
    else {
        global.zero();
        
        // now calculate the global coordinates with respect to the origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                global(j) += _T_mat(j,k)*local(k);
    }
}



