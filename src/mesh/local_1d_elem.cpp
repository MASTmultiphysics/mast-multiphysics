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
#include "mesh/local_1d_elem.h"


MAST::Local1DElem::Local1DElem(const libMesh::Elem& elem,
                               const libMesh::Point& y):
MAST::LocalElemBase(elem),
_local_y(y) {
    
    _create_local_elem();
    libmesh_assert_greater(_local_y.size(), 0.);
}


MAST::Local1DElem::~Local1DElem() {
    
    // the local element may not have been created
    // for cases where the original element lies in the xy-plane
    if (_local_elem) {
        delete _local_elem;
        for (unsigned int i=0; i<_local_nodes.size(); i++)
            delete _local_nodes[i];
    }
}


void 
MAST::Local1DElem::_create_local_elem() {
    
    libmesh_assert(_elem.dim() == 1);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    libMesh::Point v1, v2, v3, p;
    v1 = *_elem.get_node(1); v1 -= *_elem.get_node(0); v1 /= v1.size(); // local x
    v2 = _local_y;                           // vector in local x-y plane
    v3 = v1.cross(v2);                       // local z
    libmesh_assert_greater(v3.size(), 0.);   // 0. implies x == y
    v3 /= v3.size();
    v2 = v3.cross(v1); v2 /= v2.size();      // local y
    
    _T_mat  = RealMatrixX::Zero(3,3);
    
    _local_elem = libMesh::Elem::build(_elem.type()).release();
    _local_nodes.resize(_elem.n_nodes());
    for (unsigned int i=0; i<_elem.n_nodes(); i++) {
        _local_nodes[i] = new libMesh::Node;
        _local_nodes[i]->set_id() = _elem.get_node(i)->id();
        _local_elem->set_node(i) = _local_nodes[i];
    }
    
    // now the transformation matrix from old to new cs
    //        an_i vn_i = a_i v_i
    //        an_j = a_i v_i.vn_j  = a_i t_ij = T^t a_i
    //        t_ij = v_i.vn_j
    
    for (unsigned int i=0; i<3; i++) {
        _T_mat(i,0) = v1(i);
        _T_mat(i,1) = v2(i);
        _T_mat(i,2) = v3(i);
    }
    
    // now calculate the new coordinates with respect to the origin
    for (unsigned int i=0; i<_local_nodes.size(); i++) {
        p = *_elem.get_node(i);
        p -= *_elem.get_node(0); // local wrt origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                (*_local_nodes[i])(j) += _T_mat(k,j)*p(k);
    }

}


void
MAST::Local1DElem::
domain_surface_normal_in_global_coordinates(const libMesh::Point& p,
                                            RealVector3& n_global) const {
    
    for (unsigned int i=0; i<3; i++)
        n_global(i) = _local_y(i);
    
    n_global /= _local_y.size();
}




