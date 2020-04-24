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
#include "level_set/intersected_elem_homogenization_function.h"
#include "level_set/level_set_intersection.h"
#include "level_set/level_set_intersected_elem.h"


// libMesh includes
#include "libmesh/elem.h"

MAST::IntersectedElemHomogenizedDensityFunction::
IntersectedElemHomogenizedDensityFunction(const std::string& nm):
MAST::HomogenizedDensityFunctionBase(nm),
_intersection   (nullptr) {
    
}



MAST::IntersectedElemHomogenizedDensityFunction::
~IntersectedElemHomogenizedDensityFunction() {
    
}


void
MAST::IntersectedElemHomogenizedDensityFunction::
initialize_element_volume_fractions() {
    
    libmesh_assert(_analysis_mesh);
    libmesh_assert(_intersection);
    libmesh_assert(_elem_volume_fraction.empty());
        
    libMesh::MeshBase::const_element_iterator
    e_it  =  _analysis_mesh->local_elements_begin(),
    e_end  =  _analysis_mesh->local_elements_begin();
    
    for ( ; e_it != e_end; e_it++) {
        
        _intersection->init(*_level_set, **e_it, 0.,
                            _analysis_mesh->max_elem_id(),
                            _analysis_mesh->max_node_id());

        _elem_volume_fraction[*e_it] = _intersection->get_positive_phi_volume_fraction();
        _intersection->clear();
    }
}


void
MAST::IntersectedElemHomogenizedDensityFunction::
initialize_element_volume_fraction_sensitivity(const MAST::FunctionBase& f) {

    libmesh_assert(_analysis_mesh);
    libmesh_assert(_intersection);
    libmesh_assert(_elem_volume_fraction_sensitivity.empty());

    libMesh::MeshBase::const_element_iterator
    e_it  =  _analysis_mesh->local_elements_begin(),
    e_end  =  _analysis_mesh->local_elements_begin();
    
    for ( ; e_it != e_end; e_it++) {
        
        const libMesh::Elem* elem = *e_it;
        
        _intersection->init(*_level_set, *elem, 0.,
                            _analysis_mesh->max_elem_id(),
                            _analysis_mesh->max_node_id());
        
        // iterate over all elements on the positive side and
        // add compute the sensitivity of the volume
        const std::vector<const libMesh::Elem *> &
        elems_hi = _intersection->get_sub_elems_positive_phi();
        
        std::vector<const libMesh::Elem*>::const_iterator
        hi_sub_elem_it  = elems_hi.begin(),
        hi_sub_elem_end = elems_hi.end();
        
//        for (; hi_sub_elem_it != hi_sub_elem_end; hi_sub_elem_it++ ) {
//
//            const libMesh::Elem* sub_elem = *hi_sub_elem_it;
//
//            MAST::LevelSetIntersectedElem geom_elem;
//            ops.set_elem_data(elem->dim(), *elem, geom_elem);
//            geom_elem.init(*sub_elem, *_system, *_intersection);
//
//            ops.init(geom_elem);
//
//            ops.elem_calculations(J!=nullptr?true:false, sub_elem_vec, sub_elem_mat);
//
//            ops.clear_elem();
//        }

        _elem_volume_fraction_sensitivity[&f][elem] =
        _intersection->get_positive_phi_volume_fraction();
        
        _intersection->clear();
    }
}

