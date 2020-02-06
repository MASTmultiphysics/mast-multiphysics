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
#include "level_set/level_set_intersected_elem.h"
#include "level_set/level_set_intersection.h"
#include "level_set/sub_cell_fe.h"


MAST::LevelSetIntersectedElem::LevelSetIntersectedElem():
MAST::GeomElem  (),
_sub_elem       (nullptr),
_local_sub_elem (nullptr),
_intersection   (nullptr) {
    
}



MAST::LevelSetIntersectedElem::~LevelSetIntersectedElem() {
    
    if (_local_sub_elem) {
        delete _local_sub_elem;
        
        for (unsigned int i=0; i<_local_subelem_nodes.size(); i++)
            delete _local_subelem_nodes[i];
    }
}


const libMesh::Elem&
MAST::LevelSetIntersectedElem::get_quadrature_elem() const {
    
    libmesh_assert(_ref_elem);
    
    return *_sub_elem;
}


const libMesh::Elem&
MAST::LevelSetIntersectedElem::get_quadrature_local_elem() const {
    
    libmesh_assert(_local_elem);
    
    return *_local_sub_elem;
}


bool
MAST::LevelSetIntersectedElem::if_elem_has_level_set_boundary() const {
    
    libmesh_assert(_intersection);

    return  _intersection->if_elem_has_boundary();
}



bool
MAST::LevelSetIntersectedElem::if_subelem_has_side_on_level_set_boundary() const {
    
    libmesh_assert(_intersection);
 
    return _intersection->has_side_on_interface(*_sub_elem);
}



int
MAST::LevelSetIntersectedElem::get_subelem_side_on_level_set_boundary() const {
    
    libmesh_assert(_intersection);
    
    return _intersection->get_side_on_interface(*_sub_elem);
}



void
MAST::LevelSetIntersectedElem::init(const libMesh::Elem& elem,
                                    const MAST::SystemInitialization& sys_init,
                                    MAST::LevelSetIntersection& intersection) {
    
    libmesh_assert(!_intersection);
    
    
    _intersection  = &intersection;
    _sub_elem      = &elem;
    _ref_elem      = &intersection.elem();
    _sys_init      = &sys_init;
    
    // initialize the local element if needed. (not implemented yet)
    //_init_local_elem();
}


std::unique_ptr<MAST::FEBase>
MAST::LevelSetIntersectedElem::init_fe(bool init_grads,
                                       bool init_second_order_derivative,
                                       int extra_quadrature_order) const {
    
    libmesh_assert(_intersection);
    std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_sys_init, *_intersection));
    fe->set_extra_quadrature_order(extra_quadrature_order);
    fe->set_evaluate_second_order_derivatives(init_second_order_derivative);
    
    fe->init(*this, init_grads);
    return fe;
}


std::unique_ptr<MAST::FEBase>
MAST::LevelSetIntersectedElem::init_side_fe(unsigned int s,
                                            bool init_grads,
                                            bool init_second_order_derivative,
                                            int extra_quadrature_order) const {

    libmesh_assert(_intersection);
    std::unique_ptr<MAST::FEBase> fe(new MAST::SubCellFE(*_sys_init, *_intersection));
    fe->set_extra_quadrature_order(extra_quadrature_order);
    fe->set_evaluate_second_order_derivatives(init_second_order_derivative);

    fe->init_for_side(*this, s, init_grads);

    return fe;
}

