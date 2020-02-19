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
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"
#include "base/nonlinear_system.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"


MAST::GeomElem::GeomElem():
_sys_init        (nullptr),
_use_local_elem  (false),
_ref_elem        (nullptr),
_local_elem      (nullptr) {
    
}


MAST::GeomElem::~GeomElem() {
    
    if (_local_elem) {
        delete _local_elem;
        
        for (unsigned int i=0; i<_local_nodes.size(); i++)
            delete _local_nodes[i];
    }
}



const libMesh::Elem&
MAST::GeomElem::get_reference_elem() const {

    libmesh_assert(_ref_elem);
    
    return *_ref_elem;
}


const libMesh::Elem&
MAST::GeomElem::get_reference_local_elem() const {
    
    libmesh_assert(_local_elem);
    
    return *_local_elem;
}


const libMesh::Elem&
MAST::GeomElem::get_quadrature_elem() const {
    
    libmesh_assert(_ref_elem);
    
    return *_ref_elem;
}


const libMesh::Elem&
MAST::GeomElem::get_quadrature_local_elem() const {
    
    libmesh_assert(_local_elem);
    
    return *_local_elem;
}



unsigned int
MAST::GeomElem::dim() const {
    
    libmesh_assert(_ref_elem);
    
    return _ref_elem->dim();
}



unsigned int
MAST::GeomElem::n_sides_quadrature_elem() const {

    libmesh_assert(_ref_elem);

    return _ref_elem->n_sides();
}
        

libMesh::FEType
MAST::GeomElem::get_fe_type(unsigned int i) const {

    libmesh_assert(_ref_elem);
    
    return _sys_init->fetype(i);
}
        

void
MAST::GeomElem::set_local_y_vector(const RealVectorX& y_vec) {
    
    libmesh_assert_equal_to(y_vec.size(), 3);
    
    _local_y = y_vec;
}


void
MAST::GeomElem::set_bending(bool onoff) {
    _bending = onoff;;
}


void
MAST::GeomElem::init(const libMesh::Elem& elem,
                     const MAST::SystemInitialization& sys_init) {
    
    libmesh_assert(!_ref_elem);
    
    _ref_elem = &elem;
    _sys_init = &sys_init;
    
    // initialize the local element if needed. 
    _init_local_elem();
}


std::unique_ptr<MAST::FEBase>
MAST::GeomElem::init_fe(bool init_grads,
                        bool init_second_order_derivative,
                        int extra_quadrature_order) const {
    
    libmesh_assert(_ref_elem);
    
    std::unique_ptr<MAST::FEBase> fe(new MAST::FEBase(*_sys_init));
    fe->set_extra_quadrature_order(extra_quadrature_order);
    fe->set_evaluate_second_order_derivatives(init_second_order_derivative);
    
    fe->init(*this, init_grads);
    
    return fe;
}
        

std::unique_ptr<MAST::FEBase>
MAST::GeomElem::init_side_fe(unsigned int s,
                             bool init_grads,
                             bool init_second_order_derivative,
                             int extra_quadrature_order) const {
 
    std::unique_ptr<MAST::FEBase> fe(new MAST::FEBase(*_sys_init));
    fe->set_extra_quadrature_order(extra_quadrature_order);
    fe->set_evaluate_second_order_derivatives(init_second_order_derivative);

    fe->init_for_side(*this, s, init_grads);
    
    return fe;
}
        
        

void
MAST::GeomElem::external_side_loads_for_quadrature_elem
(std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc,
 std::map<unsigned int, std::vector<MAST::BoundaryConditionBase*>>& loads) const {
    
    // this assumes that the quadrature element is the same as the local element
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    std::pair<maptype::const_iterator, maptype::const_iterator> range;

    loads.clear();
    
    const libMesh::BoundaryInfo&
    binfo = *_sys_init->system().get_mesh().boundary_info;
    
    std::vector<libMesh::boundary_id_type> bids;
    
    for (unsigned int i=0; i<_ref_elem->n_sides(); i++) {
        
        bids.clear();
        binfo.boundary_ids(_ref_elem, i, bids);
        
        for (unsigned int j=0; j<bids.size(); j++) {
        
            range = bc.equal_range(bids[j]);

            maptype::const_iterator it = range.first;
            for ( ; it != range.second; it++) loads[i].push_back(it->second);
        }
    }
}
        


void
MAST::GeomElem::get_boundary_ids_on_quadrature_elem_side
(unsigned int s, std::vector<libMesh::boundary_id_type>& bc_ids) const {
    
    // this assumes that the quadrature element is the same as the local element
    bc_ids.clear();
    
    _sys_init->system().get_mesh().boundary_info->boundary_ids(_ref_elem, s, bc_ids);
}


const RealVectorX&
MAST::GeomElem::domain_surface_normal() const {

    libmesh_assert(_ref_elem);
    
    return _domain_surface_normal;
}


const RealMatrixX&
MAST::GeomElem::T_matrix() const {
 
    libmesh_assert(_local_elem);
    
    return _T_mat;
}
        

void
MAST::GeomElem::
transform_point_to_global_coordinate(const libMesh::Point& local_pt,
                                     libMesh::Point& global_pt) const {
    
    libmesh_assert(_local_elem);
    
    global_pt.zero();
    
    // now calculate the global coordinates with respect to the origin
    for (unsigned int j=0; j<3; j++)
        for (unsigned int k=0; k<3; k++)
            global_pt(j) += _T_mat(j,k)*local_pt(k);
    
    // shift to the global coordinate
    global_pt += (*_ref_elem->node_ptr(0));
}



void
MAST::GeomElem::
transform_vector_to_global_coordinate(const libMesh::Point& local_vec,
                                      libMesh::Point& global_vec) const {

    libmesh_assert(_local_elem);
    
    global_vec.zero();
    
    // now calculate the global coordinates with respect to the origin
    for (unsigned int j=0; j<3; j++)
        for (unsigned int k=0; k<3; k++)
            global_vec(j) += _T_mat(j,k)*local_vec(k);
}



void
MAST::GeomElem::
transform_vector_to_local_coordinate(const libMesh::Point& global_vec,
                                     libMesh::Point& local_vec) const {
    
    libmesh_assert(_local_elem);
    
    local_vec.zero();
    
    // now calculate the global coordinates with respect to the origin
    for (unsigned int j=0; j<3; j++)
        for (unsigned int k=0; k<3; k++)
            local_vec(j) += _T_mat(k,j) * global_vec(k);
}



void
MAST::GeomElem::
transform_vector_to_global_coordinate(const RealVectorX& local_vec,
                                      RealVectorX& global_vec) const {

    libmesh_assert(_local_elem);
    
    global_vec = _T_mat * local_vec;
}


void
MAST::GeomElem::
transform_vector_to_local_coordinate(const RealVectorX& global_vec,
                                     RealVectorX& local_vec) const {
    
    libmesh_assert(_local_elem);
    
    local_vec = _T_mat.transpose() * global_vec;
}



bool
MAST::GeomElem::use_local_elem() const {

    libmesh_assert(_ref_elem);
    return _use_local_elem;
}



void
MAST::GeomElem::_init_local_elem() {

    libmesh_assert(_ref_elem);
    libmesh_assert(!_local_elem);
    
    switch (_ref_elem->dim()) {
            
        case 1: {

            // if the y and z cooedinates of the nodes are the same then
            // the element is along the x-axis and should not need any
            // local element.
            const Real
            y  = _ref_elem->point(0)(1),
            z  = _ref_elem->point(0)(2);
            
            for (unsigned int i=1; i<_ref_elem->n_nodes(); i++) {
                if (std::fabs(y-_ref_elem->point(i)(1)) != 0. ||
                    std::fabs(z-_ref_elem->point(i)(2)) != 0.)
                    _use_local_elem = true;
            }

            if (_use_local_elem)
                _init_local_elem_1d();
        }
            break;
            
        case 2: {

            // if the z cooedinate of the nodes are the same then
            // the element is in the xy-plane and should not need any
            // local element.
            const Real
            z  = _ref_elem->point(0)(2);
            
            for (unsigned int i=1; i<_ref_elem->n_nodes(); i++) {
                if (std::fabs(z-_ref_elem->point(i)(2)) != 0.)
                    _use_local_elem = true;
            }
            
            if (_use_local_elem)
                _init_local_elem_2d();
        }
            break;
            
        case 3:
            _use_local_elem = false;
            break;
            
        default:
            libmesh_error(); // should not get here.
    }
    
}


void
MAST::GeomElem::_init_local_elem_1d() {
    
    libmesh_assert(_ref_elem);
    libmesh_assert(_use_local_elem);
    libmesh_assert(!_local_elem);
    libmesh_assert(_local_y.size());
    libmesh_assert_equal_to(_ref_elem->dim(), 1);
    
    if (_local_y.size()==0) // if the orientation vector has not been defined
    {
        if (!_bending) // if bending is not used in this element
        {   // then create a random orientation vector that is not collinear 
            // with the element's local x-axis. Added for github issue #40
            
            // Get element x-axis vector
            libMesh::Point v1;
            v1 = *_ref_elem->node_ptr(1); v1 -= *_ref_elem->node_ptr(0);

            // Perturb it by an arbitrary value and assign it to the _local_y vector.
            _local_y = RealVectorX::Zero(3);
            _local_y(0) = v1(0)+1.000;
            _local_y(1) = v1(1)-0.625;
            _local_y(2) = v1(2)+0.275;
        }
        else
        {
            libmesh_error_msg("ERROR: y_vector (for orientation of 1D element) not defined; In " << __PRETTY_FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__);
        }
    }
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    libMesh::Point v1, v2, v3, p;
    v1 = *_ref_elem->node_ptr(1); v1 -= *_ref_elem->node_ptr(0); v1 /= v1.norm(); // local x
    for (unsigned int i=0; i<3; i++) v2(i) = _local_y(i);       // vector in local x-y plane
    v3 = v1.cross(v2);                       // local z
    libmesh_assert_greater(v3.norm(), 0.);   // 0. implies x == y
    v3 /= v3.norm();
    v2 = v3.cross(v1); v2 /= v2.norm();      // local y
    
    _T_mat  = RealMatrixX::Zero(3,3);
    
    _local_elem = libMesh::Elem::build(_ref_elem->type()).release();
    _local_nodes.resize(_ref_elem->n_nodes());
    for (unsigned int i=0; i<_ref_elem->n_nodes(); i++) {
        _local_nodes[i] = new libMesh::Node;
        _local_nodes[i]->set_id() = _ref_elem->node_ptr(i)->id();
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
        p  = *_ref_elem->node_ptr(i);
        p -= *_ref_elem->node_ptr(0); // local wrt origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                (*_local_nodes[i])(j) += _T_mat(k,j)*p(k);
    }
}


void
MAST::GeomElem::_init_local_elem_2d() {
    
    libmesh_assert(_ref_elem);
    libmesh_assert(_use_local_elem);
    libmesh_assert(!_local_elem);
    libmesh_assert_equal_to(_ref_elem->dim(), 2);
    
    // first node is the origin of the new cs
    // calculate the coordinate system for the plane of the element
    libMesh::Point v1, v2, v3, p;
    v1 = *_ref_elem->node_ptr(1); v1 -= *_ref_elem->node_ptr(0); v1 /= v1.norm(); // local x
    v2 = *_ref_elem->node_ptr(2); v2 -= *_ref_elem->node_ptr(0); v2 /= v2.norm();
    v3 = v1.cross(v2); v3 /= v3.norm();      // local z
    v2 = v3.cross(v1); v2 /= v2.norm();      // local y
    
    // copy the surface normal
    _domain_surface_normal = RealVectorX::Zero(3);
    for (unsigned int i=0; i<3; i++)
        _domain_surface_normal(i) = v3(i);
    
    _T_mat = RealMatrixX::Zero(3,3);
    
    _local_elem = libMesh::Elem::build(_ref_elem->type()).release();
    _local_nodes.resize(_ref_elem->n_nodes());
    for (unsigned int i=0; i<_ref_elem->n_nodes(); i++) {
        _local_nodes[i] = new libMesh::Node;
        _local_nodes[i]->set_id() = _ref_elem->node_ptr(i)->id();
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
        p  = *_ref_elem->node_ptr(i);
        p -= *_ref_elem->node_ptr(0); // local wrt origin
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++)
                (*_local_nodes[i])(j) += _T_mat(k,j)*p(k);
    }

}
