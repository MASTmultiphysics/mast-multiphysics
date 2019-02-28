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
#include "level_set/sub_cell_fe.h"
#include "level_set/level_set_intersection.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "mesh/local_elem_base.h"


// libMesh includes
#include "libmesh/elem.h"


MAST::SubCellFE::SubCellFE(const MAST::SystemInitialization& sys,
                           const MAST::LevelSetIntersection& intersection):
MAST::FEBase    (sys),
_intersection   (intersection),
_subcell_fe     (nullptr),
_subcell_qrule  (nullptr) {
    
}



MAST::SubCellFE::~SubCellFE() {
    
    if (_subcell_fe)     delete _subcell_fe;
    if (_subcell_qrule)  delete _subcell_qrule;
}


void
MAST::SubCellFE::init(const libMesh::Elem& elem,
                      const std::vector<libMesh::Point>* pts) {

    // if there was no intersection, then move back to the parent class's
    // method
    if (!_intersection.if_intersection_through_elem() ||
        elem.parent() == nullptr) {
        MAST::FEBase::init(elem, pts);
        return;
    }
    
    libmesh_assert(!_initialized);

    const unsigned int
    nv      = _sys.n_vars();
    libMesh::FEType
    fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));

    // we first initialize the sub-cell in the nondimensional coordinate system
    // and then use the non-dimensional location of the quadrature point to
    // initialize the shape functions of the parent element. After finding the
    // locations of these points, we reinitialize the sub-cell in the physical
    // coordinates to get the JxW quantities

    // initialize the quadrature rule since we will need this later and it
    // it does not change
    if (pts == nullptr) {
        _subcell_qrule = fe_type.default_quadrature_rule
        (elem.dim(),
         _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature
    }
    else
        // Not handled yet. The specified points are in the
        // parent element and we need to initialize the sub-cell element
        // only if any of these points lie in the sub-cell.
        libmesh_assert(false);
    
    //////////////////////////////////////////////////////////////
    // create the nondimensional coordinate nodes
    //////////////////////////////////////////////////////////////
    std::vector<libMesh::Node*>
    _local_nodes(elem.n_nodes(), nullptr);
    
    std::unique_ptr<libMesh::Elem>
    local_coord_elem(libMesh::Elem::build(elem.type()).release());
    local_coord_elem->set_id(elem.id());
    
    for (unsigned int i=0; i<elem.n_nodes(); i++) {
        const libMesh::Node
        *n = elem.node_ptr(i);
        const libMesh::Point
        &p = _intersection.get_nondimensional_coordinate_for_node(*n);
        _local_nodes[i] = new libMesh::Node(p);
        // set the node id since libMesh does not allow invalid ids. This
        // should not influence the operations here
        _local_nodes[i]->set_id(n->id());
        local_coord_elem->set_node(i) = _local_nodes[i];
    }
    
    // initialize the sub-cell FE to get the JxW coordinates
    _subcell_fe = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _subcell_fe->get_JxW();

    _fe     = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    _fe->get_dphi();
    _fe->get_dphidxi();
    _fe->get_dphideta();
    _fe->get_dphidzeta();
    if (_init_second_order_derivatives) _fe->get_d2phi();

    // now, we use the xyz points of this element, which should be the location
    // of the quadrature points in the parent element, since elem was in the
    // nondimensional coordinate system
    const libMesh::Elem
    *parent = elem.parent();
    
    // make sure parent and elem have the same dimensions
    libmesh_assert_equal_to(elem.dim(), parent->dim());

    // initialize an FE object on this to get its nondimensional
    // coordinates
    std::unique_ptr<libMesh::FEBase>
    local_fe(libMesh::FEBase::build(elem.dim(), fe_type).release());
    local_fe->get_xyz();
    
    if (pts == nullptr) {
        // in the nondimensional coordinate
        local_fe->attach_quadrature_rule(_subcell_qrule);
        local_fe->reinit(local_coord_elem.get());
        
        // We use the _subcell_fe to get the
        // JxW, since the area needs to come from the element on which
        // the integration is being performed.
        _subcell_fe->attach_quadrature_rule(_subcell_qrule);
        _subcell_fe->reinit(&elem);
        
        // initialize parent FE using location of these points
        // The shape functions and their derivatives are going to come on
        // the parent elmeent, since the computations use the solution vector
        // for local computations on that element.
        _fe->reinit(parent, &local_fe->get_xyz());
        _qpoints = local_fe->get_xyz();
    }
    else
        libmesh_assert(false);

    
    // now delete the node pointers
    for (unsigned int i=0; i<elem.n_nodes(); i++)
        delete _local_nodes[i];
    
    _initialized = true;
}



void
MAST::SubCellFE::init_for_side(const libMesh::Elem& elem,
                               unsigned int s,
                               bool if_calculate_dphi) {
    
    if (!_intersection.if_intersection_through_elem() ||
        elem.parent() == nullptr) {
        
        // if there was no intersection, then move back to the parent class's
        // method.
        
        MAST::FEBase::init_for_side(elem, s, if_calculate_dphi);
        return;
    }
    
    // Otherwise, we follow the same procedure as the initialization
    // over the domain, where the quadrature points from the subcell
    // are mapped back to the original element and then the original
    // element is initialized.
    //
    // This has to be carefully done. The JxW are needed on the sub-elem,
    // however, the shape functions (and derivatives) are needed on the
    // parent element. Additionally, the surface normal are also
    // evaluated at the sub-elem.
    //
    

    libmesh_assert(!_initialized);

    const unsigned int
    nv      = _sys.n_vars();
    libMesh::FEType
    fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));

    // initialize the quadrature rule since we will need this later and it
    // it does not change
    _subcell_qrule = fe_type.default_quadrature_rule
    (elem.dim()-1,
     _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature

    //////////////////////////////////////////////////////////////
    // create the nondimensional coordinate nodes
    //////////////////////////////////////////////////////////////
    std::vector<libMesh::Node*>
    _local_nodes(elem.n_nodes(), nullptr);
    
    std::unique_ptr<libMesh::Elem>
    local_coord_elem(libMesh::Elem::build(elem.type()).release());
    local_coord_elem->set_id(elem.id());

    for (unsigned int i=0; i<elem.n_nodes(); i++) {
        
        const libMesh::Node
        *n = elem.node_ptr(i);
        const libMesh::Point
        &p = _intersection.get_nondimensional_coordinate_for_node(*n);
        _local_nodes[i] = new libMesh::Node(p);
        // set the node id since libMesh does not allow invalid ids. This
        // should not influence the operations here
        _local_nodes[i]->set_id(n->id());
        local_coord_elem->set_node(i) = _local_nodes[i];
    }
    
    // initialize the sub-cell FE to get the JxW coordinates
    _subcell_fe = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _subcell_fe->get_JxW();
    _subcell_fe->get_normals();

    _fe     = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    if (if_calculate_dphi)
        _fe->get_dphi();

    // now, we use the xyz points of this element, which should be the location
    // of the quadrature points in the parent element, since elem was in the
    // nondimensional coordinate system
    const libMesh::Elem
    *parent = elem.parent();
    
    // make sure parent and elem have the same dimensions
    libmesh_assert_equal_to(elem.dim(), parent->dim());
    
    // initialize an FE object on this to get its nondimensional
    // coordinates
    std::unique_ptr<libMesh::FEBase>
    local_fe(libMesh::FEBase::build(elem.dim(), fe_type).release());
    local_fe->get_xyz();
    
    // in the nondimensional coordinate, initialize for the side
    local_fe->attach_quadrature_rule(_subcell_qrule);
    local_fe->reinit(local_coord_elem.get(), s);
    
    // We use the _subcell_fe to get the
    // JxW and surface normals, since the area and normals need to come
    // from the element on which the integration is being performed.
    _subcell_fe->attach_quadrature_rule(_subcell_qrule);
    _subcell_fe->reinit(&elem, s);
    
    // initialize parent FE using location of these points
    // The shape functions and their derivatives are going to come on
    // the parent elmeent, since the computations use the solution vector
    // for local computations on that element.
    _fe->reinit(parent, &local_fe->get_xyz());
    _qpoints = local_fe->get_xyz();
    
    // now delete the node pointers
    for (unsigned int i=0; i<elem.n_nodes(); i++)
        delete _local_nodes[i];

    
    _initialized = true;
}


void
MAST::SubCellFE::init_for_side(const MAST::LocalElemBase& elem,
                            unsigned int s,
                            bool if_calculate_dphi) {
    
    libmesh_assert(!_initialized);
    
    _local_elem = &elem;
    
    // now that this element has been initialized, use it to initialize
    // the FE object.
    this->init_for_side(_local_elem->local_elem(), s, if_calculate_dphi);
    
    // now initialize the global xyz locations and normals
    const std::vector<libMesh::Point>
    &local_xyz     = _fe->get_xyz(),
    *local_normals = nullptr;
    
    if (_local_elem->local_elem().parent())
        local_normals  = &(_subcell_fe->get_normals());
    else
        local_normals =  &(_fe->get_normals());
    
    unsigned int
    n = (unsigned int) local_xyz.size();
    _global_xyz.resize(n);
    _global_normals.resize(n);
    
    for (unsigned int i=0; i<n; i++) {
        
        _local_elem->global_coordinates_location(local_xyz[i], _global_xyz[i]);
        _local_elem->global_coordinates_normal((*local_normals)[i], _global_normals[i]);
    }
}


const
std::vector<Real>&
MAST::SubCellFE::get_JxW() const {
    libmesh_assert(_initialized);
    
    // if there was no intersection, then move back to the parent class's
    // method
    if (_subcell_fe)
        return _subcell_fe->get_JxW();
    else
        return MAST::FEBase::get_JxW();
}


const
std::vector<libMesh::Point>&
MAST::SubCellFE::get_normals() const {
    libmesh_assert(_initialized);
    
    // if there was no intersection, then move back to the parent class's
    // method
    if (_subcell_fe)
        return _subcell_fe->get_normals();
    else
        return MAST::FEBase::get_normals();
}


const std::vector<libMesh::Point>&
MAST::SubCellFE::get_qpoints() const {
    
    libmesh_assert(_initialized);

    if (_subcell_fe)
        return _qpoints;
    else
        return MAST::FEBase::get_qpoints();
}


