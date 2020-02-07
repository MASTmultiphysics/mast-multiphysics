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
#include "level_set/sub_cell_fe.h"
#include "level_set/level_set_intersection.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "mesh/geom_elem.h"


// libMesh includes
#include "libmesh/elem.h"


MAST::SubCellFE::SubCellFE(const MAST::SystemInitialization& sys,
                           const MAST::LevelSetIntersection& intersection):
MAST::FEBase    (sys),
_intersection   (intersection) {
    
}



MAST::SubCellFE::~SubCellFE() {
    
}


void
MAST::SubCellFE::init(const MAST::GeomElem& elem,
                      bool init_grads,
                      const std::vector<libMesh::Point>* pts) {

    // if there was no intersection, then move back to the parent class's
    // method
    if (!_intersection.if_intersection_through_elem()) {
        MAST::FEBase::init(elem, init_grads, pts);
        return;
    }
    
    libmesh_assert(!_initialized);
    
    // this method does not allow quadrature points to be spcified.
    libmesh_assert(!pts);

    _elem    = &elem;
    
    // the quadrature element is the subcell on which the quadrature is to be
    // performed and the reference element is the element inside which the
    // subcell is defined
    const libMesh::Elem
    *ref_elem = nullptr,
    *q_elem   = nullptr;
    
    if (elem.use_local_elem()) {
        ref_elem = &elem.get_reference_local_elem();
        q_elem   = &elem.get_quadrature_local_elem();
    }
    else {
        ref_elem = &elem.get_reference_elem();
        q_elem   = &elem.get_quadrature_elem();
    }
    
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
    // the quadrature rule on the subcell that will be used to identify the
    // quadrature point weight in the
    std::unique_ptr<libMesh::QBase>
    subcell_qrule(fe_type.default_quadrature_rule
                  (q_elem->dim(),
                   _sys.system().extra_quadrature_order+_extra_quadrature_order).release());
    // subcell and local FE are always going to be initialized for linear
    // shape functions. This is because we create linear sub-elements
    // and since libmesh does not allow high-order FE on linear elements.
    // We will initialize the quadrature point locations on these sub elements
    // and then map it back to the high-order element
    libMesh::FEType
    linear_fetype(libMesh::FIRST, libMesh::LAGRANGE);
    std::unique_ptr<libMesh::FEBase>
    local_fe  (libMesh::FEBase::build(q_elem->dim(), linear_fetype).release()),
    subcell_fe(libMesh::FEBase::build(q_elem->dim(), linear_fetype).release());

    // quadrature points in the non-dimensional coordinate system of
    // reference element are obtained from the local_fe
    local_fe->get_xyz();

    // the JxW values are obtained from the quadrature element, and its
    // sum should be equal to the area of the quadrature element.
    // initialize the sub-cell FE to get the JxW coordinates
    subcell_fe->get_JxW();

    // create the FE object and tell what quantities are needed from it.
    _fe     = libMesh::FEBase::build(q_elem->dim(), fe_type).release();
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    if (init_grads) {
        _fe->get_dphi();
        _fe->get_dphidxi();
        _fe->get_dphideta();
        _fe->get_dphidzeta();
    }
    if (_init_second_order_derivatives) _fe->get_d2phi();
    
    
    //////////////////////////////////////////////////////////////
    // create the nondimensional coordinate nodes
    //////////////////////////////////////////////////////////////
    std::vector<libMesh::Node*>
    local_nodes(q_elem->n_nodes(), nullptr);
    
    std::unique_ptr<libMesh::Elem>
    local_coord_elem(libMesh::Elem::build(q_elem->type()).release());
    local_coord_elem->set_id(q_elem->id());
    
    // note that the quadrature element is the sub-element that was created
    // inside the reference element. We intend to use the dofs and shape
    // functions in the reference element.
    //
    // this local element is in the xi-eta coordinate system of the reference
    // element. Then, if we initialize the quadrature on this local element
    // then the physical location of these quadrature points will be the
    // locations in the xi-eta space of the reference element, which we then
    // use to initialize the FE on reference element.
    for (unsigned int i=0; i<q_elem->n_nodes(); i++) {
        const libMesh::Node
        *n = q_elem->node_ptr(i);
        const libMesh::Point
        &p = _intersection.get_nondimensional_coordinate_for_node(*n);
        local_nodes[i] = new libMesh::Node(p);
        // set the node id since libMesh does not allow invalid ids. This
        // should not influence the operations here
        local_nodes[i]->set_id(n->id());
        local_coord_elem->set_node(i) = local_nodes[i];
    }

    // in the nondimensional coordinate
    local_fe->attach_quadrature_rule(subcell_qrule.get());
    local_fe->reinit(local_coord_elem.get());
    
    // we use the xyz points of the local_elem , which should be the location
    // of the quadrature points in the reference element nondimensional c/s.
    // The shape functions and their derivatives are going to come on
    // the parent elmeent, since the computations use the solution vector
    // for local computations on that element.
    _fe->reinit(ref_elem, &local_fe->get_xyz());

    // We use the _subcell_fe to get the
    // JxW, since the area needs to come from the element on which
    // the integration is being performed.
    subcell_fe->attach_quadrature_rule(subcell_qrule.get());
    subcell_fe->reinit(q_elem);
    _subcell_JxW = subcell_fe->get_JxW();
    _qpoints     = local_fe->get_xyz();
    
    // transform the normal and
    if (elem.use_local_elem()) {
        
        // now initialize the global xyz locations and normals
        const std::vector<libMesh::Point>
        &local_xyz     = _fe->get_xyz();
        
        unsigned int
        n = (unsigned int) local_xyz.size();
        _global_xyz.resize(n);
        
        for (unsigned int i=0; i<n; i++)
            _elem->transform_point_to_global_coordinate(local_xyz[i], _global_xyz[i]);
    }
    else
        _global_xyz     = _fe->get_xyz();

    // now delete the node pointers
    for (unsigned int i=0; i<q_elem->n_nodes(); i++)
        delete local_nodes[i];
    
    _initialized = true;
}



void
MAST::SubCellFE::init_for_side(const MAST::GeomElem& elem,
                               unsigned int s,
                               bool if_calculate_dphi) {
    
    if (!_intersection.if_intersection_through_elem()) {
        
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

    _elem     = &elem;
    
    // the quadrature element is the subcell on which the quadrature is to be
    // performed and the reference element is the element inside which the
    // subcell is defined
    const libMesh::Elem
    *ref_elem = nullptr,
    *q_elem   = nullptr;
    
    if (elem.use_local_elem()) {
        ref_elem = &elem.get_reference_local_elem();
        q_elem   = &elem.get_quadrature_local_elem();
    }
    else {
        ref_elem = &elem.get_reference_elem();
        q_elem   = &elem.get_quadrature_elem();
    }

    const unsigned int
    nv      = _sys.n_vars();
    libMesh::FEType
    fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));

    // initialize the quadrature rule since we will need this later and it
    // it does not change
    std::unique_ptr<libMesh::QBase>
    subcell_qrule(fe_type.default_quadrature_rule
                  (q_elem->dim()-1,
                   _sys.system().extra_quadrature_order+_extra_quadrature_order).release());

    // subcell and local FE are always going to be initialized for linear
    // shape functions. This is because we create linear sub-elements
    // and since libmesh does not allow high-order FE on linear elements.
    // We will initialize the quadrature point locations on these sub elements
    // and then map it back to the high-order element
    libMesh::FEType
    linear_fetype(libMesh::FIRST, libMesh::LAGRANGE);

    std::unique_ptr<libMesh::FEBase>
    local_fe  (libMesh::FEBase::build(q_elem->dim(), linear_fetype).release()),
    subcell_fe(libMesh::FEBase::build(q_elem->dim(), linear_fetype).release());

    // quadrature points in the non-dimensional coordinate system of
    // reference element are obtained from the local_fe
    local_fe->get_xyz();
    
    // the JxW values are obtained from the quadrature element, and its
    // sum should be equal to the area of the quadrature element.
    // initialize the sub-cell FE to get the JxW coordinates
    subcell_fe->get_JxW();
    subcell_fe->get_normals();

    _fe     = libMesh::FEBase::build(q_elem->dim(), fe_type).release();
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    if (if_calculate_dphi)
        _fe->get_dphi();
    if (_init_second_order_derivatives) _fe->get_d2phi();
    

    //////////////////////////////////////////////////////////////
    // create the nondimensional coordinate nodes
    //////////////////////////////////////////////////////////////
    std::vector<libMesh::Node*>
    local_nodes(q_elem->n_nodes(), nullptr);
    
    std::unique_ptr<libMesh::Elem>
    local_coord_elem(libMesh::Elem::build(q_elem->type()).release());
    local_coord_elem->set_id(q_elem->id());

    for (unsigned int i=0; i<q_elem->n_nodes(); i++) {
        
        const libMesh::Node
        *n = q_elem->node_ptr(i);
        const libMesh::Point
        &p = _intersection.get_nondimensional_coordinate_for_node(*n);
        local_nodes[i] = new libMesh::Node(p);
        // set the node id since libMesh does not allow invalid ids. This
        // should not influence the operations here
        local_nodes[i]->set_id(n->id());
        local_coord_elem->set_node(i) = local_nodes[i];
    }
    
    // in the nondimensional coordinate, initialize for the side
    local_fe->attach_quadrature_rule(subcell_qrule.get());
    local_fe->reinit(local_coord_elem.get(), s);
    
    // initialize parent FE using location of these points
    // The shape functions and their derivatives are going to come on
    // the parent elmeent, since the computations use the solution vector
    // for local computations on that element.
    _fe->reinit(ref_elem, &local_fe->get_xyz());
    
    // We use the _subcell_fe to get the
    // JxW and surface normals, since the area and normals need to come
    // from the element on which the integration is being performed.
    subcell_fe->attach_quadrature_rule(subcell_qrule.get());
    try {
        subcell_fe->reinit(q_elem, s);
        _subcell_JxW = subcell_fe->get_JxW();
        _qpoints     = local_fe->get_xyz();
        // normals in the coordinate system attached to the reference element
        _local_normals = subcell_fe->get_normals();
    }
    catch (...) {
        unsigned int npts = subcell_qrule->n_points();
        _subcell_JxW.resize(npts, 0.);
        _qpoints.resize(npts);
        _local_normals.resize(npts);
    }
    
    // transform the normal and
    if (elem.use_local_elem()) {
        
        // now initialize the global xyz locations and normals
        const std::vector<libMesh::Point>
        &local_xyz     = _fe->get_xyz();
        
        unsigned int
        n = (unsigned int) local_xyz.size();
        _global_xyz.resize(n);
        _global_normals.resize(n);
        
        for (unsigned int i=0; i<n; i++) {
            
            _elem->transform_point_to_global_coordinate(local_xyz[i], _global_xyz[i]);
            _elem->transform_vector_to_global_coordinate(_local_normals[i], _global_normals[i]);
        }
    }
    else {
        _global_normals = _local_normals;
        _global_xyz     = _fe->get_xyz();
    }
    
    // now delete the node pointers
    for (unsigned int i=0; i<q_elem->n_nodes(); i++)
        delete local_nodes[i];

    
    _initialized = true;
}



const
std::vector<Real>&
MAST::SubCellFE::get_JxW() const {
    libmesh_assert(_initialized);
    
    // if there was no intersection, then move back to the parent class's
    // method
    if (_subcell_JxW.size())
        return _subcell_JxW;
    else
        return MAST::FEBase::get_JxW();
}


const std::vector<libMesh::Point>&
MAST::SubCellFE::get_qpoints() const {
    
    libmesh_assert(_initialized);

    if (_subcell_JxW.size())
        return _qpoints;
    else
        return MAST::FEBase::get_qpoints();
}


