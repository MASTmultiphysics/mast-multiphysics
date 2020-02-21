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
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"


MAST::FEBase::FEBase(const MAST::SystemInitialization& sys):
_sys                           (sys),
_extra_quadrature_order        (0),
_init_second_order_derivatives (false),
_initialized                   (false),
_elem                          (nullptr),
_fe                            (nullptr),
_qrule                         (nullptr) {
    
}


MAST::FEBase::~FEBase() {
    
    if (_fe)
        delete _fe;
    
    if (_qrule)
        delete _qrule;
}


void
MAST::FEBase::set_extra_quadrature_order(int n) {
    
    _extra_quadrature_order = n;
}


void
MAST::FEBase::set_evaluate_second_order_derivatives(bool f) {
    
    _init_second_order_derivatives = f;
}


void
MAST::FEBase::init(const MAST::GeomElem& elem,
                   bool init_grads,
                   const std::vector<libMesh::Point>* pts) {
    
    libmesh_assert(!_initialized);
    
    
    _elem   = &elem;
    const unsigned int
    nv      = _sys.n_vars();
    libMesh::FEType
    
    fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    // TODO: What effect does this assumption have for future MAST developments? What about beams with mixed Lagrange/Hermite interpolation?
    
    if (_sys.prefix().compare("warping"))
    {
        for (unsigned int i=1; i != nv; ++i)
        {
            libmesh_assert(fe_type == _sys.fetype(i));
        }
    }
    
    const libMesh::Elem*
    q_elem = nullptr;
    
    if (elem.use_local_elem())
        q_elem = &elem.get_quadrature_local_elem();
    else
        q_elem = &elem.get_quadrature_elem();
    
    // Create an adequate quadrature rule
    _fe = libMesh::FEBase::build(q_elem->dim(), fe_type).release();
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
    
    if (pts == nullptr) {
        _qrule = fe_type.default_quadrature_rule
        (q_elem->dim(),
         _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature
        _fe->attach_quadrature_rule(_qrule);
        _fe->reinit(q_elem);
    }
    else {
        _fe->reinit(q_elem, pts);
        _qpoints = *pts;
    }

    // now initialize the global xyz locations if the element uses a locally
    // defined element.
    if (elem.use_local_elem()) {
        
        const std::vector<libMesh::Point>
        local_xyz    = _fe->get_xyz();
        
        unsigned int
        n = (unsigned int) local_xyz.size();
        _global_xyz.resize(n);
        
        for (unsigned int i=0; i<n; i++)
            elem.transform_point_to_global_coordinate(local_xyz[i], _global_xyz[i]);
    }
    
    _initialized = true;
}





void
MAST::FEBase::init_for_side(const MAST::GeomElem& elem,
                            unsigned int s,
                            bool if_calculate_dphi) {

    libmesh_assert(!_initialized);

    _elem = &elem;
    
    const unsigned int
    nv    = _sys.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));
    
    const libMesh::Elem*
    q_elem = nullptr;
    
    if (elem.use_local_elem())
        q_elem = &elem.get_quadrature_local_elem();
    else
        q_elem = &elem.get_quadrature_elem();

    // Create an adequate quadrature rule
    _fe     = libMesh::FEBase::build(q_elem->dim(), fe_type).release();
    _qrule  = fe_type.default_quadrature_rule
    (q_elem->dim()-1,
     _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature
    _fe->attach_quadrature_rule(_qrule);
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    _fe->get_normals();
    if (if_calculate_dphi) {

        _fe->get_dphi();
        if (_init_second_order_derivatives) _fe->get_d2phi();
    }
    
    _fe->reinit(q_elem, s);


    // now initialize the global xyz locations and normals
    _local_normals = _fe->get_normals();
    if (elem.use_local_elem()) {

        const std::vector<libMesh::Point>
        &local_xyz     = _fe->get_xyz();
        
        unsigned int
        n = (unsigned int) local_xyz.size();
        _global_xyz.resize(n);
        _global_normals.resize(n);
        
        for (unsigned int i=0; i<n; i++) {
            
            elem.transform_point_to_global_coordinate(local_xyz[i], _global_xyz[i]);
            elem.transform_vector_to_global_coordinate(_local_normals[i], _global_normals[i]);
        }
    }
    else {
        _global_normals = _local_normals;
        _global_xyz     = _fe->get_xyz();
    }

    
    _initialized = true;
}




libMesh::FEType
MAST::FEBase::get_fe_type() const {

    libmesh_assert(_initialized);
    return _fe->get_fe_type();
}


const std::vector<Real>&
MAST::FEBase::get_JxW() const {
    
    libmesh_assert(_initialized);
    return _fe->get_JxW();
}


const std::vector<libMesh::Point>&
MAST::FEBase::get_xyz() const {
    
    libmesh_assert(_initialized);
    if (_elem->use_local_elem())
        return _global_xyz;
    else
        return _fe->get_xyz();
}


unsigned int
MAST::FEBase::n_shape_functions() const {
    
    libmesh_assert(_initialized);
    return _fe->n_shape_functions();
}


const std::vector<std::vector<Real> >&
MAST::FEBase::get_phi() const {
    
    libmesh_assert(_initialized);
    return _fe->get_phi();
}


const std::vector<std::vector<libMesh::RealVectorValue> >&
MAST::FEBase::get_dphi() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dphi();
}


const std::vector<std::vector<libMesh::RealTensorValue>>&
MAST::FEBase::get_d2phi() const {
    
    libmesh_assert(_initialized);
    return _fe->get_d2phi();
}


const std::vector<Real>&
MAST::FEBase::get_dxidx() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxidx();
}


const std::vector<Real>&
MAST::FEBase::get_dxidy() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxidy();
}


const std::vector<Real>&
MAST::FEBase::get_dxidz() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxidz();
}


const std::vector<Real>&
MAST::FEBase::get_detadx() const {
    
    libmesh_assert(_initialized);
    return _fe->get_detadx();
}


const std::vector<Real>&
MAST::FEBase::get_detady() const {
    
    libmesh_assert(_initialized);
    return _fe->get_detady();
}


const std::vector<Real>&
MAST::FEBase::get_detadz() const {
    
    libmesh_assert(_initialized);
    return _fe->get_detadz();
}


const std::vector<Real>&
MAST::FEBase::get_dzetadx() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dzetadx();
}


const std::vector<Real>&
MAST::FEBase::get_dzetady() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dzetady();
}


const std::vector<Real>&
MAST::FEBase::get_dzetadz() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dzetadz();
}


const std::vector<libMesh::RealVectorValue>&
MAST::FEBase::get_dxyzdxi() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxyzdxi();
}


const std::vector<libMesh::RealVectorValue>&
MAST::FEBase::get_dxyzdeta() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxyzdeta();
}


const std::vector<libMesh::RealVectorValue>&
MAST::FEBase::get_dxyzdzeta() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dxyzdzeta();
}


const std::vector<std::vector<Real> >&
MAST::FEBase::get_dphidxi() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dphidxi();
}


const std::vector<std::vector<Real> >&
MAST::FEBase::get_dphideta() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dphideta();
}


const std::vector<std::vector<Real> >&
MAST::FEBase::get_dphidzeta() const {
    
    libmesh_assert(_initialized);
    return _fe->get_dphidzeta();
}


const std::vector<libMesh::Point>&
MAST::FEBase::get_normals_for_reference_coordinate() const {
    
    libmesh_assert(_initialized);
    return _global_normals;
}


const std::vector<libMesh::Point>&
MAST::FEBase::get_normals_for_local_coordinate() const {

    libmesh_assert(_initialized);
    return _local_normals;
}


const std::vector<libMesh::Point>&
MAST::FEBase::get_qpoints() const {
    
    libmesh_assert(_initialized);
    if (_qrule)  // qrule was used
        return _qrule->get_points();
    else         // points were specified
        return _qpoints;
}


const libMesh::QBase&
MAST::FEBase::get_qrule() const {
    
    libmesh_assert(_initialized);
    return *_qrule;
}

