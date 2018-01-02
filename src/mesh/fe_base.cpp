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
#include "mesh/fe_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"


MAST::FEBase::FEBase(const MAST::SystemInitialization& sys):
_sys        (sys),
_extra_quadrature_order(0),
_initialized(false),
_fe         (nullptr),
_qrule      (nullptr) {
    
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
MAST::FEBase::init(const libMesh::Elem& elem,
                   const std::vector<libMesh::Point>* pts) {
    
    libmesh_assert(!_initialized);
    
    const unsigned int
    nv      = _sys.n_vars();
    libMesh::FEType
    fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));
    
    // Create an adequate quadrature rule
    _fe = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    _fe->get_dphi();
    _fe->get_dphidxi();
    _fe->get_dphideta();
    _fe->get_dphidzeta();
    
    if (pts == nullptr) {
        _qrule = fe_type.default_quadrature_rule
        (elem.dim(),
         _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature
        _fe->attach_quadrature_rule(_qrule);
        _fe->reinit(&elem);
    }
    else
        _fe->reinit(&elem, pts);

    _initialized = true;
}


void
MAST::FEBase::init_for_side(const libMesh::Elem& elem,
                   unsigned int s,
                   bool if_calculate_dphi) {

    libmesh_assert(!_initialized);

    unsigned int nv = _sys.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _sys.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _sys.fetype(i));
    
    // Create an adequate quadrature rule
    _fe     = libMesh::FEBase::build(elem.dim(), fe_type).release();
    _qrule  = fe_type.default_quadrature_rule
    (elem.dim()-1,
     _sys.system().extra_quadrature_order+_extra_quadrature_order).release();  // system extra quadrature
    _fe->attach_quadrature_rule(_qrule);
    _fe->get_phi();
    _fe->get_xyz();
    _fe->get_JxW();
    if (if_calculate_dphi)
        _fe->get_dphi();
    
    _fe->reinit(&elem, s);

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
MAST::FEBase::get_normals() const {
    
    libmesh_assert(_initialized);
    return _fe->get_normals();
}


const libMesh::QBase&
MAST::FEBase::get_qrule() const {
    
    libmesh_assert(_initialized);
    return *_qrule;
}

