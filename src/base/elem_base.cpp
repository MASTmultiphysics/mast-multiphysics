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
#include "base/elem_base.h"
#include "base/system_initialization.h"
#include "mesh/local_elem_base.h"


MAST::ElementBase::ElementBase(MAST::SystemInitialization& sys,
                               const libMesh::Elem& elem):
sensitivity_param(NULL),
_system(sys),
_elem(elem),
_active_sol_function(NULL),
_time(_system.system().time),
_fe(NULL),
_qrule(NULL) {
    
}


MAST::ElementBase::~ElementBase() {

    if (_fe)     delete _fe;
    if (_qrule)  delete _qrule;
}


libMesh::System&
MAST::ElementBase::system() {
    return _system.system();
}




const RealVectorX&
MAST::ElementBase::sol(bool if_sens) const {
    if (!if_sens)
        return _sol;
    else
        return _sol_sens;
}


void
MAST::ElementBase::set_solution(const RealVectorX &vec,
                                bool if_sens) {
    
    if (!if_sens)
        _sol = vec;
    else
        _sol_sens = vec;
}



void
MAST::ElementBase::set_complex_solution(const ComplexVectorX &vec,
                                        bool if_sens) {
    
    if (!if_sens)
        _complex_sol = vec;
    else
        _complex_sol_sens = vec;
    
}




void
MAST::ElementBase::set_velocity(const RealVectorX &vec,
                                bool if_sens) {
    
    if (!if_sens)
        _vel = vec;
    else
        _vel_sens = vec;
}



void
MAST::ElementBase::set_acceleration(const RealVectorX &vec,
                                    bool if_sens) {
    
    if (!if_sens)
        _accel = vec;
    else
        _accel_sens = vec;
}




//void
//MAST::ElementBase::set_base_solution(const RealVectorX& vec,
//                                     bool if_sens) {
//    if (!if_sens)
//        _base_sol = vec;
//    else
//        _base_sol_sens = vec;
//}




const libMesh::Elem&
MAST::ElementBase::get_elem_for_quadrature() const {
    
    return _local_elem->local_elem();
}



void
MAST::ElementBase::attach_active_solution_function(MAST::FunctionBase &f) {
    
    // make sure that this has not already been set
    libmesh_assert(!_active_sol_function);
    
    _active_sol_function = &f;
}


void
MAST::ElementBase::detach_active_solution_function() {
    
    _active_sol_function = NULL;
}



void
MAST::ElementBase::_init_fe_and_qrule(const libMesh::Elem& e,
                                      libMesh::FEBase **fe,
                                      libMesh::QBase **qrule,
                                      const std::vector<libMesh::Point>* pts) {
    
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.fetype(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.fetype(i));
    
    // Create an adequate quadrature rule
    (*fe) = libMesh::FEBase::build(e.dim(), fe_type).release();
    (*fe)->get_phi();
    (*fe)->get_JxW();
    (*fe)->get_dphi();
    (*fe)->get_dphidxi();
    (*fe)->get_dphideta();
    (*fe)->get_dphidzeta();
    
    if (pts == NULL) {
        (*qrule) = fe_type.default_quadrature_rule(e.dim(),
                                                _system.system().extra_quadrature_order).release();  // system extra quadrature
        (*fe)->attach_quadrature_rule(*qrule);
        (*fe)->reinit(&e);
    }
    else
        (*fe)->reinit(&e, pts);
}



void
MAST::ElementBase::_get_side_fe_and_qrule(const libMesh::Elem& e,
                                          unsigned int s,
                                          libMesh::FEBase **fe,
                                          libMesh::QBase **qrule,
                                          bool if_calculate_dphi) {
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.fetype(i));
    
    // Create an adequate quadrature rule
    (*fe)    = libMesh::FEBase::build(e.dim(), fe_type).release();
    (*qrule) = fe_type.default_quadrature_rule(e.dim()-1,
                                            _system.system().extra_quadrature_order).release();  // system extra quadrature
    (*fe)->attach_quadrature_rule(*qrule);
    (*fe)->get_phi();
    (*fe)->get_JxW();
    if (if_calculate_dphi)
        (*fe)->get_dphi();
    
    (*fe)->reinit(&e, s);
}

