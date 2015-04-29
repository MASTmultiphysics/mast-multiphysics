/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
_system(sys),
_elem(elem),
_active_sol_function(NULL),
_time(_system.system().time) {
    
}


MAST::ElementBase::~ElementBase()
{ }


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




void
MAST::ElementBase::set_base_solution(const RealVectorX& vec,
                                     bool if_sens) {
    if (!if_sens)
        _base_sol = vec;
    else
        _base_sol_sens = vec;
}




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
                                      const std::vector<libMesh::Point>* pts) {
    
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.fetype(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.fetype(i));
    
    // Create an adequate quadrature rule
    _fe.reset(libMesh::FEBase::build(e.dim(), fe_type).release());
    _fe->get_phi();
    _fe->get_JxW();
    _fe->get_dphi();
    _fe->get_dphidxi();
    _fe->get_dphideta();
    _fe->get_dphidzeta();
    
    if (pts == NULL) {
        _qrule.reset(fe_type.default_quadrature_rule
                     (e.dim(),
                      _system.system().extra_quadrature_order).release());  // system extra quadrature
        _fe->attach_quadrature_rule(_qrule.get());
        _fe->reinit(&e);
    }
    else
        _fe->reinit(&e, pts);
}



void
MAST::ElementBase::_get_side_fe_and_qrule(const libMesh::Elem& e,
                                          unsigned int s,
                                          std::auto_ptr<libMesh::FEBase>& fe,
                                          std::auto_ptr<libMesh::QBase>& qrule) {
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.fetype(0); // all variables are assumed to be of same type
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.fetype(i));
    
    // Create an adequate quadrature rule
    fe.reset(libMesh::FEBase::build(e.dim(), fe_type).release());
    qrule.reset(fe_type.default_quadrature_rule
                (e.dim()-1,
                 _system.system().extra_quadrature_order).release());  // system extra quadrature
    fe->attach_quadrature_rule(qrule.get());
    fe->get_phi();
    fe->get_JxW();
    
    fe->reinit(&e, s);
}

