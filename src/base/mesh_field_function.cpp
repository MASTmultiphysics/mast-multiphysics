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
#include "base/mesh_field_function.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/dof_map.h"


MAST::MeshFieldFunction::
MeshFieldFunction(MAST::SystemInitialization& sys,
                  const std::string& nm,
                  libMesh::ParallelType p_type):
MAST::FieldFunction<RealVectorX>(nm),
_use_qp_sol            (false),
_p_type                (p_type),
_qp_sol                (),
_sys                   (&sys.system()),
_function              (nullptr),
_perturbed_function    (nullptr)
{ }



MAST::MeshFieldFunction::
MeshFieldFunction(libMesh::System& sys,
                  const std::string& nm,
                  libMesh::ParallelType p_type):
MAST::FieldFunction<RealVectorX>(nm),
_use_qp_sol            (false),
_p_type                (p_type),
_qp_sol                (),
_sys                   (&sys),
_function              (nullptr),
_perturbed_function    (nullptr)
{ }




MAST::MeshFieldFunction::~MeshFieldFunction() {
 
    this->clear();
}







void
MAST::MeshFieldFunction::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealVectorX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(_function);
    
    unsigned int
    n_vars = _sys->n_vars();

    DenseRealVector v1;
    (*_function->_func)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);
    
    // now copy this to the output vector
    v = RealVectorX::Zero(n_vars);
    for (unsigned int i=0; i<n_vars; i++)
        v(i) = v1(i);
}



void
MAST::MeshFieldFunction::gradient (const libMesh::Point& p,
                                   const Real t,
                                   RealMatrixX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(_function);
    
    unsigned int
    n_vars = _sys->n_vars();
    
    std::vector<libMesh::Gradient> v1;
    _function->_func->gradient(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);
    
    // now copy this to the output vector
    v = RealMatrixX::Zero(n_vars, 3); // assume 3-dimensional by default
    for (unsigned int i=0; i<n_vars; i++)
        for (unsigned int j=0; j<3; j++)
            v(i, j) = v1[i](j);
}



void
MAST::MeshFieldFunction::perturbation(const libMesh::Point& p,
                                      const Real t,
                                      RealVectorX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(_perturbed_function);
    
    unsigned int
    n_vars = _sys->n_vars();

    DenseRealVector v1;
    (*_perturbed_function->_func)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);

    // now copy this to the output vector
    v = RealVectorX::Zero(n_vars);
    for (unsigned int i=0; i<n_vars; i++)
        v(i) = v1(i);
}


void
MAST::MeshFieldFunction::perturbation_gradient (const libMesh::Point& p,
                                                const Real t,
                                                RealMatrixX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(_function);
    
    unsigned int
    n_vars = _sys->n_vars();
    
    std::vector<libMesh::Gradient> v1;
    _perturbed_function->_func->gradient(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);
    
    // now copy this to the output vector
    v = RealMatrixX::Zero(n_vars, 3); // assume 3-dimensional by default
    for (unsigned int i=0; i<n_vars; i++)
        for (unsigned int j=0; j<3; j++)
            v(i, j) = v1[i](j);
}



void
MAST::MeshFieldFunction::derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealVectorX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the data for this function has been initialized
    std::map<const MAST::FunctionBase*, MAST::MeshFieldFunction::SolFunc*>::const_iterator
    it  = _function_sens.find(&f);
    
    // make sure that the object was initialized
    libmesh_assert(it != _function_sens.end());
    
    unsigned int
    n_vars = _sys->n_vars();

    DenseRealVector v1;
    (*it->second->_func)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);

    // now copy this to the output vector
    v = RealVectorX::Zero(n_vars);
    for (unsigned int i=0; i<n_vars; i++)
        v(i) = v1(i);
}


void
MAST::MeshFieldFunction::derivative_gradient (const MAST::FunctionBase& f,
                                              const libMesh::Point& p,
                                              const Real t,
                                              RealMatrixX& v) const {
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (_use_qp_sol) {
        v = _qp_sol;
        return;
    }
    
    // make sure that the data for this function has been initialized
    std::map<const MAST::FunctionBase*, MAST::MeshFieldFunction::SolFunc*>::const_iterator
    it  = _function_sens.find(&f);

    // make sure that the object was initialized
    libmesh_assert(it != _function_sens.end());

    unsigned int
    n_vars = _sys->n_vars();
    
    std::vector<libMesh::Gradient> v1;
    it->second->_func->gradient(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert_equal_to(v1.size(), n_vars);
    
    // now copy this to the output vector
    v = RealMatrixX::Zero(n_vars, 3); // assume 3-dimensional by default
    for (unsigned int i=0; i<n_vars; i++)
        for (unsigned int j=0; j<3; j++)
            v(i, j) = v1[i](j);
}



void
MAST::MeshFieldFunction::init(const libMesh::NumericVector<Real>& sol,
                              bool reuse_vector) {
    
    // first make sure that the object is not already initialized
    libmesh_assert(!_function);
    
    _function = new MAST::MeshFieldFunction::SolFunc;
    _init_sol_func(reuse_vector, sol, *_function);
}



void
MAST::MeshFieldFunction::init_sens(const MAST::FunctionBase& f,
                                   const libMesh::NumericVector<Real>& sol,
                                   bool reuse_vector) {
    
    // make sure the function has not already been initialized
    std::map<const MAST::FunctionBase*, MAST::MeshFieldFunction::SolFunc*>::const_iterator
    it  = _function_sens.find(&f);

    // make sure that the object was initialized
    libmesh_assert(it == _function_sens.end());
    
    MAST::MeshFieldFunction::SolFunc* func = new MAST::MeshFieldFunction::SolFunc;
    _init_sol_func(reuse_vector, sol, *func);
    
    _function_sens[&f] = func;
}




void
MAST::MeshFieldFunction::clear() {
    
    // if a pointer has been attached, then delete it and the
    // associated vector, and clear the associated system
    if (_function) {
        delete _function;
        _function = nullptr;
    }
    
    if (_perturbed_function) {
        delete _perturbed_function;
        _perturbed_function = nullptr;
    }

    // now clear all the sensitivity data
    std::map<const MAST::FunctionBase*, MAST::MeshFieldFunction::SolFunc*>::const_iterator
    it  = _function_sens.begin(),
    end = _function_sens.end();
    
    for ( ; it != end; it++)
        delete it->second;
    
    _function_sens.clear();

    // clear flags for quadrature point solution
    _use_qp_sol = false;
}




void
MAST::MeshFieldFunction::
set_element_quadrature_point_solution(RealVectorX& sol) {
    
    _use_qp_sol = true;
    _qp_sol     = sol;
}



void
MAST::MeshFieldFunction::
clear_element_quadrature_point_solution() {
    
    _use_qp_sol = false;
    _qp_sol.setZero();
}


void
MAST::MeshFieldFunction::_init_sol_func(bool reuse_sol,
                                        const libMesh::NumericVector<Real> &sol,
                                        MAST::MeshFieldFunction::SolFunc& sol_func) {
    
    // make sure it has not already been initialized
    libmesh_assert(!sol_func._func);
        
    if (reuse_sol) {
        
        // make sure the given vector is of the same type as that specified for this object
        libmesh_assert_equal_to(sol.type(), _p_type);
        sol_func._sol       = &sol;
    }
    else  {
        
        sol_func._cloned_sol = libMesh::NumericVector<Real>::build(_sys->comm()).release();
        sol_func._sol        = sol_func._cloned_sol;
        
        switch (_p_type) {
                
            case libMesh::SERIAL: {
                sol_func._cloned_sol->init(sol.size(), true, libMesh::SERIAL);
                sol.localize(*sol_func._cloned_sol);
            }
                break;
                
            case libMesh::GHOSTED: {
                
                sol_func._cloned_sol->init(_sys->n_dofs(),
                                           _sys->n_local_dofs(),
                                           _sys->get_dof_map().get_send_list(),
                                           true,
                                           libMesh::GHOSTED);
                sol.localize(*sol_func._cloned_sol, _sys->get_dof_map().get_send_list());
            }
                break;
                
            default:
                // not implemented for other types.
                libmesh_error();
                break;
        }
    }

    // finally, create the mesh interpolation function
    std::vector<unsigned int> vars;
    _sys->get_all_variable_numbers(vars);
    
    sol_func._func = new libMesh::MeshFunction(_sys->get_equation_systems(),
                                                   *sol_func._sol,
                                                   _sys->get_dof_map(),
                                                   vars);
    sol_func._func->init();
}
