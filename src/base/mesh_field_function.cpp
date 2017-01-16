/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2017  Manav Bhatia
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

// libMesh includes
#include "libmesh/dof_map.h"


MAST::MeshFieldFunction::
MeshFieldFunction(MAST::SystemInitialization& sys,
                  const std::string& nm):
MAST::FieldFunction<RealVectorX>(nm),
_use_qp_sol(false),
_qp_sol(),
_system(&sys),
_sol(nullptr),
_function(nullptr)
{ }




MAST::MeshFieldFunction::~MeshFieldFunction() {
 
    this->clear();
}







void
MAST::MeshFieldFunction::operator() (const libMesh::Point& p,
                                     const Real t,
                                     RealVectorX& v) const {
    
    // get the pointer to the mesh function from the master
    const MAST::MeshFieldFunction* master =
    dynamic_cast<const MAST::MeshFieldFunction*>(this);
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (master->_use_qp_sol) {
        v = master->_qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(master->_function);
    
    DenseRealVector v1;
    (*master->_function)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v1.size());
    
    // now copy this to the output vector
    v = RealVectorX::Zero(v1.size());
    for (unsigned int i=0; i<v1.size(); i++)
        v(i) = v1(i);
}


void
MAST::MeshFieldFunction::perturbation(const libMesh::Point& p,
                                      const Real t,
                                      RealVectorX& v) const {
    
    // get the pointer to the mesh function from the master
    const MAST::MeshFieldFunction* master =
    dynamic_cast<const MAST::MeshFieldFunction*>(this);
    
    // if the element has provided a quadrature point solution,
    // then use it
    if (master->_use_qp_sol) {
        v = master->_qp_sol;
        return;
    }
    
    // make sure that the object was initialized
    libmesh_assert(master->_perturbed_function);
    
    DenseRealVector v1;
    (*master->_perturbed_function)(p, t, v1);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v1.size());
    
    // now copy this to the output vector
    v = RealVectorX::Zero(v1.size());
    for (unsigned int i=0; i<v1.size(); i++)
        v(i) = v1(i);
}




void
MAST::MeshFieldFunction::derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     RealVectorX& v) const {
    
    
    libmesh_error_msg("To be implemented.");
}





void
MAST::MeshFieldFunction::
init(const libMesh::NumericVector<Real>& sol,
     const libMesh::NumericVector<Real>* dsol) {
    
    // this is to be done for the master only
    const MAST::MeshFieldFunction* master_const =
    dynamic_cast<const MAST::MeshFieldFunction*>(this);
    MAST::MeshFieldFunction* master =
    const_cast<MAST::MeshFieldFunction*>(master_const);
    
    
    // first make sure that the object is not already initialized
    libmesh_assert(!master->_function);
    
    libMesh::System& system = _system->system();
    
    // next, clone this solution and localize to the sendlist
    master->_sol = libMesh::NumericVector<Real>::build(system.comm()).release();

    const std::vector<libMesh::dof_id_type>& send_list =
    system.get_dof_map().get_send_list();
    
    // initialize and then localize the vector with the provided solution
    master->_sol->init(system.n_dofs(),
                       system.n_local_dofs(),
                       send_list,
                       false,
                       libMesh::GHOSTED);
    sol.localize(*master->_sol, send_list);
    
    // finally, create the mesh interpolation function
    master->_function = new libMesh::MeshFunction(system.get_equation_systems(),
                                                       *master->_sol,
                                                       system.get_dof_map(),
                                                       _system->vars());
    master->_function->init();
    
    if (dsol) {

        master->_dsol->init(system.n_dofs(),
                            system.n_local_dofs(),
                            send_list,
                            false,
                            libMesh::GHOSTED);
        dsol->localize(*master->_dsol, send_list);
        
        // finally, create the mesh interpolation function
        master->_perturbed_function =
        new libMesh::MeshFunction(system.get_equation_systems(),
                                  *master->_dsol,
                                  system.get_dof_map(),
                                  _system->vars());
        master->_perturbed_function->init();

    }
}




void
MAST::MeshFieldFunction::clear() {
    
    //only the master function will call this on itself
    if (this != this)
        return;
    
    // if a pointer has been attached, then delete it and the
    // associated vector, and clear the associated system
    if (_function) {
        delete _function;
        _function = nullptr;
        
        delete _sol;
        _sol = nullptr;
        
        _system = nullptr;
    }

    if (_perturbed_function) {
        delete _perturbed_function;
        _perturbed_function = nullptr;
        
        delete _dsol;
        _dsol = nullptr;
    }

    
    // clear flags for quadrature point solution
    _use_qp_sol = false;
}




void
MAST::MeshFieldFunction::
set_element_quadrature_point_solution(RealVectorX& sol) {
    
    // this is to be done for the master only
    const MAST::MeshFieldFunction* master_const =
    dynamic_cast<const MAST::MeshFieldFunction*>(this);
    MAST::MeshFieldFunction* master =
    const_cast<MAST::MeshFieldFunction*>(master_const);

    master->_use_qp_sol = true;
    master->_qp_sol     = sol;
}



void
MAST::MeshFieldFunction::
clear_element_quadrature_point_solution() {

    // this is to be done for the master only
    const MAST::MeshFieldFunction* master_const =
    dynamic_cast<const MAST::MeshFieldFunction*>(this);
    MAST::MeshFieldFunction* master =
    const_cast<MAST::MeshFieldFunction*>(master_const);
    
    master->_use_qp_sol = false;
    master->_qp_sol.setZero();
}

