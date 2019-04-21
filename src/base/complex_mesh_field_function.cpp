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
#include "base/complex_mesh_field_function.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"


MAST::ComplexMeshFieldFunction::
ComplexMeshFieldFunction(MAST::SystemInitialization& sys,
                         const std::string& nm):
MAST::FieldFunction<ComplexVectorX>(nm),
_system(&sys),
_sol_re(nullptr),
_sol_im(nullptr),
_perturbed_sol_re(nullptr),
_perturbed_sol_im(nullptr),
_function_re(nullptr),
_function_im(nullptr),
_perturbed_function_re(nullptr),
_perturbed_function_im(nullptr)
{ }




MAST::ComplexMeshFieldFunction::~ComplexMeshFieldFunction() {
    
    this->clear();
}





void
MAST::ComplexMeshFieldFunction::
init(const libMesh::NumericVector<Real>& sol_re,
     const libMesh::NumericVector<Real>& sol_im) {
    
    // first make sure that the object is not already initialized
    libmesh_assert(!_function_re);
    
    MAST::NonlinearSystem& system = _system->system();
    
    // next, clone this solution and localize to the sendlist
    _sol_re = libMesh::NumericVector<Real>::build(system.comm()).release();
    _sol_im = libMesh::NumericVector<Real>::build(system.comm()).release();
    
    const std::vector<libMesh::dof_id_type>& send_list =
    system.get_dof_map().get_send_list();
    
    // initialize and then localize the vector with the provided solution
    /*_sol_re->init(system.n_dofs(),
                  system.n_local_dofs(),
                  send_list,
                  false,
                  libMesh::GHOSTED);
    sol_re.localize(*_sol_re, send_list);
    
    _sol_im->init(system.n_dofs(),
                  system.n_local_dofs(),
                  send_list,
                  false,
                  libMesh::GHOSTED);
    sol_im.localize(*_sol_im, send_list);*/

    _sol_re->init(sol_re.size(), true, libMesh::SERIAL);
    sol_re.localize(*_sol_re);
    _sol_im->init(sol_im.size(), true, libMesh::SERIAL);
    sol_im.localize(*_sol_im);

    
    // finally, create the mesh interpolation function
    _function_re = new libMesh::MeshFunction(system.get_equation_systems(),
                                             *_sol_re,
                                             system.get_dof_map(),
                                             _system->vars());
    _function_re->init();
    
    _function_im = new libMesh::MeshFunction(system.get_equation_systems(),
                                             *_sol_im,
                                             system.get_dof_map(),
                                             _system->vars());
    _function_im->init();
}




void
MAST::ComplexMeshFieldFunction::
init_perturbation(const libMesh::NumericVector<Real>& sol_re,
                  const libMesh::NumericVector<Real>& sol_im) {
    
    // first make sure that the object is not already initialized
    libmesh_assert(!_perturbed_function_re);
    
    MAST::NonlinearSystem& system = _system->system();
    
    // next, clone this solution and localize to the sendlist
    _perturbed_sol_re = libMesh::NumericVector<Real>::build(system.comm()).release();
    _perturbed_sol_im = libMesh::NumericVector<Real>::build(system.comm()).release();
    
    const std::vector<libMesh::dof_id_type>& send_list =
    system.get_dof_map().get_send_list();
    
    // initialize and then localize the vector with the provided solution
    /*_perturbed_sol_re->init(system.n_dofs(),
                            system.n_local_dofs(),
                            send_list,
                            false,
                            libMesh::GHOSTED);
    sol_re.localize(*_sol_re, send_list);
    
    _perturbed_sol_im->init(system.n_dofs(),
                            system.n_local_dofs(),
                            send_list,
                            false,
                            libMesh::GHOSTED);
    sol_im.localize(*_perturbed_sol_im, send_list);*/
    _perturbed_sol_re->init(sol_re.size(), true, libMesh::SERIAL);
    sol_re.localize(*_perturbed_sol_re);
    _perturbed_sol_im->init(sol_im.size(), true, libMesh::SERIAL);
    sol_im.localize(*_perturbed_sol_im);
    
    
    // finally, create the mesh interpolation function
    _perturbed_function_re = new libMesh::MeshFunction(system.get_equation_systems(),
                                                       *_perturbed_sol_re,
                                                       system.get_dof_map(),
                                                       _system->vars());
    _perturbed_function_re->init();
    
    _perturbed_function_im = new libMesh::MeshFunction(system.get_equation_systems(),
                                                       *_perturbed_sol_im,
                                                       system.get_dof_map(),
                                                       _system->vars());
    _perturbed_function_im->init();
}




void
MAST::ComplexMeshFieldFunction::operator() (const libMesh::Point& p,
                                            const Real t,
                                            ComplexVectorX& v) const {
    
    // make sure that the object was initialized
    libmesh_assert(_function_re);
    
    DenseRealVector v_re, v_im;
    (*_function_re)(p, t, v_re);
    (*_function_im)(p, t, v_im);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v_re.size());
    
    // now copy this to the output vector
    v = ComplexVectorX::Zero(v_re.size());
    for (unsigned int i=0; i<v_re.size(); i++)
        v(i) = std::complex<Real>(v_re(i), v_im(i));
}





void
MAST::ComplexMeshFieldFunction::perturbation(const libMesh::Point& p,
                                             const Real t,
                                             ComplexVectorX& v) const {
    
    // make sure that the object was initialized
    libmesh_assert(_perturbed_function_re);
    
    DenseRealVector v_re, v_im;
    (*_perturbed_function_re)(p, t, v_re);
    (*_perturbed_function_im)(p, t, v_im);
    
    // make sure that the mesh function was able to find the element
    // and a solution
    libmesh_assert(v_re.size());
    
    // now copy this to the output vector
    v = ComplexVectorX::Zero(v_re.size());
    for (unsigned int i=0; i<v_re.size(); i++)
        v(i) = std::complex<Real>(v_re(i), v_im(i));
}








void
MAST::ComplexMeshFieldFunction::clear() {
    
    // if a pointer has been attached, then delete it and the
    // associated vector, and clear the associated system
    if (_function_re) {
        delete _function_re;
        delete _function_im;
        _function_re = nullptr;
        _function_im = nullptr;
        
        delete _sol_re;
        delete _sol_im;
        _sol_re = nullptr;
        _sol_im = nullptr;
        
    }
    
    if (_perturbed_function_re) {
        delete _perturbed_function_re;
        delete _perturbed_function_im;
        _perturbed_function_re = nullptr;
        _perturbed_function_im = nullptr;
        
        delete _perturbed_sol_re;
        delete _perturbed_sol_im;
        _perturbed_sol_re = nullptr;
        _perturbed_sol_im = nullptr;
    }
}




