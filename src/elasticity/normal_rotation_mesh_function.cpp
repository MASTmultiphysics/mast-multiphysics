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
#include "elasticity/normal_rotation_mesh_function.h"
#include "base/mesh_field_function.h"


MAST::NormalRotationMeshFunction::
NormalRotationMeshFunction(const std::string& nm,
                           MAST::MeshFieldFunction& func):
MAST::NormalRotationFunctionBase<RealVectorX>(nm),
_func(func)
{ }



void
MAST::NormalRotationMeshFunction::operator()(const libMesh::Point& p,
                                             const libMesh::Point& n,
                                             const Real t,
                                             RealVectorX& dn_rot) const {
    
    
    dn_rot.setZero();
    dn_rot.setZero();
    
    libMesh::MeshFunction&
    function = _func.get_function();
    
    // translation is obtained by direct interpolation of the u,v,w vars
    
    DenseRealVector v;
    function(p, 0., v);
    
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient> gradients;
    function.gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    RealVectorX
    rot = RealVectorX::Zero(3);
    
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}



void
MAST::NormalRotationMeshFunction::perturbation(const libMesh::Point& p,
                                               const libMesh::Point& n,
                                               const Real t,
                                               RealVectorX& dn_rot) const {
    
    
    dn_rot.setZero();
    dn_rot.setZero();
    
    libMesh::MeshFunction&
    perturbed_function = _func.get_perturbed_function();
    
    // translation is obtained by direct interpolation of the u,v,w vars
    
    DenseRealVector v;
    perturbed_function(p, 0., v);
    
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient> gradients;
    perturbed_function.gradient(p, 0., gradients);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    RealVectorX
    rot = RealVectorX::Zero(3);
    
    rot(0) = gradients[2](1) - gradients[1](2); // dwz/dy - dwy/dz
    rot(1) = gradients[0](2) - gradients[2](0); // dwx/dz - dwz/dx
    rot(2) = gradients[1](0) - gradients[0](1); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}


