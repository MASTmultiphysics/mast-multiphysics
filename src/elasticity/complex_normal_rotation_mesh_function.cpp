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
#include "elasticity/complex_normal_rotation_mesh_function.h"
#include "base/complex_mesh_field_function.h"


MAST::ComplexNormalRotationMeshFunction::
ComplexNormalRotationMeshFunction(const std::string& nm,
                                  MAST::ComplexMeshFieldFunction& func):
MAST::NormalRotationFunctionBase<ComplexVectorX>(nm),
_func(func)
{ }



void
MAST::ComplexNormalRotationMeshFunction::operator()(const libMesh::Point& p,
                                                    const libMesh::Point& n,
                                                    const Real t,
                                                    ComplexVectorX& dn_rot) const {
    
    
    dn_rot.setZero();
    
    libMesh::MeshFunction
    &function_re = *_func.get_function().first,
    &function_im = *_func.get_function().second;
    
    // translation is obtained by direct interpolation of the u,v,w vars
    
    DenseRealVector
    v_re,
    v_im;
    
    function_re(p, 0., v_re);
    function_re(p, 0., v_im);
    
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient>
    gradients_re,
    gradients_im;
    function_re.gradient(p, 0., gradients_re);
    function_im.gradient(p, 0., gradients_im);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    ComplexVectorX
    rot = ComplexVectorX::Zero(3);
    
    rot(0) = std::complex<Real>(gradients_re[2](1) - gradients_re[1](2),
                                gradients_im[2](1) - gradients_im[1](2)); // dwz/dy - dwy/dz
    rot(1) = std::complex<Real>(gradients_re[0](2) - gradients_re[2](0),
                                gradients_im[0](2) - gradients_im[2](0)); // dwx/dz - dwz/dx
    rot(2) = std::complex<Real>(gradients_re[1](0) - gradients_re[0](1),
                                gradients_im[1](0) - gradients_im[0](1)); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}



void
MAST::ComplexNormalRotationMeshFunction::perturbation(const libMesh::Point& p,
                                                      const libMesh::Point& n,
                                                      const Real t,
                                                      ComplexVectorX& dn_rot) const {
    
    
    dn_rot.setZero();
    
    libMesh::MeshFunction
    &perturbed_function_re = *_func.get_perturbed_function().first,
    &perturbed_function_im = *_func.get_perturbed_function().second;
    
    // translation is obtained by direct interpolation of the u,v,w vars
    
    DenseRealVector
    v_re,
    v_im;
    perturbed_function_re(p, 0., v_re);
    perturbed_function_im(p, 0., v_im);
    
    // perturbation of the normal requires calculation of the curl of
    // displacement at the given point
    std::vector<libMesh::Gradient>
    gradients_re,
    gradients_im;
    perturbed_function_re.gradient(p, 0., gradients_re);
    perturbed_function_im.gradient(p, 0., gradients_im);
    
    // TODO: these need to be mapped from local 2D to 3D space
    
    // now prepare the rotation vector
    ComplexVectorX
    rot = ComplexVectorX::Zero(3);
    
    rot(0) = std::complex<Real>(gradients_re[2](1) - gradients_re[1](2),
                                gradients_im[2](1) - gradients_im[1](2)); // dwz/dy - dwy/dz
    rot(1) = std::complex<Real>(gradients_re[0](2) - gradients_re[2](0),
                                gradients_im[0](2) - gradients_im[2](0)); // dwx/dz - dwz/dx
    rot(2) = std::complex<Real>(gradients_re[1](0) - gradients_re[0](1),
                                gradients_im[1](0) - gradients_im[0](1)); // dwy/dx - dwx/dy
    
    // now do the cross-products
    dn_rot(0) =   rot(1) * n(2) - rot(2) * n(1);
    dn_rot(1) = -(rot(0) * n(2) - rot(2) * n(0));
    dn_rot(2) =   rot(0) * n(1) - rot(1) * n(0);
}


