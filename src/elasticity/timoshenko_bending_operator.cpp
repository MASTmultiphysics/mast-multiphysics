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
#include "elasticity/timoshenko_bending_operator.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "numerics/fem_operator_matrix.h"
#include "mesh/fe_base.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/field_function_base.h"


MAST::TimoshenkoBendingOperator::
TimoshenkoBendingOperator(MAST::StructuralElementBase& elem):
MAST::BendingOperator1D(elem),
_shear_quadrature_reduction(2)
{ }


MAST::TimoshenkoBendingOperator::~TimoshenkoBendingOperator() { }


void
MAST::TimoshenkoBendingOperator::
initialize_bending_strain_operator(const MAST::FEBase& fe,
                                   const unsigned int qp,
                                   MAST::FEMOperatorMatrix& Bmat_bend_v,
                                   MAST::FEMOperatorMatrix& Bmat_bend_w) {
    
    this->initialize_bending_strain_operator_for_yz(fe, qp, 1., 1.,
                                                    Bmat_bend_v,
                                                    Bmat_bend_w);
}




void
MAST::TimoshenkoBendingOperator::
initialize_bending_strain_operator_for_yz(const MAST::FEBase& fe,
                                          const unsigned int qp,
                                          const Real y,
                                          const Real z,
                                          MAST::FEMOperatorMatrix& Bmat_bend_v,
                                          MAST::FEMOperatorMatrix& Bmat_bend_w) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const std::vector<std::vector<Real> >& phi = fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi.size();
    
    RealVectorX phi_vec = RealVectorX::Zero(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    phi_vec *= -y;
    Bmat_bend_v.set_shape_function(0, 5, phi_vec); // v-bending: thetaz
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    phi_vec *= z;
    Bmat_bend_w.set_shape_function(0, 4, phi_vec); // w-bending : thetay
}



void
MAST::TimoshenkoBendingOperator::
calculate_transverse_shear_residual(bool request_jacobian,
                                    RealVectorX& local_f,
                                    RealMatrixX& local_jac) {
    
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::unique_ptr<MAST::FEBase>
    fe(_structural_elem.assembly().build_fe(_elem));
    fe->set_extra_quadrature_order(-_shear_quadrature_reduction);
    fe->init(_elem);
    
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >&                      phi = fe->get_phi();
    const std::vector<Real>&                                    JxW = fe->get_JxW();
    const std::vector<libMesh::Point>&                          xyz = fe->get_xyz();
    
    const unsigned int n_phi = (unsigned int)phi.size(), n2 = 6*n_phi;
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    
    RealVectorX
    vec_n2    = RealVectorX::Zero(n2),
    vec_2     = RealVectorX::Zero(2);
    RealMatrixX
    material_trans_shear_mat,
    mat_n2n2  = RealMatrixX::Zero(n2,n2),
    mat_2n2   = RealMatrixX::Zero(2,n2);
    
    
    FEMOperatorMatrix
    Bmat_trans_v,
    Bmat_trans_w;
    Bmat_trans_v.reinit(2, 6, n_phi); // one shear stress for v-bending
    Bmat_trans_w.reinit(2, 6, n_phi); // one shear stress for w-bending
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX>>
    mat_stiff = property.transverse_shear_stiffness_matrix(_structural_elem);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*mat_stiff)(xyz[qp],
                     _structural_elem.system().time,
                     material_trans_shear_mat);
        
        _transverse_shear_operations(phi,
                                     dphi,
                                     JxW,
                                     qp,
                                     material_trans_shear_mat,
                                     Bmat_trans_v,
                                     Bmat_trans_w,
                                     phi_vec,
                                     vec_n2,
                                     vec_2,
                                     mat_n2n2,
                                     mat_2n2,
                                     request_jacobian,
                                     local_f,
                                     local_jac);
    }
}



void
MAST::TimoshenkoBendingOperator::
calculate_transverse_shear_residual_sensitivity(const MAST::FunctionBase &p,
                                                bool request_jacobian,
                                                RealVectorX& local_f,
                                                RealMatrixX& local_jac) {
    
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::unique_ptr<MAST::FEBase>
    fe(_structural_elem.assembly().build_fe(_elem));
    fe->set_extra_quadrature_order(-_shear_quadrature_reduction);
    fe->init(_elem);
    
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >&                      phi = fe->get_phi();
    const std::vector<Real>&                                    JxW = fe->get_JxW();
    const std::vector<libMesh::Point>&                          xyz = fe->get_xyz();
    
    const unsigned int n_phi = (unsigned int)phi.size(), n2 = 6*n_phi;
    RealVectorX phi_vec  = RealVectorX::Zero(n_phi);
    
    RealVectorX
    vec_n2    = RealVectorX::Zero(n2),
    vec_2     = RealVectorX::Zero(2);
    RealMatrixX
    material_trans_shear_mat,
    mat_n2n2  = RealMatrixX::Zero(n2,n2),
    mat_2n2   = RealMatrixX::Zero(2,n2);
    
    
    FEMOperatorMatrix
    Bmat_trans_v,
    Bmat_trans_w;
    Bmat_trans_v.reinit(2, 6, n_phi); // one shear stress for v-bending
    Bmat_trans_w.reinit(2, 6, n_phi); // one shear stress for w-bending
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff = property.transverse_shear_stiffness_matrix(_structural_elem);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        mat_stiff->derivative(p,
                              xyz[qp],
                              _structural_elem.system().time,
                              material_trans_shear_mat);
        
        _transverse_shear_operations(phi,
                                     dphi,
                                     JxW,
                                     qp,
                                     material_trans_shear_mat,
                                     Bmat_trans_v,
                                     Bmat_trans_w,
                                     phi_vec,
                                     vec_n2,
                                     vec_2,
                                     mat_n2n2,
                                     mat_2n2,
                                     request_jacobian,
                                     local_f,
                                     local_jac);
    }
}




void
MAST::TimoshenkoBendingOperator::
_transverse_shear_operations(const std::vector<std::vector<Real> >& phi,
                             const std::vector<std::vector<libMesh::RealVectorValue> >& dphi,
                             const std::vector<Real>& JxW,
                             const unsigned int     qp,
                             const RealMatrixX&     material,
                             FEMOperatorMatrix&     Bmat_v,
                             FEMOperatorMatrix&     Bmat_w,
                             RealVectorX&           phi_vec,
                             RealVectorX&           vec_n2,
                             RealVectorX&           vec_2,
                             RealMatrixX&           mat_n2n2,
                             RealMatrixX&           mat_2n2,
                             bool                   request_jacobian,
                             RealVectorX&           local_f,
                             RealMatrixX&           local_jac) {

    // initialize the strain operator
    for ( unsigned int i_nd=0; i_nd<phi.size(); i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    Bmat_v.set_shape_function(0, 1, phi_vec); // gamma-xy:  v
    Bmat_w.set_shape_function(0, 2, phi_vec); // gamma-xz : w
    
    for ( unsigned int i_nd=0; i_nd<phi.size(); i_nd++ )
        phi_vec(i_nd) = phi[i_nd][qp];  // phi
    Bmat_w.set_shape_function(0, 4, phi_vec); // gamma-xz:  thetay
    phi_vec *= -1.0;
    Bmat_v.set_shape_function(0, 5, phi_vec); // gamma-xy : thetaz
    
    
    // now add the transverse shear component
    Bmat_v.vector_mult(vec_2, _structural_elem.local_solution());
    vec_2 = material * vec_2;
    Bmat_v.vector_mult_transpose(vec_n2, vec_2);
    local_f += JxW[qp] * vec_n2;
    
    Bmat_w.vector_mult(vec_2, _structural_elem.local_solution());
    vec_2 = material * vec_2;
    Bmat_w.vector_mult_transpose(vec_n2, vec_2);
    local_f += JxW[qp] * vec_n2;
    
    if (request_jacobian) {
        
        // now add the transverse shear component
        Bmat_v.left_multiply(mat_2n2, material);
        Bmat_v.right_multiply_transpose(mat_n2n2, mat_2n2);
        local_jac += JxW[qp] * mat_n2n2;
        
        Bmat_w.left_multiply(mat_2n2, material);
        Bmat_w.right_multiply_transpose(mat_n2n2, mat_2n2);
        local_jac += JxW[qp] * mat_n2n2;
    }
}


