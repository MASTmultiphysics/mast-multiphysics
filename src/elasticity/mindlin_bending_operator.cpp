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
#include "elasticity/mindlin_bending_operator.h"
#include "elasticity/structural_element_base.h"
#include "property_cards/element_property_card_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/nonlinear_system.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"
#include "base/assembly_base.h"
#include "base/field_function_base.h"



MAST::MindlinBendingOperator::
MindlinBendingOperator(MAST::StructuralElementBase& elem):
MAST::BendingOperator2D(elem),
_shear_quadrature_reduction(2)
{ }



MAST::MindlinBendingOperator::~MindlinBendingOperator() { }



void
MAST::MindlinBendingOperator::
initialize_bending_strain_operator(const MAST::FEBase& fe,
                                   const unsigned int qp,
                                   MAST::FEMOperatorMatrix& Bmat_bend) {
    
    this->initialize_bending_strain_operator_for_z(fe, qp, 1., Bmat_bend);
}



void
MAST::MindlinBendingOperator::
initialize_bending_strain_operator_for_z(const MAST::FEBase& fe,
                                         const unsigned int qp,
                                         const Real z,
                                         MAST::FEMOperatorMatrix& Bmat_bend) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const std::vector<std::vector<Real> >& phi = fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi.size();
    
    RealVectorX phi_vec = RealVectorX::Zero(n_phi);
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](0);  // dphi/dx
    
    phi_vec   *= z;
    Bmat_bend.set_shape_function(0, 4, phi_vec); // epsilon-x: thetay
    phi_vec   *= -1.0;
    Bmat_bend.set_shape_function(2, 3, phi_vec); // gamma-xy : thetax
    
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
    
    phi_vec   *= z;
    Bmat_bend.set_shape_function(2, 4, phi_vec); // gamma-xy : thetay
    //Bmat_trans.set_shape_function(1, 2, phi_vec); // gamma-yz : w
    phi_vec   *= -1.0;
    Bmat_bend.set_shape_function(1, 3, phi_vec); // epsilon-y: thetax
}




void
MAST::MindlinBendingOperator::
calculate_transverse_shear_residual(bool request_jacobian,
                                    RealVectorX& local_f,
                                    RealMatrixX& local_jac) {
    
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false, -_shear_quadrature_reduction));

    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n2    = 6*n_phi;
    
    RealVectorX
    phi_vec   = RealVectorX::Zero(n_phi),
    vec_n2    = RealVectorX::Zero(n2),
    vec_2     = RealVectorX::Zero(2);
    RealMatrixX
    material_trans_shear_mat,
    mat_n2n2    = RealMatrixX::Zero(n2,n2),
    mat_2n2     = RealMatrixX::Zero(2,n2);
    
    
    FEMOperatorMatrix Bmat_trans;
    Bmat_trans.reinit(2, 6, n_phi); // only two shear stresses
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
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
                                     Bmat_trans,
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
MAST::MindlinBendingOperator::
calculate_transverse_shear_residual_sensitivity(const MAST::FunctionBase& p,
                                                bool request_jacobian,
                                                RealVectorX& local_f,
                                                RealMatrixX& local_jac) {
    
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false, -_shear_quadrature_reduction));
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n2    = 6*n_phi;
    
    RealVectorX
    phi_vec   = RealVectorX::Zero(n_phi),
    vec_n2    = RealVectorX::Zero(n2),
    vec_2     = RealVectorX::Zero(2);
    RealMatrixX
    material_trans_shear_mat,
    mat_n2n2    = RealMatrixX::Zero(n2,n2),
    mat_2n2     = RealMatrixX::Zero(2,n2);

    
    FEMOperatorMatrix Bmat_trans;
    Bmat_trans.reinit(2, 6, n_phi); // only two shear stresses
    
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
                                     Bmat_trans,
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
MAST::MindlinBendingOperator::
calculate_transverse_shear_residual_boundary_velocity
(const MAST::FunctionBase& p,
 const unsigned int s,
 const MAST::FieldFunction<RealVectorX>& vel_f,
 bool request_jacobian,
 RealVectorX& local_f,
 RealMatrixX& local_jac) {
    
    const MAST::ElementPropertyCardBase& property = _structural_elem.elem_property();
    
    // make an fe and quadrature object for the requested order for integrating
    // transverse shear
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_side_fe(s, true, -_shear_quadrature_reduction));

    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe->get_dphi();
    const std::vector<std::vector<Real> >& phi                      = fe->get_phi();
    const std::vector<libMesh::Point>& xyz                          = fe->get_xyz();
    const std::vector<libMesh::Point>& face_normals                 = fe->get_normals_for_local_coordinate();
    std::vector<Real> JxW_Vn                                        = fe->get_JxW();

    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n2    = 6*n_phi,
    dim   = 2;
    
    RealVectorX
    phi_vec     = RealVectorX::Zero(n_phi),
    vec_n2      = RealVectorX::Zero(n2),
    vec_2       = RealVectorX::Zero(2),
    vel         = RealVectorX::Zero(dim);
    RealMatrixX
    material_trans_shear_mat,
    dmaterial_trans_shear_mat_dp,
    mat_n2n2    = RealMatrixX::Zero(n2,n2),
    mat_2n2     = RealMatrixX::Zero(2,n2);

    
    FEMOperatorMatrix Bmat_trans;
    Bmat_trans.reinit(2, 6, n_phi); // only two shear stresses
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff = property.transverse_shear_stiffness_matrix(_structural_elem);
    
    
    Real
    vn  = 0.;
    
    // modify the JxW_Vn by multiplying the normal velocity to it
    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        vel_f.derivative(p, xyz[qp], _structural_elem.system().time, vel);
        vn = 0.;
        for (unsigned int i=0; i<dim; i++)
            vn += vel(i)*face_normals[qp](i);
        JxW_Vn[qp] *= vn;
    }

    for (unsigned int qp=0; qp<JxW_Vn.size(); qp++) {
        
        (*mat_stiff)(xyz[qp],
                     _structural_elem.system().time,
                     material_trans_shear_mat);

        _transverse_shear_operations(phi,
                                     dphi,
                                     JxW_Vn,
                                     qp,
                                     material_trans_shear_mat,
                                     Bmat_trans,
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
MAST::MindlinBendingOperator::
_transverse_shear_operations(const std::vector<std::vector<Real> >& phi,
                             const std::vector<std::vector<libMesh::RealVectorValue> >& dphi,
                             const std::vector<Real>& JxW,
                             const unsigned int     qp,
                             const RealMatrixX&     material,
                             FEMOperatorMatrix&     Bmat,
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
    
    Bmat.set_shape_function(0, 2, phi_vec); // gamma-xz:  w
    
    for ( unsigned int i_nd=0; i_nd<phi.size(); i_nd++ )
        phi_vec(i_nd) = dphi[i_nd][qp](1);  // dphi/dy
    
    Bmat.set_shape_function(1, 2, phi_vec); // gamma-yz : w
    
    for ( unsigned int i_nd=0; i_nd<phi.size(); i_nd++ )
        phi_vec(i_nd) = phi[i_nd][qp];  // phi
    
    Bmat.set_shape_function(0, 4, phi_vec); // gamma-xz:  thetay
    phi_vec  *= -1.0;
    Bmat.set_shape_function(1, 3, phi_vec); // gamma-yz : thetax
    
    
    // now add the transverse shear component
    Bmat.vector_mult(vec_2, _structural_elem.local_solution());
    vec_2 = material * vec_2;
    Bmat.vector_mult_transpose(vec_n2, vec_2);
    local_f += JxW[qp] * vec_n2;
    
    if (request_jacobian) {
        
        // now add the transverse shear component
        Bmat.left_multiply(mat_2n2, material);
        Bmat.right_multiply_transpose(mat_n2n2, mat_2n2);
        local_jac += JxW[qp] * mat_n2n2;
    }
}
