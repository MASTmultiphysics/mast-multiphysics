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
#include "elasticity/solid_element_3d.h"
#include "elasticity/stress_output_base.h"
#include "base/boundary_condition_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"
#include "base/assembly_base.h"
#include "base/field_function_base.h"
#include "numerics/fem_operator_matrix.h"
#include "mesh/geom_elem.h"
#include "mesh/fe_base.h"
#include "property_cards/element_property_card_base.h"


MAST::StructuralElement3D::
StructuralElement3D(MAST::SystemInitialization& sys,
                    const MAST::GeomElem& elem,
                    const MAST::ElementPropertyCardBase& p):
MAST::StructuralElementBase(sys, elem, p) {
    
}


MAST::StructuralElement3D::~StructuralElement3D() {
    
}


bool
MAST::StructuralElement3D::inertial_residual (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac_xddot,
                                              RealMatrixX& jac_xdot,
                                              RealMatrixX& jac) {
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false));
    
    const std::vector<Real>& JxW               = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz     = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    
    const unsigned int
    n_phi    = (unsigned int)phi.size(),
    n1       =3,
    n2       =3*n_phi;
    
    RealMatrixX
    material_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1, n2),
    mat2_n2n2     = RealMatrixX::Zero(n2, n2);
    RealVectorX
    phi_vec    = RealVectorX::Zero(n_phi),
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n2    = RealVectorX::Zero(n2),
    local_acc  = RealVectorX::Zero(n2);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    mat_inertia  = _property.inertia_matrix(*this);
    
    MAST::FEMOperatorMatrix Bmat;
    
    local_acc.topRows(n2) = _local_accel.topRows(n2);
    
    if (_property.if_diagonal_mass_matrix()) {
        
        (*mat_inertia)(xyz[0], _time, material_mat);
        
        Real vol = 0.;
        const unsigned int nshp = fe->n_shape_functions();
        for (unsigned int i=0; i<JxW.size(); i++)
            vol += JxW[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<3; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                jac_xddot(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        
        f.topRows(n2) =  jac_xddot.topLeftCorner(n2, n2) * _local_accel.topRows(n2);
    }
    else {
        libMesh::Point p;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            (*mat_inertia)(xyz[0], _time, material_mat);
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(3, phi_vec);
            
            Bmat.left_multiply(mat1_n1n2, material_mat);
            
            vec1_n1 = mat1_n1n2 * local_acc;
            Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
            
            f.topRows(n2) += JxW[qp] * vec2_n2;
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(mat2_n2n2,
                                              mat1_n1n2);
                jac_xddot.topLeftCorner(n2, n2) += JxW[qp]*mat2_n2n2;
            }
        }
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::internal_residual(bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac) {

    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW            = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = fe->get_xyz();
    const unsigned int
    n_phi              = (unsigned int)fe->n_shape_functions(),
    n1                 =6,
    n2                 =3*n_phi,
    n3                 =30,
    n_nodes            =_elem.get_reference_elem().n_nodes();
    
    RealMatrixX
    material_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    mat1_n1n2    = RealMatrixX::Zero(n1, n2),
    mat2_n2n2    = RealMatrixX::Zero(n2, n2),
    mat3_3n2     = RealMatrixX::Zero(3, n2),
    mat4_33      = RealMatrixX::Zero(3, 3),
    mat5_n1n3    = RealMatrixX::Zero(n1, n3),
    mat6_n2n3    = RealMatrixX::Zero(n2, n3),
    mat7_3n3     = RealMatrixX::Zero(3, n3),
    Gmat         = RealMatrixX::Zero(6, n3),
    K_alphaalpha = RealMatrixX::Zero(n3, n3),
    K_ualpha     = RealMatrixX::Zero(n2, n3),
    K_corr       = RealMatrixX::Zero(n2, n2);
    RealVectorX
    strain    = RealVectorX::Zero(6),
    stress    = RealVectorX::Zero(6),
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    vec3_3    = RealVectorX::Zero(3),
    local_disp= RealVectorX::Zero(n2),
    f_alpha   = RealVectorX::Zero(n3),
    alpha     = RealVectorX::Zero(n3);//*_incompatible_sol;
    
    // copy the values from the global to the local element
    local_disp.topRows(n2) = _local_sol.topRows(n2);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w,
    Bmat_inc;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, n_nodes);
    Bmat_nl_x.reinit(3, 3, n_nodes);
    Bmat_nl_y.reinit(3, 3, n_nodes);
    Bmat_nl_z.reinit(3, 3, n_nodes);
    Bmat_nl_u.reinit(3, 3, n_nodes);
    Bmat_nl_v.reinit(3, 3, n_nodes);
    Bmat_nl_w.reinit(3, 3, n_nodes);
    Bmat_inc.reinit(n1, n3, 1);            // six stress-strain components

    /*
    // initialize the incompatible mode mapping at element mid-point
    _init_incompatible_fe_mapping(_elem);
    
    ///////////////////////////////////////////////////////////////////
    // first for loop to evaluate alpha
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *_fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        this->initialize_incompatible_strain_operator(qp, *_fe, Bmat_inc, Gmat);
        
        // calculate the incompatible mode matrices
        // incompatible mode diagonal stiffness matrix
        mat5_n1n3    =  material_mat * Gmat;
        K_alphaalpha += JxW[qp] * ( Gmat.transpose() * mat5_n1n3);
        
        // off-diagonal coupling matrix
        // linear strain term
        Bmat_lin.right_multiply_transpose(mat6_n2n3, mat5_n1n3);
        K_ualpha  += JxW[qp] * mat6_n2n3;
        
        if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
            
            // nonlinear component
            // along x
            mat7_3n3  = mat_x.transpose() * mat5_n1n3;
            Bmat_nl_x.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
            
            // along y
            mat7_3n3  = mat_y.transpose() * mat5_n1n3;
            Bmat_nl_y.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
            
            // along z
            mat7_3n3  = mat_z.transpose() * mat5_n1n3;
            Bmat_nl_z.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
        }
    }
    
    
    // incompatible mode corrections
    K_alphaalpha = K_alphaalpha.inverse();
    K_corr = K_ualpha * K_alphaalpha * K_ualpha.transpose();
    
    if (request_jacobian)
        jac.topLeftCorner(n2, n2) -= K_corr;
    */
    
    ///////////////////////////////////////////////////////////////////////
    // second for loop to calculate the residual and stiffness contributions
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        //this->initialize_incompatible_strain_operator(qp, *_fe, Bmat_inc, Gmat);
        
        // calculate the stress
        stress = material_mat * (strain + Gmat * alpha);
        
        // residual from incompatible modes
        f_alpha += JxW[qp] * Gmat.transpose() * stress;
        
        // calculate contribution to the residual
        // linear strain operator
        Bmat_lin.vector_mult_transpose(vec2_n2, stress);
        f.topRows(n2) += JxW[qp] * vec2_n2;
        
        if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
            
            // nonlinear strain operator
            // x
            vec3_3 = mat_x.transpose() * stress;
            Bmat_nl_x.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
            
            // y
            vec3_3 = mat_y.transpose() * stress;
            Bmat_nl_y.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
            
            // z
            vec3_3 = mat_z.transpose() * stress;
            Bmat_nl_z.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
        }
        
        if (request_jacobian) {
            
            // the strain includes the following expansion
            // delta_epsilon = B_lin + mat_x B_x + mat_y B_y + mat_z B_z
            // Hence, the tangent stiffness matrix will include
            // components from epsilon^T C epsilon
            
            ////////////////////////////////////////////////////////
            // B_lin^T C B_lin
            Bmat_lin.left_multiply(mat1_n1n2, material_mat);
            Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            
            if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
                
                // B_x^T mat_x^T C B_lin
                mat3_3n2 = mat_x.transpose() * mat1_n1n2;
                Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // B_y^T mat_y^T C B_lin
                mat3_3n2 = mat_y.transpose() * mat1_n1n2;
                Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // B_z^T mat_z^T C B_lin
                mat3_3n2 = mat_z.transpose() * mat1_n1n2;
                Bmat_nl_z.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                ///////////////////////////////////////////////////////
                for (unsigned int i_dim=0; i_dim<3; i_dim++) {
                    switch (i_dim) {
                        case 0:
                            Bmat_nl_x.left_multiply(mat1_n1n2, mat_x);
                            break;
                            
                        case 1:
                            Bmat_nl_y.left_multiply(mat1_n1n2, mat_y);
                            break;
                            
                        case 2:
                            Bmat_nl_z.left_multiply(mat1_n1n2, mat_z);
                            break;
                    }
                    
                    // B_lin^T C mat_x_i B_x_i
                    mat1_n1n2 =  material_mat * mat1_n1n2;
                    Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_x^T mat_x^T C mat_x B_x
                    mat3_3n2 = mat_x.transpose() * mat1_n1n2;
                    Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_y^T mat_y^T C mat_x B_x
                    mat3_3n2 = mat_y.transpose() * mat1_n1n2;
                    Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_z^T mat_z^T C mat_x B_x
                    mat3_3n2 = mat_z.transpose() * mat1_n1n2;
                    Bmat_nl_z.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                }
                
                // use the stress to calculate the final contribution
                // to the Jacobian stiffness matrix
                mat4_33(0,0) = stress(0);
                mat4_33(1,1) = stress(1);
                mat4_33(2,2) = stress(2);
                mat4_33(0,1) = mat4_33(1,0) = stress(3);
                mat4_33(1,2) = mat4_33(2,1) = stress(4);
                mat4_33(0,2) = mat4_33(2,0) = stress(5);
                
                // u-disp
                Bmat_nl_u.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_u.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // v-disp
                Bmat_nl_v.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_v.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // w-disp
                Bmat_nl_w.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_w.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    // if jacobian is requested, add a small diagonal value for the
    // rotational dofs
    if (request_jacobian)
        jac.bottomRightCorner(n2, n2) += RealMatrixX::Identity(n2, n2) *
        1.0e-20 * jac.diagonal().maxCoeff();

    // correction to the residual from incompatible mode
    //f.topRows(n2) -= c * (K_alphaalpha * f_alpha);
    
    return request_jacobian;
}




void
MAST::StructuralElement3D::
update_incompatible_mode_solution(const RealVectorX& dsol) {
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW            = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = fe->get_xyz();
    const unsigned int
    n_phi              = (unsigned int)fe->n_shape_functions(),
    n1                 =6,
    n2                 =3*n_phi,
    n3                 =30,
    n_nodes            = _elem.get_reference_elem().n_nodes();
    
    RealMatrixX
    material_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    mat1_n1n2    = RealMatrixX::Zero(n1, n2),
    mat2_n2n2    = RealMatrixX::Zero(n2, n2),
    mat3_3n2     = RealMatrixX::Zero(3, n2),
    mat4_33      = RealMatrixX::Zero(3, 3),
    mat5_n1n3    = RealMatrixX::Zero(n1, n3),
    mat6_n2n3    = RealMatrixX::Zero(n2, n3),
    mat7_3n3     = RealMatrixX::Zero(3, n3),
    Gmat         = RealMatrixX::Zero(6, n3),
    K_alphaalpha = RealMatrixX::Zero(n3, n3),
    K_ualpha     = RealMatrixX::Zero(n2, n3),
    K_corr       = RealMatrixX::Zero(n2, n2);
    RealVectorX
    strain    = RealVectorX::Zero(6),
    stress    = RealVectorX::Zero(6),
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    vec3_3    = RealVectorX::Zero(3),
    local_disp= RealVectorX::Zero(n2),
    f         = RealVectorX::Zero(n3),
    alpha     = RealVectorX::Zero(n3);//*_incompatible_sol;

    // copy the values from the global to the local element
    local_disp.topRows(n2) = _local_sol.topRows(n2);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w,
    Bmat_inc;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, n_nodes);
    Bmat_nl_x.reinit(3, 3, n_nodes);
    Bmat_nl_y.reinit(3, 3, n_nodes);
    Bmat_nl_z.reinit(3, 3, n_nodes);
    Bmat_nl_u.reinit(3, 3, n_nodes);
    Bmat_nl_v.reinit(3, 3, n_nodes);
    Bmat_nl_w.reinit(3, 3, n_nodes);
    Bmat_inc.reinit(n1, n3, 1);            // six stress-strain components
    
    
    // initialize the incompatible mode mapping at element mid-point
    _init_incompatible_fe_mapping(_elem.get_reference_elem());

    ///////////////////////////////////////////////////////////////////
    // first for loop to evaluate alpha
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        //this->initialize_incompatible_strain_operator(qp, *_fe, Bmat_inc, Gmat);

        // calculate the stress
        stress = material_mat * (strain + Gmat * alpha);
        
        // residual of the incompatible strains
        f += JxW[qp] * Gmat.transpose() * stress;
        
        // calculate the incompatible mode matrices
        // incompatible mode diagonal stiffness matrix
        mat5_n1n3    =  material_mat * Gmat;
        K_alphaalpha += JxW[qp] * ( Gmat.transpose() * mat5_n1n3);

        // off-diagonal coupling matrix
        // linear strain term
        Bmat_lin.right_multiply_transpose(mat6_n2n3, mat5_n1n3);
        K_ualpha  += JxW[qp] * mat6_n2n3;
        
        if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
            
            // nonlinear component
            // along x
            mat7_3n3  = mat_x.transpose() * mat5_n1n3;
            Bmat_nl_x.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
            
            // along y
            mat7_3n3  = mat_y.transpose() * mat5_n1n3;
            Bmat_nl_y.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
            
            // along z
            mat7_3n3  = mat_z.transpose() * mat5_n1n3;
            Bmat_nl_z.right_multiply_transpose(mat6_n2n3, mat7_3n3);
            K_ualpha  += JxW[qp] * mat6_n2n3;
        }
    }
    
    
    // incompatible mode Jacobian inverse
    K_alphaalpha = K_alphaalpha.inverse();
    
    // update the alpha values
    alpha += K_alphaalpha * (-f - K_ualpha.transpose() * dsol.topRows(n2));
}



bool
MAST::StructuralElement3D::internal_residual_sensitivity(const MAST::FunctionBase& p,
                                                         bool request_jacobian,
                                                         RealVectorX& f,
                                                         RealMatrixX& jac) {
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW            = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = fe->get_xyz();
    const unsigned int
    n_phi              = (unsigned int)fe->n_shape_functions(),
    n1                 =6,
    n2                 =3*n_phi,
    n3                 =30,
    n_nodes            =_elem.get_reference_elem().n_nodes();
    
    RealMatrixX
    material_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    mat1_n1n2    = RealMatrixX::Zero(n1, n2),
    mat2_n2n2    = RealMatrixX::Zero(n2, n2),
    mat3_3n2     = RealMatrixX::Zero(3, n2),
    mat4_33      = RealMatrixX::Zero(3, 3),
    mat5_n1n3    = RealMatrixX::Zero(n1, n3),
    mat6_n2n3    = RealMatrixX::Zero(n2, n3),
    mat7_3n3     = RealMatrixX::Zero(3, n3),
    Gmat         = RealMatrixX::Zero(6, n3),
    K_alphaalpha = RealMatrixX::Zero(n3, n3),
    K_ualpha     = RealMatrixX::Zero(n2, n3),
    K_corr       = RealMatrixX::Zero(n2, n2);
    RealVectorX
    strain    = RealVectorX::Zero(6),
    stress    = RealVectorX::Zero(6),
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    vec3_3    = RealVectorX::Zero(3),
    local_disp= RealVectorX::Zero(n2),
    f_alpha   = RealVectorX::Zero(n3),
    alpha     = RealVectorX::Zero(n3);//*_incompatible_sol;
    
    // copy the values from the global to the local element
    local_disp.topRows(n2) = _local_sol.topRows(n2);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w,
    Bmat_inc;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, n_nodes);
    Bmat_nl_x.reinit(3, 3, n_nodes);
    Bmat_nl_y.reinit(3, 3, n_nodes);
    Bmat_nl_z.reinit(3, 3, n_nodes);
    Bmat_nl_u.reinit(3, 3, n_nodes);
    Bmat_nl_v.reinit(3, 3, n_nodes);
    Bmat_nl_w.reinit(3, 3, n_nodes);
    Bmat_inc.reinit(n1, n3, 1);            // six stress-strain components
    
    ///////////////////////////////////////////////////////////////////////
    // second for loop to calculate the residual and stiffness contributions
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        mat_stiff->derivative(p, xyz[qp], _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        //this->initialize_incompatible_strain_operator(qp, *_fe, Bmat_inc, Gmat);
        
        // calculate the stress
        stress = material_mat * (strain + Gmat * alpha);
        
        // residual from incompatible modes
        f_alpha += JxW[qp] * Gmat.transpose() * stress;
        
        // calculate contribution to the residual
        // linear strain operator
        Bmat_lin.vector_mult_transpose(vec2_n2, stress);
        f.topRows(n2) += JxW[qp] * vec2_n2;
        
        if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
            
            // nonlinear strain operator
            // x
            vec3_3 = mat_x.transpose() * stress;
            Bmat_nl_x.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
            
            // y
            vec3_3 = mat_y.transpose() * stress;
            Bmat_nl_y.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
            
            // z
            vec3_3 = mat_z.transpose() * stress;
            Bmat_nl_z.vector_mult_transpose(vec2_n2, vec3_3);
            f.topRows(n2) += JxW[qp] * vec2_n2;
        }
        
        if (request_jacobian) {
            
            // the strain includes the following expansion
            // delta_epsilon = B_lin + mat_x B_x + mat_y B_y + mat_z B_z
            // Hence, the tangent stiffness matrix will include
            // components from epsilon^T C epsilon
            
            ////////////////////////////////////////////////////////
            // B_lin^T C B_lin
            Bmat_lin.left_multiply(mat1_n1n2, material_mat);
            Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            
            if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
                
                // B_x^T mat_x^T C B_lin
                mat3_3n2 = mat_x.transpose() * mat1_n1n2;
                Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // B_y^T mat_y^T C B_lin
                mat3_3n2 = mat_y.transpose() * mat1_n1n2;
                Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // B_z^T mat_z^T C B_lin
                mat3_3n2 = mat_z.transpose() * mat1_n1n2;
                Bmat_nl_z.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                ///////////////////////////////////////////////////////
                for (unsigned int i_dim=0; i_dim<3; i_dim++) {
                    switch (i_dim) {
                        case 0:
                            Bmat_nl_x.left_multiply(mat1_n1n2, mat_x);
                            break;
                            
                        case 1:
                            Bmat_nl_y.left_multiply(mat1_n1n2, mat_y);
                            break;
                            
                        case 2:
                            Bmat_nl_z.left_multiply(mat1_n1n2, mat_z);
                            break;
                    }
                    
                    // B_lin^T C mat_x_i B_x_i
                    mat1_n1n2 =  material_mat * mat1_n1n2;
                    Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_x^T mat_x^T C mat_x B_x
                    mat3_3n2 = mat_x.transpose() * mat1_n1n2;
                    Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_y^T mat_y^T C mat_x B_x
                    mat3_3n2 = mat_y.transpose() * mat1_n1n2;
                    Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                    // B_z^T mat_z^T C mat_x B_x
                    mat3_3n2 = mat_z.transpose() * mat1_n1n2;
                    Bmat_nl_z.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                    jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                    
                }
                
                // use the stress to calculate the final contribution
                // to the Jacobian stiffness matrix
                mat4_33(0,0) = stress(0);
                mat4_33(1,1) = stress(1);
                mat4_33(2,2) = stress(2);
                mat4_33(0,1) = mat4_33(1,0) = stress(3);
                mat4_33(1,2) = mat4_33(2,1) = stress(4);
                mat4_33(0,2) = mat4_33(2,0) = stress(5);
                
                // u-disp
                Bmat_nl_u.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_u.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // v-disp
                Bmat_nl_v.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_v.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // w-disp
                Bmat_nl_w.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_w.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    // if jacobian is requested, add a small diagonal value for the
    // rotational dofs
    if (request_jacobian)
        jac.bottomRightCorner(n2, n2) += RealMatrixX::Identity(n2, n2) *
        1.0e-20 * jac.diagonal().maxCoeff();
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::prestress_residual (bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac) {
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::prestress_residual_sensitivity (const MAST::FunctionBase& p,
                                                           bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac) {
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::
surface_pressure_residual(bool request_jacobian,
                          RealVectorX &f,
                          RealMatrixX &jac,
                          const unsigned int side,
                          MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_side_fe(side, false, false));

    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_reference_coordinate();
    const unsigned int
    n_phi  = (unsigned int)phi.size(),
    n1     = 3,
    n2     = 3*n_phi;
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    
    
    FEMOperatorMatrix Bmat;
    Real press;
    
    RealVectorX
    phi_vec     = RealVectorX::Zero(n_phi),
    force       = RealVectorX::Zero(n1),
    local_f     = RealVectorX::Zero(n2),
    vec_n2      = RealVectorX::Zero(n2);
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(n1, phi_vec);
        
        // get pressure value
        func(qpoint[qp], _time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    f.topRows(n2) -= local_f;
    
    return (request_jacobian);
}





bool
MAST::StructuralElement3D::
surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                      bool request_jacobian,
                                      RealVectorX &f,
                                      RealMatrixX &jac,
                                      const unsigned int side,
                                      MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_side_fe(side, false, false));

    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_reference_coordinate();
    const unsigned int
    n_phi  = (unsigned int)phi.size(),
    n1     = 3,
    n2     = 6*n_phi;
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    
    
    FEMOperatorMatrix Bmat;
    Real press;
    
    RealVectorX
    phi_vec     = RealVectorX::Zero(n_phi),
    force       = RealVectorX::Zero(2*n1),
    local_f     = RealVectorX::Zero(n2),
    vec_n2      = RealVectorX::Zero(n2);
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func.derivative(p, qpoint[qp], _time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    f -= local_f;
    
    return (request_jacobian);
}





bool
MAST::StructuralElement3D::thermal_residual(bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac,
                                            MAST::BoundaryConditionBase& bc) {
    
    std::unique_ptr<MAST::FEBase>
    fe(_elem.init_fe(true, false));

    const std::vector<Real>& JxW            = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = fe->get_xyz();
    
    const unsigned int
    n_phi   = (unsigned int)fe->get_phi().size(),
    n1      = 6,
    n2      = 3*n_phi,
    n_nodes = _elem.get_reference_elem().n_nodes();
    
    RealMatrixX
    material_exp_A_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    mat2_n2n2    = RealMatrixX::Zero(n2,n2),
    mat3_3n2     = RealMatrixX::Zero(3,n2),
    mat4_33      = RealMatrixX::Zero(3,3);
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_3    = RealVectorX::Zero(3),
    vec3_n2   = RealVectorX::Zero(n2),
    delta_t   = RealVectorX::Zero(1),
    local_disp= RealVectorX::Zero(n2),
    strain    = RealVectorX::Zero(6);

    // copy the values from the global to the local element
    local_disp.topRows(n2) = _local_sol.topRows(n2);
    
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, n_nodes);
    Bmat_nl_x.reinit(3, 3, n_nodes);
    Bmat_nl_y.reinit(3, 3, n_nodes);
    Bmat_nl_z.reinit(3, 3, n_nodes);
    Bmat_nl_u.reinit(3, 3, n_nodes);
    Bmat_nl_v.reinit(3, 3, n_nodes);
    Bmat_nl_w.reinit(3, 3, n_nodes);

    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    mat = _property.thermal_expansion_A_matrix(*this);
    
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real t, t0;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*mat)       (xyz[qp], _time, material_exp_A_mat);
        temp_func    (xyz[qp], _time, t);
        ref_temp_func(xyz[qp], _time, t0);
        delta_t(0) = t-t0;
        
        vec1_n1 = material_exp_A_mat * delta_t; // [C]{alpha (T - T0)}
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        
        // linear strain operotor
        Bmat_lin.vector_mult_transpose(vec3_n2, vec1_n1);
        f.topRows(n2) -= JxW[qp] * vec3_n2;
        
        if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
            
            // nonlinear strain operotor
            // x
            vec2_3 = mat_x.transpose() * vec1_n1;
            Bmat_nl_x.vector_mult_transpose(vec3_n2, vec2_3);
            f.topRows(n2) -= JxW[qp] * vec3_n2;
            
            // y
            vec2_3 = mat_y.transpose() * vec1_n1;
            Bmat_nl_y.vector_mult_transpose(vec3_n2, vec2_3);
            f.topRows(n2) -= JxW[qp] * vec3_n2;
            
            // z
            vec2_3 = mat_z.transpose() * vec1_n1;
            Bmat_nl_z.vector_mult_transpose(vec3_n2, vec2_3);
            f.topRows(n2) -= JxW[qp] * vec3_n2;
            
            // Jacobian for the nonlinear case
            if (request_jacobian) {
                
                // use the stress to calculate the final contribution
                // to the Jacobian stiffness matrix
                mat4_33(0,0) = vec1_n1(0);
                mat4_33(1,1) = vec1_n1(1);
                mat4_33(2,2) = vec1_n1(2);
                mat4_33(0,1) = mat4_33(1,0) = vec1_n1(3);
                mat4_33(1,2) = mat4_33(2,1) = vec1_n1(4);
                mat4_33(0,2) = mat4_33(2,0) = vec1_n1(5);
                
                // u-disp
                Bmat_nl_u.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_u.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) -= JxW[qp] * mat2_n2n2;
                
                // v-disp
                Bmat_nl_v.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_v.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) -= JxW[qp] * mat2_n2n2;
                
                // w-disp
                Bmat_nl_w.left_multiply(mat3_3n2, mat4_33);
                Bmat_nl_w.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) -= JxW[qp] * mat2_n2n2;
            }
        }
    }

    // Jacobian contribution from von Karman strain
    return request_jacobian;
}



bool
MAST::StructuralElement3D::thermal_residual_sensitivity(const MAST::FunctionBase& p,
                                                        bool request_jacobian,
                                                        RealVectorX& f,
                                                        RealMatrixX& jac,
                                                        MAST::BoundaryConditionBase& bc) {
    
    // to be implemented
    libmesh_error();
    
    return false;
}




bool
MAST::StructuralElement3D::
piston_theory_residual(bool request_jacobian,
                       RealVectorX &f,
                       RealMatrixX& jac_xdot,
                       RealMatrixX& jac,
                       const unsigned int side,
                       MAST::BoundaryConditionBase& bc) {
    
    
    libmesh_error(); // to be implemented
    
    return (request_jacobian);
}



bool
MAST::StructuralElement3D::
piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                   bool request_jacobian,
                                   RealVectorX &f,
                                   RealMatrixX& jac_xdot,
                                   RealMatrixX& jac,
                                   const unsigned int side,
                                   MAST::BoundaryConditionBase& bc) {
    
    
    libmesh_error(); // to be implemented
    
    return (request_jacobian);
}





bool
MAST::StructuralElement3D::calculate_stress(bool request_derivative,
                                            const MAST::FunctionBase* p,
                                            MAST::StressStrainOutputBase& output) {
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true, false));
    std::vector<libMesh::Point>     qp_loc = fe->get_qpoints();


    // now that the FE object has been initialized, evaluate the stress values
    
    
    const std::vector<Real> &JxW              = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz    = fe->get_xyz();
    const unsigned int
    n_phi              = (unsigned int)fe->n_shape_functions(),
    n1                 =6,
    n2                 =3*n_phi,
    n3                 =30,
    n_nodes            = _elem.get_reference_elem().n_nodes() ;
    
    RealMatrixX
    material_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    Gmat         = RealMatrixX::Zero(6, n3);
    
    RealVectorX
    strain    = RealVectorX::Zero(6),
    stress    = RealVectorX::Zero(6),
    local_disp= RealVectorX::Zero(n2),
    alpha     = RealVectorX::Zero(n3);//*_incompatible_sol;

    // copy the values from the global to the local element
    local_disp.topRows(n2) = _local_sol.topRows(n2);
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w,
    Bmat_inc;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, n_nodes);
    Bmat_nl_x.reinit(3, 3, n_nodes);
    Bmat_nl_y.reinit(3, 3, n_nodes);
    Bmat_nl_z.reinit(3, 3, n_nodes);
    Bmat_nl_u.reinit(3, 3, n_nodes);
    Bmat_nl_v.reinit(3, 3, n_nodes);
    Bmat_nl_w.reinit(3, 3, n_nodes);
    Bmat_inc.reinit(n1, n3, 1);            // six stress-strain components
    
    // a reference to the stress output data structure
    MAST::StressStrainOutputBase& stress_output =
    dynamic_cast<MAST::StressStrainOutputBase&>(output);
    
    // initialize the incompatible mode mapping at element mid-point
    _init_incompatible_fe_mapping(_elem.get_reference_elem());
    
    ///////////////////////////////////////////////////////////////////////
    // second for loop to calculate the residual and stiffness contributions
    for (unsigned int qp=0; qp<qp_loc.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff)(xyz[qp], _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
                                                        *fe,
                                                        local_disp,
                                                        strain,
                                                        mat_x, mat_y, mat_z,
                                                        Bmat_lin,
                                                        Bmat_nl_x,
                                                        Bmat_nl_y,
                                                        Bmat_nl_z,
                                                        Bmat_nl_u,
                                                        Bmat_nl_v,
                                                        Bmat_nl_w);
        //this->initialize_incompatible_strain_operator(qp, *fe, Bmat_inc, Gmat);
        
        // calculate the stress
        strain += Gmat * alpha;
        stress = material_mat * strain;
        
        // set the stress and strain data
        MAST::StressStrainOutputBase::Data*
        data = nullptr;
        
        // if neither the derivative nor sensitivity is requested, then
        // we assume that a new data entry is to be provided. Otherwise,
        // we assume that the stress at this quantity already
        // exists, and we only need to append sensitivity/derivative
        // data to it
        if (!request_derivative && !p)
            data = &(stress_output.add_stress_strain_at_qp_location(_elem,
                                                                    qp,
                                                                    qp_loc[qp],
                                                                    xyz[qp],
                                                                    stress,
                                                                    strain,
                                                                    JxW[qp]));
        else
            data = &(stress_output.get_stress_strain_data_for_elem_at_qp(_elem, qp));

        
        if (request_derivative) {
            // to be implemented
            libmesh_error();
        }
    }
    
    return request_derivative;
}






void
MAST::StructuralElement3D::initialize_strain_operator (const unsigned int qp,
                                                       const MAST::FEBase& fe,
                                                       FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >&
    dphi = fe.get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi  = RealVectorX::Zero(n_phi);
    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(3, 1, phi); //  gamma_xy = dv/dx + ...
    Bmat.set_shape_function(5, 2, phi); //  gamma_zx = dw/dx + ...
    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    Bmat.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat.set_shape_function(3, 0, phi); //  gamma_xy = du/dy + ...
    Bmat.set_shape_function(4, 2, phi); //  gamma_yz = dw/dy + ...
    
    // dN/dz
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](2);
    Bmat.set_shape_function(2, 2, phi); //  epsilon_xx = dw/dz
    Bmat.set_shape_function(4, 1, phi); //  gamma_xy = dv/dz + ...
    Bmat.set_shape_function(5, 0, phi); //  gamma_zx = du/dz + ...
}



void
MAST::StructuralElement3D::
initialize_green_lagrange_strain_operator(const unsigned int qp,
                                          const MAST::FEBase& fe,
                                          const RealVectorX& local_disp,
                                          RealVectorX& epsilon,
                                          RealMatrixX& mat_x,
                                          RealMatrixX& mat_y,
                                          RealMatrixX& mat_z,
                                          MAST::FEMOperatorMatrix& Bmat_lin,
                                          MAST::FEMOperatorMatrix& Bmat_nl_x,
                                          MAST::FEMOperatorMatrix& Bmat_nl_y,
                                          MAST::FEMOperatorMatrix& Bmat_nl_z,
                                          MAST::FEMOperatorMatrix& Bmat_nl_u,
                                          MAST::FEMOperatorMatrix& Bmat_nl_v,
                                          MAST::FEMOperatorMatrix& Bmat_nl_w) {
    
    epsilon.setZero();
    mat_x.setZero();
    mat_y.setZero();
    mat_z.setZero();
    
    const std::vector<std::vector<libMesh::RealVectorValue> >&
    dphi = fe.get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi  = RealVectorX::Zero(n_phi);

    // make sure all matrices are the right size
    libmesh_assert_equal_to(epsilon.size(), 6);
    libmesh_assert_equal_to(mat_x.rows(), 6);
    libmesh_assert_equal_to(mat_x.cols(), 3);
    libmesh_assert_equal_to(mat_y.rows(), 6);
    libmesh_assert_equal_to(mat_y.cols(), 3);
    libmesh_assert_equal_to(mat_z.rows(), 6);
    libmesh_assert_equal_to(mat_z.cols(), 3);
    libmesh_assert_equal_to(Bmat_lin.m(), 6);
    libmesh_assert_equal_to(Bmat_lin.n(), 3*n_phi);
    libmesh_assert_equal_to(Bmat_nl_x.m(), 3);
    libmesh_assert_equal_to(Bmat_nl_x.n(), 3*n_phi);
    libmesh_assert_equal_to(Bmat_nl_y.m(), 3);
    libmesh_assert_equal_to(Bmat_nl_y.n(), 3*n_phi);
    libmesh_assert_equal_to(Bmat_nl_z.m(), 3);
    libmesh_assert_equal_to(Bmat_nl_z.n(), 3*n_phi);

    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    // linear strain operator
    Bmat_lin.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat_lin.set_shape_function(3, 1, phi); //  gamma_xy = dv/dx + ...
    Bmat_lin.set_shape_function(5, 2, phi); //  gamma_zx = dw/dx + ...
    
    if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
        
        // nonlinear strain operator in x
        Bmat_nl_x.set_shape_function(0, 0, phi); // du/dx
        Bmat_nl_x.set_shape_function(1, 1, phi); // dv/dx
        Bmat_nl_x.set_shape_function(2, 2, phi); // dw/dx
        
        // nonlinear strain operator in u
        Bmat_nl_u.set_shape_function(0, 0, phi); // du/dx
        Bmat_nl_v.set_shape_function(0, 1, phi); // dv/dx
        Bmat_nl_w.set_shape_function(0, 2, phi); // dw/dx
    }
    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    // linear strain operator
    Bmat_lin.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat_lin.set_shape_function(3, 0, phi); //  gamma_xy = du/dy + ...
    Bmat_lin.set_shape_function(4, 2, phi); //  gamma_yz = dw/dy + ...
    
    if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
        
        // nonlinear strain operator in y
        Bmat_nl_y.set_shape_function(0, 0, phi); // du/dy
        Bmat_nl_y.set_shape_function(1, 1, phi); // dv/dy
        Bmat_nl_y.set_shape_function(2, 2, phi); // dw/dy
        
        // nonlinear strain operator in v
        Bmat_nl_u.set_shape_function(1, 0, phi); // du/dy
        Bmat_nl_v.set_shape_function(1, 1, phi); // dv/dy
        Bmat_nl_w.set_shape_function(1, 2, phi); // dw/dy
    }
    
    // dN/dz
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](2);
    Bmat_lin.set_shape_function(2, 2, phi); //  epsilon_xx = dw/dz
    Bmat_lin.set_shape_function(4, 1, phi); //  gamma_xy = dv/dz + ...
    Bmat_lin.set_shape_function(5, 0, phi); //  gamma_zx = du/dz + ...
    
    if (_property.strain_type() == MAST::NONLINEAR_STRAIN) {
        
        // nonlinear strain operator in z
        Bmat_nl_z.set_shape_function(0, 0, phi); // du/dz
        Bmat_nl_z.set_shape_function(1, 1, phi); // dv/dz
        Bmat_nl_z.set_shape_function(2, 2, phi); // dw/dz
        
        // nonlinear strain operator in w
        Bmat_nl_u.set_shape_function(2, 0, phi); // du/dz
        Bmat_nl_v.set_shape_function(2, 1, phi); // dv/dz
        Bmat_nl_w.set_shape_function(2, 2, phi); // dw/dz
        
        
        // calculate the displacement gradient to create the GL strain
        RealVectorX
        ddisp_dx = RealVectorX::Zero(3),
        ddisp_dy = RealVectorX::Zero(3),
        ddisp_dz = RealVectorX::Zero(3);
        
        Bmat_nl_x.vector_mult(ddisp_dx, local_disp);  // {du/dx, dv/dx, dw/dx}
        Bmat_nl_y.vector_mult(ddisp_dy, local_disp);  // {du/dy, dv/dy, dw/dy}
        Bmat_nl_z.vector_mult(ddisp_dz, local_disp);  // {du/dz, dv/dz, dw/dz}
        
        // prepare the displacement gradient matrix: F = grad(u)
        RealMatrixX
        F = RealMatrixX::Zero(3,3),
        E = RealMatrixX::Zero(3,3);
        F.col(0) = ddisp_dx;
        F.col(1) = ddisp_dy;
        F.col(2) = ddisp_dz;
        
        // this calculates the Green-Lagrange strain in the reference config
        E = 0.5*(F + F.transpose() + F.transpose() * F);
        
        // now, add this to the strain vector
        epsilon(0) = E(0,0);
        epsilon(1) = E(1,1);
        epsilon(2) = E(2,2);
        epsilon(3) = E(0,1) + E(1,0);
        epsilon(4) = E(1,2) + E(2,1);
        epsilon(5) = E(0,2) + E(2,0);
        
        // now initialize the matrices with strain components
        // that multiply the Bmat_nl terms
        mat_x.row(0) =     ddisp_dx;
        mat_x.row(3) =     ddisp_dy;
        mat_x.row(5) =     ddisp_dz;
        
        mat_y.row(1) =     ddisp_dy;
        mat_y.row(3) =     ddisp_dx;
        mat_y.row(4) =     ddisp_dz;
        
        
        mat_z.row(2) =     ddisp_dz;
        mat_z.row(4) =     ddisp_dy;
        mat_z.row(5) =     ddisp_dx;
    }
    else
        Bmat_lin.vector_mult(epsilon, local_disp);
}




void
MAST::StructuralElement3D::
initialize_incompatible_strain_operator(const unsigned int  qp,
                                        const MAST::FEBase& fe,
                                        FEMOperatorMatrix&  Bmat,
                                        RealMatrixX&        G_mat) {
    
    RealVectorX phi_vec = RealVectorX::Zero(1);
    
    // get the location of element coordinates
    const std::vector<libMesh::Point>& q_point = fe.get_qrule().get_points();
    const Real
    xi  = q_point[qp](0),
    eta = q_point[qp](1),
    phi = q_point[qp](2);
    
    const std::vector<std::vector<Real> >&
    dshapedxi  = fe.get_dphidxi(),
    dshapedeta = fe.get_dphideta(),
    dshapedphi = fe.get_dphidzeta();
    
    const libMesh::Elem& e = _elem.get_reference_elem();
    
    // calculate the deformed xyz coordinates
    const unsigned int n_nodes = e.n_nodes();
    RealVectorX
    xdef     = RealVectorX::Zero(n_nodes),
    ydef     = RealVectorX::Zero(n_nodes),
    zdef     = RealVectorX::Zero(n_nodes),
    phivec   = RealVectorX::Zero(n_nodes);
    
    // set the current values of nodal coordinates
    for (unsigned int i_node=0; i_node<n_nodes; i_node++) {
        xdef(i_node) = e.point(i_node)(0);// + _sol(i_node*3+0);
        ydef(i_node) = e.point(i_node)(1);// + _sol(i_node*3+1);
        zdef(i_node) = e.point(i_node)(2);// + _sol(i_node*3+2);
    }
    
    // calculate dxyz/dxi
    // make sure that the number of shape functions is the same as the number
    // of nodes. Meaning that this formulation is limited to Lagrange
    // elemnts only.
    libmesh_assert_equal_to(dshapedxi.size(), n_nodes);

    RealMatrixX
    jac = RealMatrixX::Zero(3,3);
    
    // first derivatives wrt xi
    for (unsigned int i_node=0; i_node<n_nodes; i_node++)
        phivec(i_node)  =  dshapedxi[i_node][qp];
    
    jac(0,0) =  phivec.dot(xdef);
    jac(0,1) =  phivec.dot(ydef);
    jac(0,2) =  phivec.dot(zdef);
    
    
    // second, derivatives wrt eta
    for (unsigned int i_node=0; i_node<n_nodes; i_node++)
        phivec(i_node)  =  dshapedeta[i_node][qp];
    
    jac(1,0) =  phivec.dot(xdef);
    jac(1,1) =  phivec.dot(ydef);
    jac(1,2) =  phivec.dot(zdef);
    
    // lastly, derivatives wrt phi
    for (unsigned int i_node=0; i_node<n_nodes; i_node++)
        phivec(i_node)  =  dshapedphi[i_node][qp];
    
    jac(2,0) =  phivec.dot(xdef);
    jac(2,1) =  phivec.dot(ydef);
    jac(2,2) =  phivec.dot(zdef);

    
    // now set the shape function values
    // epsilon_xx
    phi_vec(0) =         xi;   Bmat.set_shape_function(0,  0, phi_vec);
    phi_vec(0) =     xi*eta;   Bmat.set_shape_function(0, 15, phi_vec);
    phi_vec(0) =     xi*phi;   Bmat.set_shape_function(0, 16, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(0, 24, phi_vec);

    
    // epsilon_yy
    phi_vec(0) =        eta;   Bmat.set_shape_function(1,  1, phi_vec);
    phi_vec(0) =     xi*eta;   Bmat.set_shape_function(1, 17, phi_vec);
    phi_vec(0) =    eta*phi;   Bmat.set_shape_function(1, 18, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(1, 25, phi_vec);

    // epsilon_zz
    phi_vec(0) =        phi;   Bmat.set_shape_function(2,  2, phi_vec);
    phi_vec(0) =     xi*phi;   Bmat.set_shape_function(2, 19, phi_vec);
    phi_vec(0) =    eta*phi;   Bmat.set_shape_function(2, 20, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(2, 26, phi_vec);

    // epsilon_xy
    phi_vec(0) =         xi;   Bmat.set_shape_function(3,  3, phi_vec);
    phi_vec(0) =        eta;   Bmat.set_shape_function(3,  4, phi_vec);
    phi_vec(0) =     xi*phi;   Bmat.set_shape_function(3,  9, phi_vec);
    phi_vec(0) =    eta*phi;   Bmat.set_shape_function(3, 10, phi_vec);
    phi_vec(0) =     xi*eta;   Bmat.set_shape_function(3, 21, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(3, 27, phi_vec);

    // epsilon_yz
    phi_vec(0) =        eta;   Bmat.set_shape_function(4,  5, phi_vec);
    phi_vec(0) =        phi;   Bmat.set_shape_function(4,  6, phi_vec);
    phi_vec(0) =     xi*eta;   Bmat.set_shape_function(4, 11, phi_vec);
    phi_vec(0) =     xi*phi;   Bmat.set_shape_function(4, 12, phi_vec);
    phi_vec(0) =    eta*phi;   Bmat.set_shape_function(4, 22, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(4, 28, phi_vec);

    // epsilon_xz
    phi_vec(0) =         xi;   Bmat.set_shape_function(5,  7, phi_vec);
    phi_vec(0) =        phi;   Bmat.set_shape_function(5,  8, phi_vec);
    phi_vec(0) =     xi*eta;   Bmat.set_shape_function(5, 13, phi_vec);
    phi_vec(0) =    eta*phi;   Bmat.set_shape_function(5, 14, phi_vec);
    phi_vec(0) =     xi*phi;   Bmat.set_shape_function(5, 23, phi_vec);
    phi_vec(0) = xi*eta*phi;   Bmat.set_shape_function(5, 29, phi_vec);
    
    
    Bmat.left_multiply(G_mat, _T0_inv_tr);
    G_mat /= jac.determinant();
}




void
MAST::StructuralElement3D::_init_incompatible_fe_mapping( const libMesh::Elem& e) {
    
    libmesh_assert(e.type() == libMesh::HEX8);
    
    unsigned int nv = _system.system().n_vars();
    
    libmesh_assert (nv);
    libMesh::FEType fe_type = _system.system().variable_type(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.system().variable_type(i));
    
    // Create an adequate quadrature rule
    std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(e.dim(), fe_type).release());
    const std::vector<libMesh::Point> pts(1); // initializes point to (0,0,0)
    
    fe->get_dphidxi();
    fe->get_dphideta();
    fe->get_dphidzeta();
    
    fe->reinit(&e, &pts); // reinit at (0,0,0)
    
    _T0_inv_tr = RealMatrixX::Zero(6,6);
    
    const std::vector<std::vector<Real> >&
    dshapedxi  = fe->get_dphidxi(),
    dshapedeta = fe->get_dphideta(),
    dshapedphi = fe->get_dphidzeta();

    // calculate the deformed xyz coordinates
    RealVectorX
    xdef     = RealVectorX::Zero(e.n_nodes()),
    ydef     = RealVectorX::Zero(e.n_nodes()),
    zdef     = RealVectorX::Zero(e.n_nodes()),
    phi      = RealVectorX::Zero(e.n_nodes());
    
    // set the current values of nodal coordinates
    for (unsigned int i_node=0; i_node<e.n_nodes(); i_node++) {
        xdef(i_node) = e.point(i_node)(0);// + _local_sol(i_node*3+0);
        ydef(i_node) = e.point(i_node)(1);// + _local_sol(i_node*3+1);
        zdef(i_node) = e.point(i_node)(2);// + _local_sol(i_node*3+2);
    }
    
    // calculate dxyz/dxi
    // make sure that the number of shape functions is the same as the number
    // of nodes. Meaning that this formulation is limited to Lagrange
    // elemnts only.
    libmesh_assert_equal_to(dshapedxi.size(), e.n_nodes());
    
    RealMatrixX
    jac = RealMatrixX::Zero(3,3),
    T0  = RealMatrixX::Zero(6,6);

    // first derivatives wrt xi
    for (unsigned int i_node=0; i_node<e.n_nodes(); i_node++)
        phi(i_node)  =  dshapedxi[i_node][0];

    jac(0,0) =  phi.dot(xdef);
    jac(0,1) =  phi.dot(ydef);
    jac(0,2) =  phi.dot(zdef);
    

    // second, derivatives wrt eta
    for (unsigned int i_node=0; i_node<e.n_nodes(); i_node++)
        phi(i_node)  =  dshapedeta[i_node][0];

    jac(1,0) =  phi.dot(xdef);
    jac(1,1) =  phi.dot(ydef);
    jac(1,2) =  phi.dot(zdef);
    
    // lastly, derivatives wrt phi
    for (unsigned int i_node=0; i_node<e.n_nodes(); i_node++)
        phi(i_node)  =  dshapedphi[i_node][0];

    jac(2,0) =  phi.dot(xdef);
    jac(2,1) =  phi.dot(ydef);
    jac(2,2) =  phi.dot(zdef);
    
    
    // we first set the values of the T0 matrix and then get its inverse
    T0(0,0) =   jac(0,0)*jac(0,0);
    T0(0,1) =   jac(1,0)*jac(1,0);
    T0(0,2) =   jac(2,0)*jac(2,0);
    T0(0,3) = 2*jac(0,0)*jac(1,0);
    T0(0,4) = 2*jac(1,0)*jac(2,0);
    T0(0,5) = 2*jac(2,0)*jac(0,0);
    
    T0(1,0) =   jac(0,1)*jac(0,1);
    T0(1,1) =   jac(1,1)*jac(1,1);
    T0(1,2) =   jac(2,1)*jac(2,1);
    T0(1,3) = 2*jac(0,1)*jac(1,1);
    T0(1,4) = 2*jac(1,1)*jac(2,1);
    T0(1,5) = 2*jac(2,1)*jac(0,1);

    T0(2,0) =   jac(0,2)*jac(0,2);
    T0(2,1) =   jac(1,2)*jac(1,2);
    T0(2,2) =   jac(2,2)*jac(2,2);
    T0(2,3) = 2*jac(0,2)*jac(1,2);
    T0(2,4) = 2*jac(1,2)*jac(2,2);
    T0(2,5) = 2*jac(2,2)*jac(0,2);

    T0(3,0) =   jac(0,0)*jac(0,1);
    T0(3,1) =   jac(1,0)*jac(1,1);
    T0(3,2) =   jac(2,0)*jac(2,1);
    T0(3,3) =   jac(0,0)*jac(1,1)+jac(1,0)*jac(0,1);
    T0(3,4) =   jac(1,0)*jac(2,1)+jac(2,0)*jac(1,1);
    T0(3,5) =   jac(2,1)*jac(0,0)+jac(2,0)*jac(0,1);

    T0(4,0) =   jac(0,1)*jac(0,2);
    T0(4,1) =   jac(1,1)*jac(1,2);
    T0(4,2) =   jac(2,1)*jac(2,2);
    T0(4,3) =   jac(0,1)*jac(1,2)+jac(1,1)*jac(0,2);
    T0(4,4) =   jac(1,1)*jac(2,2)+jac(2,1)*jac(1,2);
    T0(4,5) =   jac(2,2)*jac(0,1)+jac(2,1)*jac(0,2);

    T0(5,0) =   jac(0,0)*jac(0,2);
    T0(5,1) =   jac(1,0)*jac(1,2);
    T0(5,2) =   jac(2,0)*jac(2,2);
    T0(5,3) =   jac(0,0)*jac(1,2)+jac(1,0)*jac(0,2);
    T0(5,4) =   jac(1,0)*jac(2,2)+jac(2,0)*jac(1,2);
    T0(5,5) =   jac(2,2)*jac(0,0)+jac(2,0)*jac(0,2);
    
    _T0_inv_tr = jac.determinant() * T0.inverse().transpose();
}


