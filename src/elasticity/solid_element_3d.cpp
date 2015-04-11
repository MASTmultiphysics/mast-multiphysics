
// MAST includes
#include "elasticity/solid_element_3d.h"
#include "numerics/fem_operator_matrix.h"
#include "mesh/local_elem_base.h"
#include "property_cards/element_property_card_base.h"
#include "base/boundary_condition_base.h"



bool
MAST::StructuralElement3D::internal_residual(bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac,
                                             bool if_ignore_ho_jac) {
    
    const std::vector<Real>& JxW            = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = _fe->get_xyz();
    const unsigned int
    n_phi              = (unsigned int)_fe->n_shape_functions(),
    n1                 =6,
    n2                 =3*n_phi;
    
    RealMatrixX
    material_mat,
    mat_x        = RealMatrixX::Zero(6,3),
    mat_y        = RealMatrixX::Zero(6,3),
    mat_z        = RealMatrixX::Zero(6,3),
    mat1_n1n2    = RealMatrixX::Zero(n1, n2),
    mat2_n2n2    = RealMatrixX::Zero(n2, n2),
    mat3_3n2     = RealMatrixX::Zero(3, n2),
    mat4_33      = RealMatrixX::Zero(3, 3);
    RealVectorX
    strain    = RealVectorX::Zero(6),
    stress    = RealVectorX::Zero(6),
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    vec3_3    = RealVectorX::Zero(3),
    local_disp= RealVectorX::Zero(n2);
    
    // copy the values from the global to the local element
    for (unsigned int i=0; i<n2; i++) local_disp(i) = _local_sol(i);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    libMesh::Point p;
    MAST::FEMOperatorMatrix
    Bmat_lin,
    Bmat_nl_x,
    Bmat_nl_y,
    Bmat_nl_z,
    Bmat_nl_u,
    Bmat_nl_v,
    Bmat_nl_w;
    // six stress components, related to three displacements
    Bmat_lin.reinit(n1, 3, _elem.n_nodes());
    Bmat_nl_x.reinit(3, 3, _elem.n_nodes());
    Bmat_nl_y.reinit(3, 3, _elem.n_nodes());
    Bmat_nl_z.reinit(3, 3, _elem.n_nodes());
    Bmat_nl_u.reinit(3, 3, _elem.n_nodes());
    Bmat_nl_v.reinit(3, 3, _elem.n_nodes());
    Bmat_nl_w.reinit(3, 3, _elem.n_nodes());
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff)(p, _time, material_mat);
        
        this->initialize_green_lagrange_strain_operator(qp,
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
        
        // calculate the stress
        stress = material_mat * strain;
        
        // calculate contribution to the residual
        // linear strain operator
        Bmat_lin.vector_mult_transpose(vec2_n2, stress);
        f.topRows(n2) += JxW[qp] * vec2_n2;
        
        // nonlinear strain operator
        // x
        vec3_3 = 0.5 * mat_x.transpose() * stress;
        Bmat_nl_x.vector_mult_transpose(vec2_n2, vec3_3);
        f.topRows(n2) += JxW[qp] * vec2_n2;

        // y
        vec3_3 = 0.5 * mat_y.transpose() * stress;
        Bmat_nl_y.vector_mult_transpose(vec2_n2, vec3_3);
        f.topRows(n2) += JxW[qp] * vec2_n2;

        // z
        vec3_3 = 0.5 * mat_z.transpose() * stress;
        Bmat_nl_z.vector_mult_transpose(vec2_n2, vec3_3);
        f.topRows(n2) += JxW[qp] * vec2_n2;

        
        if (request_jacobian) {
            
            // the strain includes the following expansion
            // epsilon = (B_lin + .5 (mat_x B_x + mat_y B_y + mat_z B_z ) )
            // Hence, the tangent stiffness matrix will include
            // components from epsilon^T C epsilon
            
            ////////////////////////////////////////////////////////
            // B_lin^T C B_lin
            Bmat_lin.left_multiply(mat1_n1n2, material_mat);
            Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            
            // 1/2 B_x^T mat_x^T C B_lin
            mat3_3n2 = 0.5 * mat_x.transpose() * mat1_n1n2;
            Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
            jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            
            // 1/2 B_y^T mat_y^T C B_lin
            mat3_3n2 = 0.5 * mat_y.transpose() * mat1_n1n2;
            Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
            jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;

            // 1/2 B_z^T mat_z^T C B_lin
            mat3_3n2 = 0.5 * mat_z.transpose() * mat1_n1n2;
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
                
                // 1/2 B_lin^T C mat_x_i B_x_i
                mat1_n1n2 = 0.5 * material_mat * mat1_n1n2;
                Bmat_lin.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // 1/4 B_x^T mat_x^T C mat_x B_x
                mat3_3n2 = 0.5 * mat_x.transpose() * mat1_n1n2;
                Bmat_nl_x.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // 1/4 B_y^T mat_y^T C mat_x B_x
                mat3_3n2 = 0.5 * mat_y.transpose() * mat1_n1n2;
                Bmat_nl_y.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
                
                // 1/4 B_z^T mat_z^T C mat_x B_x
                mat3_3n2 = 0.5 * mat_z.transpose() * mat1_n1n2;
                Bmat_nl_z.right_multiply_transpose(mat2_n2n2, mat3_3n2);
                jac.topLeftCorner(n2, n2) += JxW[qp] * mat2_n2n2;
            }
            
            // use the stress to calculate the final contribution
            // to the Jacobian stiffness matrix
            mat4_33(0,0) = stress(0);
            mat4_33(1,1) = stress(1);
            mat4_33(2,2) = stress(2);
            mat4_33(0,1) = mat4_33(1,0) = 0.5*stress(3);
            mat4_33(1,2) = mat4_33(2,1) = 0.5*stress(4);
            mat4_33(0,2) = mat4_33(2,0) = 0.5*stress(5);

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
    
//    // place small values at diagonal of the rotational terms
//    if (request_jacobian) {
//        for (unsigned int i=n2/2; i<n2; i++)
//            jac(i,i) += 1.0e-12;
//    }
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::internal_residual_sensitivity(bool request_jacobian,
                                                         RealVectorX& f,
                                                         RealMatrixX& jac,
                                                         bool if_ignore_ho_jac) {
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::prestress_residual (bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac) {
    
    return request_jacobian;
}



bool
MAST::StructuralElement3D::prestress_residual_sensitivity (bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac) {
    
    return request_jacobian;
}




bool
MAST::StructuralElement3D::thermal_residual(bool request_jacobian,
                                            RealVectorX& f,
                                            RealMatrixX& jac,
                                            MAST::BoundaryConditionBase& bc) {
    
    const std::vector<Real>& JxW            = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = _fe->get_xyz();
    
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = 6,
    n2    = 3*n_phi;
    
    RealMatrixX
    material_exp_A_mat;
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec3_n2   = RealVectorX::Zero(n2),
    delta_t   = RealVectorX::Zero(1);
    
    
    libMesh::Point p;
    FEMOperatorMatrix Bmat;
    Bmat.reinit(n1, 3, n_phi); // three stress-strain components
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    mat = _property.thermal_expansion_A_matrix(*this);
    
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real t, t0;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        (*mat)       (p, _time, material_exp_A_mat);
        temp_func    (p, _time, t);
        ref_temp_func(p, _time, t0);
        delta_t(0) = t-t0;
        
        vec1_n1 = material_exp_A_mat * delta_t; // [C]{alpha (T - T0)}
        
        this->initialize_strain_operator(qp, Bmat);
        
        Bmat.vector_mult_transpose(vec3_n2, vec1_n1);
        
        f.topRows(n2) -= JxW[qp] * vec3_n2;
    }
    
    // Jacobian contribution from von Karman strain
    return false;
}



bool
MAST::StructuralElement3D::thermal_residual_sensitivity(bool request_jacobian,
                                                        RealVectorX& f,
                                                        RealMatrixX& jac,
                                                        MAST::BoundaryConditionBase& p) {
    
    return false;
}


void
MAST::StructuralElement3D::initialize_strain_operator (const unsigned int qp,
                                                       FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >&
    dphi = _fe->get_dphi();
    
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
    

    
    const std::vector<std::vector<libMesh::RealVectorValue> >&
    dphi = _fe->get_dphi();
    
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
    
    // nonlinear strain operator in x
    Bmat_nl_x.set_shape_function(0, 0, phi); // du/dx
    Bmat_nl_x.set_shape_function(1, 1, phi); // dv/dx
    Bmat_nl_x.set_shape_function(2, 2, phi); // dw/dx

    // nonlinear strain operator in u
    Bmat_nl_u.set_shape_function(0, 0, phi); // du/dx
    Bmat_nl_u.set_shape_function(1, 0, phi); // du/dy
    Bmat_nl_u.set_shape_function(2, 0, phi); // du/dz

    
    // dN/dy
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](1);
    // linear strain operator
    Bmat_lin.set_shape_function(1, 1, phi); //  epsilon_yy = dv/dy
    Bmat_lin.set_shape_function(3, 0, phi); //  gamma_xy = du/dy + ...
    Bmat_lin.set_shape_function(4, 2, phi); //  gamma_yz = dw/dy + ...
    
    // nonlinear strain operator in y
    Bmat_nl_y.set_shape_function(0, 0, phi); // du/dy
    Bmat_nl_y.set_shape_function(1, 1, phi); // dv/dy
    Bmat_nl_y.set_shape_function(2, 2, phi); // dw/dy

    // nonlinear strain operator in v
    Bmat_nl_v.set_shape_function(0, 1, phi); // dv/dx
    Bmat_nl_v.set_shape_function(1, 1, phi); // dv/dy
    Bmat_nl_v.set_shape_function(2, 1, phi); // dv/dz

    // dN/dz
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](2);
    Bmat_lin.set_shape_function(2, 2, phi); //  epsilon_xx = dw/dz
    Bmat_lin.set_shape_function(4, 1, phi); //  gamma_xy = dv/dz + ...
    Bmat_lin.set_shape_function(5, 0, phi); //  gamma_zx = du/dz + ...

    // nonlinear strain operator in z
    Bmat_nl_z.set_shape_function(0, 0, phi); // du/dz
    Bmat_nl_z.set_shape_function(1, 1, phi); // dv/dz
    Bmat_nl_z.set_shape_function(2, 2, phi); // dw/dz

    // nonlinear strain operator in w
    Bmat_nl_w.set_shape_function(0, 2, phi); // dw/dx
    Bmat_nl_w.set_shape_function(1, 2, phi); // dw/dy
    Bmat_nl_w.set_shape_function(2, 2, phi); // dw/dz

    
    // calculate the displacement gradient to create the
    RealVectorX
    ddisp_dx = RealVectorX::Zero(3),
    ddisp_dy = RealVectorX::Zero(3),
    ddisp_dz = RealVectorX::Zero(3);
    
    Bmat_nl_x.vector_mult(ddisp_dx, local_disp);
    Bmat_nl_y.vector_mult(ddisp_dy, local_disp);
    Bmat_nl_z.vector_mult(ddisp_dz, local_disp);

    // prepare the deformation gradient matrix
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
    mat_x.row(3) = 0.5*ddisp_dy;
    mat_x.row(5) = 0.5*ddisp_dz;
    
    mat_y.row(1) =     ddisp_dx;
    mat_y.row(3) = 0.5*ddisp_dx;
    mat_y.row(4) = 0.5*ddisp_dz;

    
    mat_z.row(2) =     ddisp_dx;
    mat_z.row(4) = 0.5*ddisp_dy;
    mat_z.row(5) = 0.5*ddisp_dx;

}




