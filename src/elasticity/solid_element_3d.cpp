
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
    n_phi              = (unsigned int)JxW.size(),
    n1                 =6,
    n2                 =3*n_phi;
    
    RealMatrixX
    material_mat,
    mat1_n1n2(n1, n2),
    mat2_n2n2(n2, n2);
    RealVectorX
    vec1_n1(n1),
    vec2_n2(n2);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> > mat_stiff =
    _property.stiffness_A_matrix(*this);
    
    libMesh::Point p;
    MAST::FEMOperatorMatrix Bmat;
    // six stress components, related to three displacements
    Bmat.reinit(n1, 3, _elem.n_nodes());
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff)(p, _time, material_mat);
        
        this->initialize_strain_operator(qp, Bmat);
        
        // calculate the stress
        Bmat.left_multiply(mat1_n1n2, material_mat);
        vec1_n1  = mat1_n1n2*_local_sol; // this is stress
        
        // now calculate the internal residual vector
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f = JxW[qp] * vec2_n2;
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            jac = JxW[qp] * mat2_n2n2;
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
    material_exp_A_mat,
    mat1_n1n2(n1, n2),
    mat2_n2n2(n2, n2),
    stress;
    RealVectorX
    vec1_n1(n1),
    vec2_n1(n1),
    vec3_n2(n2),
    delta_t(1);
    
    
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
        
        f += JxW[qp] * vec3_n2;
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
    RealVectorX phi(n_phi);
    
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


