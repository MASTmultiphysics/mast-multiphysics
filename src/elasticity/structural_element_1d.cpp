/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
#include "elasticity/structural_element_1d.h"
#include "numerics/fem_operator_matrix.h"
#include "mesh/local_elem_base.h"
#include "base/system_initialization.h"
#include "property_cards/element_property_card_base.h"
#include "base/boundary_condition_base.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/bending_operator.h"


MAST::StructuralElement1D::StructuralElement1D(MAST::SystemInitialization& sys,
                                               const libMesh::Elem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, elem, p) {
    
    // now initialize the finite element data structures
    _init_fe_operators(get_elem_for_quadrature(), &_fe, &_qrule, &_bending_operator);
}






void
MAST::StructuralElement1D::
initialize_direct_strain_operator(const unsigned int qp,
                                  const libMesh::FEBase& fe,
                                  MAST::FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi   = RealVectorX::Zero(n_phi);
    
    libmesh_assert_equal_to(Bmat.m(), 2);
    libmesh_assert_equal_to(Bmat.n(), 6*n_phi);
    
    // now set the shape function values
    // dN/dx
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = dphi[i_nd][qp](0);
    Bmat.set_shape_function(0, 0, phi); //  epsilon_xx = du/dx
    Bmat.set_shape_function(1, 3, phi); //  torsion operator = dtheta_x/dx
}



void
MAST::StructuralElement1D::
initialize_von_karman_strain_operator(const unsigned int qp,
                                      const libMesh::FEBase& fe,
                                      RealVectorX& vk_strain,
                                      RealMatrixX& vk_dvdxi_mat,
                                      RealMatrixX& vk_dwdxi_mat,
                                      MAST::FEMOperatorMatrix& Bmat_v_vk,
                                      MAST::FEMOperatorMatrix& Bmat_w_vk) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_strain.size(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat.rows(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat.cols(), 2);
    libmesh_assert_equal_to(Bmat_v_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_v_vk.n(), 6*n_phi);
    libmesh_assert_equal_to(vk_dwdxi_mat.rows(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat.cols(), 2);
    libmesh_assert_equal_to(Bmat_w_vk.m(), 2);
    libmesh_assert_equal_to(Bmat_w_vk.n(), 6*n_phi);
    
    Real dv=0., dw=0.;
    vk_strain.setConstant(0.);
    vk_dvdxi_mat.setConstant(0.);
    vk_dwdxi_mat.setConstant(0.);
    
    RealVectorX phi_vec   = RealVectorX::Zero(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);                // dphi/dx
        dv += phi_vec(i_nd)*_local_sol(n_phi+i_nd);   // dv/dx
        dw += phi_vec(i_nd)*_local_sol(2*n_phi+i_nd); // dw/dx
    }
    
    Bmat_v_vk.set_shape_function(0, 1, phi_vec); // dv/dx
    Bmat_w_vk.set_shape_function(0, 2, phi_vec); // dw/dx
    vk_dvdxi_mat(0, 0) = dv;                   // epsilon-xx : dv/dx
    vk_dwdxi_mat(0, 0) = dw;                   // epsilon-xx : dw/dx
    vk_strain(0) = 0.5*(dv*dv+dw*dw);          // 1/2 * [(dv/dx)^2 + (dw/dx)^2]
}




void
MAST::StructuralElement1D::
initialize_von_karman_strain_operator_sensitivity(const unsigned int qp,
                                                  const libMesh::FEBase& fe,
                                                  RealMatrixX& vk_dvdxi_mat_sens,
                                                  RealMatrixX& vk_dwdxi_mat_sens) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.rows(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.cols(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.rows(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.cols(), 2);
    
    Real dv=0., dw=0.;
    vk_dvdxi_mat_sens.setConstant(0.);
    vk_dwdxi_mat_sens.setConstant(0.);
    
    RealVectorX phi_vec   = RealVectorX::Zero(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);                // dphi/dx
        dv += phi_vec(i_nd)*_local_sol_sens(n_phi+i_nd);   // dv/dx
        dw += phi_vec(i_nd)*_local_sol_sens(2*n_phi+i_nd); // dw/dx
    }
    
    vk_dvdxi_mat_sens(0, 0) = dv;                   // epsilon-xx : dv/dx
    vk_dwdxi_mat_sens(0, 0) = dw;                   // epsilon-xx : dw/dx
}




bool
MAST::StructuralElement1D::calculate_stress(bool request_derivative,
                                            bool request_sensitivity,
                                            MAST::OutputFunctionBase& output) {
    
    // prepare the side finite element

    // ask the output object about the quadrature points at which
    // the stress evaluations need to be peformed
    MAST::PointwiseOutputEvaluationMode mode = output.evaluation_mode();
    
    std::vector<libMesh::Point> qp_loc;

    libMesh::FEBase         *fe_ptr     = NULL;
    libMesh::QBase          *qrule_ptr  = NULL;
    MAST::BendingOperator   *bend_ptr   = NULL;
    
    switch (mode) {
        case MAST::CENTROID: {
            qp_loc.resize(1);
            qp_loc[0] = libMesh::Point();
            _init_fe_operators(get_elem_for_quadrature(),
                               &fe_ptr,
                               &qrule_ptr,
                               &bend_ptr,
                               &qp_loc);
        }
            break;
            
        case MAST::SPECIFIED_QP: {
            qp_loc = output.get_qp_for_evaluation();
            _init_fe_operators(get_elem_for_quadrature(),
                               &fe_ptr,
                               &qrule_ptr,
                               &bend_ptr,
                               &qp_loc);
        }
            break;
            
        case MAST::ELEM_QP:
            // this will initialize the FE object at the points specified
            // by the quadrature rule
            _init_fe_operators(get_elem_for_quadrature(),
                               &fe_ptr,
                               &qrule_ptr,
                               &bend_ptr);
            qp_loc = qrule_ptr->get_points();
            break;
            
        defalut:
            // should not get here
            libmesh_error();
    }

    std::auto_ptr<libMesh::FEBase>         fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>          qrule(qrule_ptr);
    std::auto_ptr<MAST::BendingOperator1D>
    bend(dynamic_cast<MAST::BendingOperator1D*>(bend_ptr));

    
    // now that the FE object has been initialized, evaluate the stress values

    
    const std::vector<Real> &JxW              = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz    = fe->get_xyz();
    const unsigned int
    n_phi    = (unsigned int)fe->n_shape_functions(),
    n1       = this->n_direct_strain_components(),
    n2       = 6*n_phi,
    n3       = this->n_von_karman_strain_components();

    Real
    y = 0.,
    z = 0.;
    
    RealMatrixX
    material_A_mat,
    vk_dvdxi_mat = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat = RealMatrixX::Zero(n1,n3),
    dstrain_dX   = RealMatrixX::Zero(n1,n2),
    dstress_dX   = RealMatrixX::Zero(n1,n2),
    mat_n1n2     = RealMatrixX::Zero(n1,n2),
    eye          = RealMatrixX::Identity(n1, n1),
    dstrain_dX_3D= RealMatrixX::Zero(6,n2),
    dstress_dX_3D= RealMatrixX::Zero(6,n2);

    RealVectorX
    strain      = RealVectorX::Zero(n1),
    stress      = RealVectorX::Zero(n1),
    strain_bend = RealVectorX::Zero(n1),
    strain_vk   = RealVectorX::Zero(n3),
    strain_3D   = RealVectorX::Zero(6),
    stress_3D   = RealVectorX::Zero(6),
    dstrain_dp  = RealVectorX::Zero(n1),
    dstress_dp  = RealVectorX::Zero(n1);

    
    libMesh::Point p;
    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1,  _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy

    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _property.stiffness_A_matrix(*this);

    // get the thickness values for the bending strain calculation
    const MAST::FieldFunction<Real>
    &hy     =  _property.get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  _property.get<MAST::FieldFunction<Real> >("hz");

    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, fe->get_fe_type()) != MAST::NO_BENDING);
    
    // a reference to the stress output data structure
    MAST::StressStrainOutputBase& stress_output =
    dynamic_cast<MAST::StressStrainOutputBase&>(output);

    // TODO: improve the stress calculation by including shear stress due
    // to torsion and shear flow due to beam bending.
    
    ///////////////////////////////////////////////////////////////////////
    // second for loop to calculate the residual and stiffness contributions
    for (unsigned int qp=0; qp<qp_loc.size(); qp++) {
        
        _local_elem->global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff_A)(p, _time, material_A_mat);
        
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        Bmat_mem.vector_mult(strain, _local_sol);
        
        if (if_bending) {
            // von Karman strain
            if (if_vk) {  // get the vonKarman strain operator if needed
                
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
                                                            strain_vk,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                strain += strain_vk;
            }
            
            // add to this the bending strain
            // TODO: add coupling due to h_offset
            hy(p, _time, y);
            hz(p, _time, z);
            // TODO: this assumes isotropic section. Multilayered sections need
            // special considerations
            bend->initialize_bending_strain_operator_for_yz(*fe,
                                                            qp,
                                                            qp_loc[qp](1) * y/2.,
                                                            qp_loc[qp](2) * z/2.,
                                                            Bmat_bend);
            Bmat_bend.vector_mult(strain_bend, _local_sol);
            
            
            // add stress due to bending.
            strain += strain_bend;
        }

        
        // note that this assumes linear material laws
        stress = material_A_mat * strain;

        // now set the data for the 3D stress-strain vector
        // this is using only the direct strain/stress.
        // this can be improved by estimating the shear stresses from
        // torsion and shear flow from bending.
        strain_3D(0)  =   strain(0);
        stress_3D(0)  =   stress(0);
        
        // set the stress and strain data
        MAST::StressStrainOutputBase::Data&
        data = stress_output.add_stress_strain_at_qp_location(&_elem,
                                                              qp_loc[qp],
                                                              xyz[qp],
                                                              stress_3D,
                                                              strain_3D,
                                                              JxW[qp]);

        
        // calculate the derivative if requested
        if (request_derivative || request_sensitivity) {
            
            Bmat_mem.left_multiply(dstrain_dX, eye);  // membrane strain is linear
            
            if (if_bending) {
                
                // von Karman strain
                if (if_vk) {
                    
                    Bmat_v_vk.left_multiply(mat_n1n2, vk_dvdxi_mat);
                    dstrain_dX   +=  mat_n1n2;

                    Bmat_w_vk.left_multiply(mat_n1n2, vk_dwdxi_mat);
                    dstrain_dX   +=  mat_n1n2;
                }
                
                // bending strain
                Bmat_bend.left_multiply(mat_n1n2, eye);
                dstrain_dX  +=   mat_n1n2;
            }
            
            // note: this assumes lienar material laws
            dstress_dX  = material_A_mat * dstrain_dX;
            
            // copy to the 3D structure
            dstress_dX_3D.row(0)  = dstress_dX.row(0);
            dstrain_dX_3D.row(0)  = dstrain_dX.row(0);
            
            if (request_derivative)
                data.set_derivatives(dstress_dX_3D, dstrain_dX_3D);
                

            if (request_sensitivity) {
                // sensitivity of the response, s, is
                //   ds/dp   = partial s/partial p  +
                //             partial s/partial X   dX/dp
                //   the first part of the sensitivity is obtained from
                //
                // the first term includes direct sensitivity of the stress
                // with respect to the parameter, while holding the solution
                // constant. This should include influence of shape changes,
                // if the parameter is shape-dependent.
                // TODO: include shape sensitivity.
                // presently, only material parameter is included

                // get the material matrix sensitivity
                mat_stiff_A->derivative(MAST::PARTIAL_DERIVATIVE,
                                        *sensitivity_param,
                                        p,
                                        _time,
                                        material_A_mat);

                // partial sensitivity of strain is zero unless it is a
                // shape parameter.
                // TODO: shape sensitivity of strain operator
                
                // now use this to calculate the stress sensitivity.
                dstress_dp  =  material_A_mat * strain;
                dstrain_dp  =  RealVectorX::Zero(n1);
                
                //
                // use the derivative data to evaluate the second term in the
                // sensitivity
                //
                dstress_dp  += dstress_dX * _local_sol_sens;
                dstrain_dp  += dstrain_dX * _local_sol_sens;
                
                // copy the 3D object
                stress_3D(0) = dstress_dp(0);
                strain_3D(0) = dstrain_dp(0);
                
                // tell the data object about the sensitivity values
                data.set_sensitivity(sensitivity_param,
                                     stress_3D,
                                     strain_3D);
            }
        }
    }

    // if either derivative or sensitivity was requested, it was provided
    // by this routine
    return request_derivative || request_sensitivity;
}





bool
MAST::StructuralElement1D::internal_residual (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac,
                                              bool if_ignore_ho_jac)
{
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi    = (unsigned int)_fe->get_phi().size(),
    n1       = this->n_direct_strain_components(),
    n2       = 6*n_phi,
    n3       = this->n_von_karman_strain_components();
    
    RealMatrixX
    material_A_mat,
    material_B_mat,
    material_D_mat,
    mat1_n1n2    = RealMatrixX::Zero(n1,n2),
    mat2_n2n2    = RealMatrixX::Zero(n2,n2),
    mat3,
    mat4_n3n2    = RealMatrixX::Zero(n3,2),
    vk_dvdxi_mat = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat = RealMatrixX::Zero(n1,n3),
    stress       = RealMatrixX::Zero(2,2),
    stress_l     = RealMatrixX::Zero(2,2),
    local_jac    = RealMatrixX::Zero(n2,n2);
    
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2);
    
    local_f.setZero();
    local_jac.setZero();
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _property.stiffness_A_matrix(*this),
    mat_stiff_B  = _property.stiffness_B_matrix(*this),
    mat_stiff_D  = _property.stiffness_D_matrix(*this);
    
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff_A)(p, _time, material_A_mat);
        
        if (if_bending) {
            (*mat_stiff_B)(p, _time, material_B_mat);
            (*mat_stiff_D)(p, _time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_residual_operation(if_bending, if_vk, n2, qp, *_fe, JxW,
                                     request_jacobian, if_ignore_ho_jac,
                                     local_f, local_jac,
                                     Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
                                     stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                     material_A_mat,
                                     material_B_mat, material_D_mat, vec1_n1,
                                     vec2_n1, vec3_n2, vec4_n3,
                                     vec5_n3, mat1_n1n2, mat2_n2n2,
                                     mat3, mat4_n3n2);
        
    }
    
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_residual(request_jacobian,
                                                               local_f,
                                                               local_jac,
                                                               NULL);
    
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f += vec3_n2;

    if (request_jacobian) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac += mat2_n2n2;
    }
    
    return request_jacobian;
}





bool
MAST::StructuralElement1D::internal_residual_sensitivity (bool request_jacobian,
                                                          RealVectorX& f,
                                                          RealMatrixX& jac,
                                                          bool if_ignore_ho_jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi,
    n3    = this->n_von_karman_strain_components();
    
    RealMatrixX
    material_A_mat,
    material_B_mat,
    material_D_mat,
    material_trans_shear_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1,n2),
    mat2_n2n2     = RealMatrixX::Zero(n2,n2),
    mat3,
    mat4_n3n2     = RealMatrixX::Zero(n3,n2),
    vk_dvdxi_mat  = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat  = RealMatrixX::Zero(n1,n3),
    stress        = RealMatrixX::Zero(2,2),
    stress_l      = RealMatrixX::Zero(2,2),
    local_jac     = RealMatrixX::Zero(n2,n2);
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2);
    
    local_f.setZero();
    local_jac.setZero();

    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A = _property.stiffness_A_matrix(*this),
    mat_stiff_B = _property.stiffness_B_matrix(*this),
    mat_stiff_D = _property.stiffness_D_matrix(*this);
    
    libMesh::Point p;
    
    // first calculate the sensitivity due to the parameter
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        mat_stiff_A->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                p,
                                _time,
                                material_A_mat);
        
        if (if_bending) {
            mat_stiff_B->derivative(MAST::PARTIAL_DERIVATIVE,
                                    *this->sensitivity_param,
                                    p, _time, material_B_mat);
            mat_stiff_D->derivative(MAST::PARTIAL_DERIVATIVE,
                                    *this->sensitivity_param,
                                    p, _time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_residual_operation(if_bending, if_vk, n2, qp, *_fe, JxW,
                                     request_jacobian, if_ignore_ho_jac,
                                     local_f, local_jac,
                                     Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
                                     stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                     material_A_mat,
                                     material_B_mat, material_D_mat, vec1_n1,
                                     vec2_n1, vec3_n2, vec4_n3,
                                     vec5_n3, mat1_n1n2, mat2_n2n2,
                                     mat3, mat4_n3n2);
        
        // this accounts for the sensitivity of the linear stress as a result of
        // static solution. This is needed only for cases that require linearized
        // geometric stiffness matrix, for example in buckling analysis
        if (if_bending && if_vk && if_ignore_ho_jac && request_jacobian) {
            (*mat_stiff_A)(p, _time, material_A_mat);
            (*mat_stiff_B)(p, _time, material_B_mat);
            
            _linearized_geometric_stiffness_sensitivity_with_static_solution
            (n2, qp, *_fe, JxW, local_jac,
             Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk,
             stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
             material_A_mat, material_B_mat, vec1_n1,
             vec2_n1, mat1_n1n2, mat2_n2n2,
             mat3);
        }
        
    }
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (if_bending && _bending_operator->include_transverse_shear_energy())
        _bending_operator->calculate_transverse_shear_residual(request_jacobian,
                                                               local_f, local_jac,
                                                               this->sensitivity_param);
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f += vec3_n2;
    if (request_jacobian) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac += mat2_n2n2;
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElement1D::
internal_residual_jac_dot_state_sensitivity (RealMatrixX& jac) {
    
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW            = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi,
    n3    = this->n_von_karman_strain_components();
    
    RealMatrixX
    material_A_mat,
    material_B_mat,
    material_D_mat,
    mat1_n1n2     = RealMatrixX::Zero(n1,n2),
    mat2_n2n2     = RealMatrixX::Zero(n2,n2),
    mat3,
    vk_dvdxi_mat_sens = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat_sens = RealMatrixX::Zero(n1,n3),
    mat4_n3n2         = RealMatrixX::Zero(n3,n2),
    vk_dvdxi_mat      = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat      = RealMatrixX::Zero(n1,n3),
    stress            = RealMatrixX::Zero(2,2),
    stress_l          = RealMatrixX::Zero(2,2),
    local_jac         = RealMatrixX::Zero(n2,n2);
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2);

    local_f.setZero();
    local_jac.setZero();


    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    // without the nonlinear strain, this matrix is zero.
    if (!if_vk)
        return false;
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _property.stiffness_A_matrix(*this),
    mat_stiff_B  = _property.stiffness_B_matrix(*this),
    mat_stiff_D  = _property.stiffness_D_matrix(*this);
    
    
    libMesh::Point p;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], p);
        
        // get the material matrix
        (*mat_stiff_A)(p, _time, material_A_mat);
        
        (*mat_stiff_B)(p, _time, material_B_mat);
        (*mat_stiff_D)(p, _time, material_D_mat);
        
        // now calculte the quantity for these matrices
        this->initialize_direct_strain_operator(qp, *_fe, Bmat_mem);
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        Bmat_mem.vector_mult(vec1_n1, _local_sol_sens);
        vec2_n1 = material_A_mat * vec1_n1; // linear direct stress
        
        // copy the stress values to a matrix
        stress(0,0)   = vec2_n1(0); // sigma_xx
        
        // get the bending strain operator
        vec2_n1.setConstant(0.); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(*_fe, qp, Bmat_bend);
            
            //  evaluate the bending stress and add that to the stress vector
            // for evaluation in the nonlinear stress term
            Bmat_bend.vector_mult(vec2_n1, _local_sol_sens);
            vec1_n1 = material_B_mat * vec2_n1;
            stress(0,0)   += vec1_n1(0);
            
            if (if_vk) {  // get the vonKarman strain operator if needed
                
                this->initialize_von_karman_strain_operator(qp,
                                                            *_fe,
                                                            vec2_n1,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                this->initialize_von_karman_strain_operator_sensitivity(qp,
                                                                        *_fe,
                                                                        vk_dvdxi_mat_sens,
                                                                        vk_dwdxi_mat_sens);
                // sensitivity of von Karman strain
                vec2_n1.setConstant(0.);
                vec2_n1(0) = (vk_dvdxi_mat(0,0)*vk_dvdxi_mat_sens(0,0) +
                              vk_dwdxi_mat(0,0)*vk_dwdxi_mat_sens(0,0));
                vec1_n1 = material_A_mat * vec2_n1;
                stress(0,0) += vec1_n1(0);
            }
        }
        
        // copy the stress to use here.
        vec1_n1.setConstant(0.);
        
        // now calculate the matrix
        // membrane - vk: v-displacement
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat_sens);
        mat3 = material_A_mat * mat3;
        Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // membrane - vk: w-displacement
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat_sens);
        mat3 = material_A_mat * mat3;
        Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - membrane: v-displacement
        Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - membrane: w-displacement
        Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
        mat3 = vk_dwdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - vk: v-displacement
        mat3 = RealMatrixX::Zero(2, n2);
        Bmat_v_vk.left_multiply(mat3, stress);
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * material_A_mat * mat3;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat_sens);
        mat3 = vk_dvdxi_mat.transpose() * material_A_mat * mat3;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - vk: w-displacement
        mat3 = RealMatrixX::Zero(2, n2);
        Bmat_w_vk.left_multiply(mat3, stress);
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
        mat3 = vk_dwdxi_mat_sens.transpose() * material_A_mat * mat3;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat_sens);
        mat3 = vk_dwdxi_mat.transpose() * material_A_mat * mat3;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // coupling of v, w-displacements
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat_sens);
        mat3 = vk_dvdxi_mat.transpose() * material_A_mat * mat3;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * material_A_mat * mat3;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat_sens);
        mat3 = vk_dwdxi_mat.transpose() * material_A_mat * mat3;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
        mat3 = vk_dwdxi_mat_sens.transpose() * material_A_mat * mat3;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // bending - vk: v-displacement
        mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
        Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat_sens);
        mat3 = material_B_mat.transpose() * mat3;
        Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // bending - vk: w-displacement
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat_sens);
        mat3 = material_B_mat.transpose() * mat3;
        Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - bending: v-displacement
        Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        // vk - bending: w-displacement
        Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
        mat3 = vk_dwdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
    }
    
    transform_matrix_to_global_system(local_jac, mat2_n2n2);
    jac += mat2_n2n2;
    
    return true;
}




void
MAST::StructuralElement1D::
_internal_residual_operation(bool if_bending,
                             bool if_vk,
                             const unsigned int n2,
                             const unsigned int qp,
                             const libMesh::FEBase& fe,
                             const std::vector<Real>& JxW,
                             bool request_jacobian,
                             bool if_ignore_ho_jac,
                             RealVectorX& local_f,
                             RealMatrixX& local_jac,
                             MAST::FEMOperatorMatrix& Bmat_mem,
                             MAST::FEMOperatorMatrix& Bmat_bend,
                             MAST::FEMOperatorMatrix& Bmat_v_vk,
                             MAST::FEMOperatorMatrix& Bmat_w_vk,
                             RealMatrixX& stress,
                             RealMatrixX& stress_l,
                             RealMatrixX& vk_dvdxi_mat,
                             RealMatrixX& vk_dwdxi_mat,
                             RealMatrixX& material_A_mat,
                             RealMatrixX& material_B_mat,
                             RealMatrixX& material_D_mat,
                             RealVectorX& vec1_n1,
                             RealVectorX& vec2_n1,
                             RealVectorX& vec3_n2,
                             RealVectorX& vec4_2,
                             RealVectorX& vec5_2,
                             RealMatrixX& mat1_n1n2,
                             RealMatrixX& mat2_n2n2,
                             RealMatrixX& mat3,
                             RealMatrixX& mat4_2n2)
{
    this->initialize_direct_strain_operator(qp, fe, Bmat_mem);
    
    // first handle constant throught the thickness stresses: membrane and vonKarman
    Bmat_mem.vector_mult(vec1_n1, _local_sol);
    vec2_n1 = material_A_mat * vec1_n1; // linear direct stress
    
    // copy the stress values to a matrix
    stress_l(0,0) = vec2_n1(0); // sigma_xx
    stress(0,0)   = vec2_n1(0);
    stress(1,1)   = vec1_n1(0); // temporary storage of the axial strain
    stress(0,1)   = vec2_n1(1); // temporary storage of the torsion force
    
    // get the bending strain operator
    vec2_n1.setConstant(0.); // used to store vk strain, if applicable
    if (if_bending) {
        _bending_operator->initialize_bending_strain_operator(fe, qp, Bmat_bend);
        
        //  evaluate the bending stress and add that to the stress vector
        // for evaluation in the nonlinear stress term
        Bmat_bend.vector_mult(vec2_n1, _local_sol);
        vec1_n1 = material_B_mat * vec2_n1;
        stress_l(0,0) += vec1_n1(0);
        stress(0,0)   += vec1_n1(0);
        
        if (if_vk) {  // get the vonKarman strain operator if needed
            
            this->initialize_von_karman_strain_operator(qp,
                                                        fe,
                                                        vec2_n1, // epsilon_vk
                                                        vk_dvdxi_mat,
                                                        vk_dwdxi_mat,
                                                        Bmat_v_vk,
                                                        Bmat_w_vk);
            vec1_n1 = material_A_mat * vec2_n1;
            stress(0,0) += vec1_n1(0); // total strain that multiplies with the membrane strain
            stress(1,1) += vec2_n1(0); // add the two strains to get the direct strain
        }
    }
    
    // copy the stress to use here.
    vec1_n1.setConstant(0.);
    vec1_n1(0)   = stress(0,0);
    vec1_n1(1)   = stress(0,1); // use the torsion strain from the temporary location
    stress(0, 1) = 0.;   // zero out the temporary value storing the torsion strain
    
    // now the internal force vector
    // this includes the membrane strain operator with all A and B material operators
    Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
    local_f += JxW[qp] * vec3_n2;
    
    if (if_bending) {
        if (if_vk) {
            // von Karman strain: direct stress
            vec4_2 = vk_dvdxi_mat.transpose() * vec1_n1;
            Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
            local_f += JxW[qp] * vec3_n2;
            
            // von Karman strain: direct stress
            vec4_2 = vk_dwdxi_mat.transpose() * vec1_n1;
            Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
            local_f += JxW[qp] * vec3_n2;
        }
        
        // use the direct strain from the temprary storage
        vec2_n1(0)  = stress(1,1);
        stress(1,1) = 0.;
        // now coupling with the bending strain
        // B_bend^T [B] B_mem
        vec1_n1 = material_B_mat * vec2_n1;
        Bmat_bend.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        
        // now bending stress
        Bmat_bend.vector_mult(vec2_n1, _local_sol);
        vec1_n1 = material_D_mat * vec2_n1;
        Bmat_bend.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
    }
    
    if (request_jacobian) {
        // membrane - membrane
        Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
        Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
        local_jac += JxW[qp] * mat2_n2n2;
                
        if (if_bending) {
            if (if_vk) {
                // membrane - vk: v-displacement
                mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
                Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                mat3 = material_A_mat * mat3;
                Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // membrane - vk: w-displacement
                mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
                Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3 = material_A_mat * mat3;
                Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - membrane: v-displacement
                Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
                mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - membrane: w-displacement
                Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
                mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // if only the first order term of the Jacobian is needed, for
                // example for linearized buckling analysis, then the linear
                // stress combined with the variation of the von Karman strain
                // is included. Otherwise, all terms are included
                if (if_ignore_ho_jac) {
                    // vk - vk: v-displacement: first order term
                    mat3 = RealMatrixX::Zero(2, n2);
                    Bmat_v_vk.left_multiply(mat3, stress_l);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    // vk - vk: v-displacement: first order term
                    mat3 = RealMatrixX::Zero(2, n2);
                    Bmat_w_vk.left_multiply(mat3, stress_l);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                }
                else {
                    // vk - vk: v-displacement
                    mat3 = RealMatrixX::Zero(2, n2);
                    Bmat_v_vk.left_multiply(mat3, stress);
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
                    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                    mat3 = vk_dvdxi_mat.transpose() * material_A_mat * mat3;
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    // vk - vk: w-displacement
                    mat3 = RealMatrixX::Zero(2, n2);
                    Bmat_w_vk.left_multiply(mat3, stress);
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
                    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                    mat3 = vk_dwdxi_mat.transpose() * material_A_mat * mat3;
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    // coupling of v, w-displacements
                    mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
                    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                    mat3 = vk_dvdxi_mat.transpose() * material_A_mat * mat3;
                    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                    mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
                    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                    mat3 = vk_dwdxi_mat.transpose() * material_A_mat * mat3;
                    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                    local_jac += JxW[qp] * mat2_n2n2;
                    
                }
                
                // bending - vk: v-displacement
                mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
                Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
                mat3 = material_B_mat.transpose() * mat3;
                Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // bending - vk: w-displacement
                mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
                Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3 = material_B_mat.transpose() * mat3;
                Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - bending: v-displacement
                Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - bending: w-displacement
                Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
            
            // bending - membrane
            Bmat_mem.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_bend.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;
            
            // membrane - bending
            Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;
            
            // bending - bending
            Bmat_bend.left_multiply(mat1_n1n2, material_D_mat);
            Bmat_bend.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;
        }
    }
}





void
MAST::StructuralElement1D::
_linearized_geometric_stiffness_sensitivity_with_static_solution
(const unsigned int n2,
 const unsigned int qp,
 const libMesh::FEBase& fe,
 const std::vector<Real>& JxW,
 RealMatrixX& local_jac,
 MAST::FEMOperatorMatrix& Bmat_mem,
 MAST::FEMOperatorMatrix& Bmat_bend,
 MAST::FEMOperatorMatrix& Bmat_v_vk,
 MAST::FEMOperatorMatrix& Bmat_w_vk,
 RealMatrixX& stress_l,
 RealMatrixX& vk_dvdxi_mat,
 RealMatrixX& vk_dwdxi_mat,
 RealMatrixX& material_A_mat,
 RealMatrixX& material_B_mat,
 RealVectorX& vec1_n1,
 RealVectorX& vec2_n1,
 RealMatrixX& mat1_n1n2,
 RealMatrixX& mat2_n2n2,
 RealMatrixX& mat3) {
    
    this->initialize_direct_strain_operator(qp, fe, Bmat_mem);
    _bending_operator->initialize_bending_strain_operator(fe, qp, Bmat_bend);
    
    // first handle constant throught the thickness stresses: membrane and vonKarman
    Bmat_mem.vector_mult(vec1_n1, _local_sol);
    vec2_n1 = material_A_mat * vec1_n1; // linear direct stress
    
    // copy the stress values to a matrix
    stress_l(0,0) = vec2_n1(0); // sigma_xx
    
    // get the von Karman operator matrix
    this->initialize_von_karman_strain_operator(qp,
                                                fe,
                                                vec2_n1, // epsilon_vk
                                                vk_dvdxi_mat,
                                                vk_dwdxi_mat,
                                                Bmat_v_vk,
                                                Bmat_w_vk);
    
    // sensitivity of the vk_dwdxi matrix due to solution sensitivity
    this->initialize_von_karman_strain_operator_sensitivity(qp,
                                                            fe,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat);
    
    
    // membrane - vk: v-displacement
    mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
    mat3 = material_A_mat * mat3;
    Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // membrane - vk: w-displacement
    mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3 = material_A_mat * mat3;
    Bmat_mem.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - membrane: v-displacement
    Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
    mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - membrane: w-displacement
    Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
    mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - vk: v-displacement: first order term
    mat3 = RealMatrixX::Zero(2, n2);
    Bmat_v_vk.left_multiply(mat3, stress_l);
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - vk: v-displacement: first order term
    mat3 = RealMatrixX::Zero(2, n2);
    Bmat_w_vk.left_multiply(mat3, stress_l);
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // bending - vk: v-displacement
    mat3 = RealMatrixX::Zero(vk_dvdxi_mat.rows(), n2);
    Bmat_v_vk.left_multiply(mat3, vk_dvdxi_mat);
    mat3 = material_B_mat.transpose() * mat3;
    Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // bending - vk: w-displacement
    mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
    Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
    mat3 = material_B_mat.transpose() * mat3;
    Bmat_bend.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - bending: v-displacement
    Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
    mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
    Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
    
    // vk - bending: w-displacement
    Bmat_bend.left_multiply(mat1_n1n2, material_B_mat);
    mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
    Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
    local_jac += JxW[qp] * mat2_n2n2;
}


void
MAST::StructuralElement1D::_convert_prestress_A_mat_to_vector(const RealMatrixX& mat,
                                                              RealVectorX& vec) const {
    
    libmesh_assert_equal_to(mat.rows(), 2);
    libmesh_assert_equal_to(mat.cols(), 2);
    vec = RealVectorX::Zero(2);
    vec(0) = mat(0,0);
}


void
MAST::StructuralElement1D::_convert_prestress_B_mat_to_vector(const RealMatrixX& mat,
                                                              RealVectorX& vec) const {
    
    libmesh_assert_equal_to(mat.rows(), 2);
    libmesh_assert_equal_to(mat.cols(), 2);
    vec = RealVectorX::Zero(2);
    vec(0) = mat(0,0);
    vec(1) = mat(0,1);
}



bool
MAST::StructuralElement1D::prestress_residual (bool request_jacobian,
                                               RealVectorX& f,
                                               RealMatrixX& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    =6*n_phi,
    n3    = this->n_von_karman_strain_components();
    
    RealMatrixX
    mat2_n2n2     = RealMatrixX::Zero(n2, n2),
    mat3,
    vk_dvdxi_mat  = RealMatrixX::Zero(n1, n3),
    vk_dwdxi_mat  = RealMatrixX::Zero(n1, n3),
    local_jac     = RealMatrixX::Zero(n2, n2),
    prestress_mat_A,
    prestress_mat_B;
    RealVectorX
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2),
    prestress_vec_A,
    prestress_vec_B;
    
    local_f.setZero();
    local_jac.setZero();
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    prestress_A = _property.prestress_A_matrix(*this),
    prestress_B = _property.prestress_B_matrix(*this);
    
    libMesh::Point p;
    
    // now calculate the quantity
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], p);
        
        (*prestress_A)(p, _time, prestress_mat_A);
        (*prestress_B)(p, _time, prestress_mat_B);
        _convert_prestress_A_mat_to_vector(prestress_mat_A, prestress_vec_A);
        _convert_prestress_B_mat_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, *_fe, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.setConstant(0.); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(*_fe, qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *_fe,
                                                            vec2_n1,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f += JxW[qp] * vec3_n2; // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain: v-displacement
                vec4_n3 = vk_dvdxi_mat.transpose() * prestress_vec_A;
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f += JxW[qp] * vec3_n2; // epsilon_vk * sigma_0
                
                // von Karman strain: w-displacement
                vec4_n3 = vk_dwdxi_mat.transpose() * prestress_vec_A;
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f += JxW[qp] * vec3_n2; // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                // v-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_v_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // w-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_w_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f += vec3_n2;
    if (request_jacobian && if_vk) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac += mat2_n2n2;
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::StructuralElement1D::prestress_residual_sensitivity (bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi,
    n3    = this->n_von_karman_strain_components();

    RealMatrixX
    mat2_n2n2     = RealMatrixX::Zero(n2,n2),
    mat3,
    vk_dwdxi_mat  = RealMatrixX::Zero(n1,n3),
    vk_dvdxi_mat  = RealMatrixX::Zero(n1,n3),
    local_jac     = RealMatrixX::Zero(n2,n2),
    prestress_mat_A,
    prestress_mat_B;
    RealVectorX
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2),
    prestress_vec_A,
    prestress_vec_B;
    
    local_f.setZero();
    local_jac.setZero();

    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX> >
    prestress_A = _property.prestress_A_matrix(*this),
    prestress_B = _property.prestress_B_matrix(*this);
    
    libMesh::Point p;
    
    // transform to the local coordinate system
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], p);
        
        prestress_A->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                p, _time, prestress_mat_A);
        prestress_B->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                p, _time, prestress_mat_B);
        _convert_prestress_A_mat_to_vector(prestress_mat_A, prestress_vec_A);
        _convert_prestress_B_mat_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, *_fe, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.setConstant(0.); // used to store vk strain, if applicable
        if (if_bending) {
            _bending_operator->initialize_bending_strain_operator(*_fe, qp, Bmat_bend);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *_fe,
                                                            vec2_n1,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
        }
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        // multiply this with the constant through the thickness strain
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, prestress_vec_A);
        local_f += JxW[qp] * vec3_n2; // epsilon_mem * sigma_0
        
        if (if_bending) {
            if (if_vk) {
                // von Karman strain: v-displacement
                vec4_n3 = vk_dvdxi_mat.transpose() * prestress_vec_A;
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f += JxW[qp] * vec3_n2; // epsilon_vk * sigma_0
                
                // von Karman strain: w-displacement
                vec4_n3 = vk_dwdxi_mat.transpose() * prestress_vec_A;
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_n3);
                local_f += JxW[qp] * vec3_n2; // epsilon_vk * sigma_0
            }
            
            // now coupling with the bending strain
            Bmat_bend.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (if_bending && if_vk) {
                // v-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_v_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // w-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_w_vk.left_multiply(mat3, prestress_mat_A);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f += vec3_n2;
    if (request_jacobian && if_vk) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac += mat2_n2n2;
    }
    
    // only the nonlinear strain returns a Jacobian for prestressing
    return (request_jacobian && if_vk);
}





bool
MAST::StructuralElement1D::thermal_residual (bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac,
                                             MAST::BoundaryConditionBase& bc)
{
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi,
    n3    = this->n_von_karman_strain_components();
    
    RealMatrixX
    material_exp_A_mat,
    material_exp_B_mat,
    mat1_n1n2    = RealMatrixX::Zero(n1,n2),
    mat2_n2n2    = RealMatrixX::Zero(n2,n2),
    mat3,
    mat4_n3n2    = RealMatrixX::Zero(n3,n2),
    vk_dvdxi_mat = RealMatrixX::Zero(n1,n3),
    vk_dwdxi_mat = RealMatrixX::Zero(n1,n3),
    stress       = RealMatrixX::Zero(2,2),
    local_jac    = RealMatrixX::Zero(n2, n2);
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_2     = RealVectorX::Zero(2),
    vec5_n3    = RealVectorX::Zero(n3),
    local_f    = RealVectorX::Zero(n2),
    delta_t    = RealVectorX::Zero(1);
    
    local_f.setZero();
    local_jac.setZero();

    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    expansion_A = _property.thermal_expansion_A_matrix(*this),
    expansion_B = _property.thermal_expansion_B_matrix(*this);
    
    // temperature function
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real t, t0;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], pt);
        
        // get the material property
        (*expansion_A)(pt, _time, material_exp_A_mat);
        (*expansion_B)(pt, _time, material_exp_B_mat);
        
        // get the temperature function
        temp_func    (pt, _time, t);
        ref_temp_func(pt, _time, t0);
        delta_t(0) = t-t0;
        
        vec1_n1 = material_exp_A_mat * delta_t; // [C]{alpha (T - T0)} (with membrane strain)
        stress(0,0) = vec1_n1(0); // sigma_xx
        vec2_n1 = material_exp_B_mat * delta_t; // [C]{alpha (T - T0)} (with bending strain)
        
        this->initialize_direct_strain_operator(qp, *_fe, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        
        if (if_bending) {
            // bending strain
            _bending_operator->initialize_bending_strain_operator(*_fe, qp, Bmat_bend);
            Bmat_bend.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;
            
            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *_fe,
                                                            vec2_n1, // epsilon_vk
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                // von Karman strain: v-displacement
                vec4_2 = vk_dvdxi_mat.transpose() * vec1_n1;
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f += JxW[qp] * vec3_n2;
                
                // von Karman strain: w-displacement
                vec4_2 = vk_dwdxi_mat.transpose() * vec1_n1;
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f += JxW[qp] * vec3_n2;
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                
                // vk - vk: v-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_v_vk.left_multiply(mat3, stress);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - vk: w-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_w_vk.left_multiply(mat3, stress);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f -= vec3_n2;
    if (request_jacobian && if_vk) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac -= mat2_n2n2;
    }
    
    // Jacobian contribution from von Karman strain
    return request_jacobian && if_vk;
}




bool
MAST::StructuralElement1D::thermal_residual_sensitivity (bool request_jacobian,
                                                         RealVectorX& f,
                                                         RealMatrixX& jac,
                                                         MAST::BoundaryConditionBase& bc)
{
    MAST::FEMOperatorMatrix Bmat_mem, Bmat_bend, Bmat_v_vk, Bmat_w_vk;
    
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = _fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)_fe->get_phi().size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi,
    n3    = this->n_von_karman_strain_components();
    
    RealMatrixX
    material_exp_A_mat,
    material_exp_B_mat,
    material_exp_A_mat_sens,
    material_exp_B_mat_sens,
    mat1_n1n2     = RealMatrixX::Zero(n1,n2),
    mat2_n2n2     = RealMatrixX::Zero(n2,n2),
    mat3,
    mat4_n3n2     = RealMatrixX::Zero(n3,n2),
    vk_dvdxi_mat  = RealMatrixX::Zero(2,2),
    vk_dwdxi_mat  = RealMatrixX::Zero(2,2),
    stress        = RealMatrixX::Zero(2,2),
    local_jac     = RealMatrixX::Zero(n2,n2);
    RealVectorX
    vec1_n1      = RealVectorX::Zero(n1),
    vec2_n1      = RealVectorX::Zero(n1),
    vec3_n2      = RealVectorX::Zero(n2),
    vec4_2       = RealVectorX::Zero(2),
    vec5_n1      = RealVectorX::Zero(n1),
    local_f      = RealVectorX::Zero(n2),
    delta_t      = RealVectorX::Zero(1),
    delta_t_sens = RealVectorX::Zero(1);
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::VON_KARMAN_STRAIN),
    if_bending = (_property.bending_model(_elem, _fe->get_fe_type()) != MAST::NO_BENDING);
    
    std::auto_ptr<MAST::FieldFunction<RealMatrixX > >
    expansion_A = _property.thermal_expansion_A_matrix(*this),
    expansion_B = _property.thermal_expansion_B_matrix(*this);
    
    // temperature function
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real t, t0, t_sens;
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        this->local_elem().global_coordinates_location(xyz[qp], pt);
        
        // get the material property
        (*expansion_A)(pt, _time, material_exp_A_mat);
        (*expansion_B)(pt, _time, material_exp_B_mat);
        expansion_A->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                pt, _time, material_exp_A_mat_sens);
        expansion_B->derivative(MAST::PARTIAL_DERIVATIVE,
                                *this->sensitivity_param,
                                pt, _time, material_exp_B_mat_sens);
        
        // get the temperature function
        temp_func(pt, _time, t);
        temp_func.derivative(MAST::PARTIAL_DERIVATIVE,
                             *this->sensitivity_param,
                             pt, _time, t_sens);
        ref_temp_func(pt, _time, t0);
        delta_t(0)      = t-t0;
        delta_t_sens(0) = t_sens;
        
        // now prepare the membrane force sensitivity
        vec1_n1 = material_exp_A_mat * delta_t_sens; // [C]{alpha (dT/dp)} (with membrane strain)
        vec2_n1 = material_exp_A_mat_sens * delta_t; // d([C].{alpha})/dp (T - T0)} (with membrane
        vec1_n1 += vec2_n1;  // sensitivity of the thermal membrane force
        stress(0,0) = vec1_n1(0); // sigma_xx
        
        // now prepare the membrane-bending coupling force sensitivity
        vec2_n1 = material_exp_B_mat * delta_t_sens; // [C]{alpha dT/dp} (with bending strain)
        vec5_n1 = material_exp_B_mat_sens * delta_t; // d([C].{alpha})/dp (T - T0)} (with bending strain)
        vec2_n1 += vec5_n1;  // sensitivity of the thermal membrane force
        
        
        this->initialize_direct_strain_operator(qp, *_fe, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        
        if (if_bending) {
            // bending strain
            _bending_operator->initialize_bending_strain_operator(*_fe, qp, Bmat_bend);
            Bmat_bend.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;
            
            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *_fe,
                                                            vec2_n1, // epsilon_vk
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                // von Karman strain: v-displacement
                vec4_2 = vk_dvdxi_mat.transpose() * vec1_n1;
                Bmat_v_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f += JxW[qp] * vec3_n2;
                
                // von Karman strain: w-displacement
                vec4_2 = vk_dwdxi_mat.transpose() * vec1_n1;
                Bmat_w_vk.vector_mult_transpose(vec3_n2, vec4_2);
                local_f += JxW[qp] * vec3_n2;
            }
            
            if (request_jacobian && if_vk) { // Jacobian only for vk strain
                                             // vk - vk: v-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_v_vk.left_multiply(mat3, stress);
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - vk: w-displacement
                mat3 = RealMatrixX::Zero(2, n2);
                Bmat_w_vk.left_multiply(mat3, stress);
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
        }
    }
    
    
    // now transform to the global coorodinate system
    transform_vector_to_global_system(local_f, vec3_n2);
    f -= vec3_n2;
    if (request_jacobian && if_vk) {
        transform_matrix_to_global_system(local_jac, mat2_n2n2);
        jac -= mat2_n2n2;
    }
    
    // Jacobian contribution from von Karman strain
    return request_jacobian && if_vk;
}







bool
MAST::StructuralElement1D::
piston_theory_residual(bool request_jacobian,
                       RealVectorX &f,
                       RealMatrixX& jac_xdot,
                       RealMatrixX& jac,
                       MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    
    const std::vector<Real> &JxW                = _fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = _fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = _fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi;
    
    
    
    // convert to piston theory boundary condition so that the necessary
    // flow properties can be obtained
    const MAST::PistonTheoryBoundaryCondition& piston_bc =
    dynamic_cast<MAST::PistonTheoryBoundaryCondition&>(bc);
    
    const Real
    mach      =   piston_bc.mach(),
    a_inf     =   piston_bc.a_inf(),
    gamma     =   piston_bc.gamma(),
    rho       =   piston_bc.rho(),
    q_dyn     =   .5*rho*pow(mach*a_inf ,2);
    
    const unsigned int
    order     =   piston_bc.order();
    
    FEMOperatorMatrix
    Bmat_v,         // operator matrix for the v-displacement
    Bmat_w,         // operator matrix for the w-displacement
    Bmat_dvdx,      // operator matrix to calculate the derivativ of v wrt x
    Bmat_dwdx;      // operator matrix to calculate the derivativ of w wrt x

    Bmat_dvdx.reinit(2, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_dwdx.reinit(2, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    RealVectorX
    phi_vec  = RealVectorX::Zero(n_phi),
    force    = RealVectorX::Zero(n1),
    local_f  = RealVectorX::Zero(n2),
    vec_n1   = RealVectorX::Zero(n1),
    vec_n2   = RealVectorX::Zero(n2),
    vel      = RealVectorX::Zero(3),
    cp       = RealVectorX::Zero(3),
    normal   = RealVectorX::Zero(3),
    dummy    = RealVectorX::Zero(2);
    

    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = -1.;
    
    
    RealMatrixX
    dvdx             = RealMatrixX::Zero(2,2),
    dwdx             = RealMatrixX::Zero(2,2),
    local_jac        = RealMatrixX::Zero(n2,n2),
    local_jac_xdot   = RealMatrixX::Zero(n2,n2),
    mat_n2n2         = RealMatrixX::Zero(n2,n2);
    
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        
        //_local_elem->global_coordinates_location(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat_v.reinit(n1, _system.n_vars(), n_phi);
        Bmat_v.set_shape_function(0, 1, phi_vec);
        Bmat_w.reinit(n1, _system.n_vars(), n_phi);
        Bmat_w.set_shape_function(1, 2, phi_vec);

        // use the Bmat to calculate the velocity vector
        Bmat_v.right_multiply(vec_n1, _local_vel);
        vel(1)     = vec_n1(0);   // v_dot
        Bmat_w.right_multiply(vec_n1, _local_vel);
        vel(2)     = vec_n1(1);   // w_dot
        
        // get the operators for dv/dx and dw/dx. We will use the
        // von Karman strain operators for this
        initialize_von_karman_strain_operator(qp,
                                              *_fe,
                                              dummy,
                                              dvdx,
                                              dwdx,
                                              Bmat_dvdx,
                                              Bmat_dwdx);
        
        // add the normal velocity due to the deformation slope
        vel(1)  += (mach*a_inf)*dvdx(0,0);  // this now is v_dot+dv/dx*(mach*a_inf)
        vel(2)  += (mach*a_inf)*dwdx(0,0);  // this now is w_dot+dw/dx*(mach*a_inf)
        
        
        // now use this information to evaluate the normal cp due to
        // both the v and w motions
        for (unsigned int i=1; i<3; i++)
            cp(i)  = piston_theory_cp(order, vel(i), a_inf, gamma, mach);
        
        // calculate force and discrete vector for v-disp
        force.setZero();
        force(0) = cp(1) * normal(1);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;

        // now the w-disp
        force.setZero();
        force(1) = cp(2) * normal(2);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;
        
        
        // calculate the Jacobian if requested
        if (request_jacobian) {
            
            // we need the derivative of cp wrt normal velocity
            for (unsigned int i=1; i<3; i++)
                cp(i)  = piston_theory_dcp_dvn(order, vel(i), a_inf, gamma, mach);
         
            // calculate the component of Jacobian due to v-velocity
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_v);
            local_jac_xdot += (JxW[qp] * cp(1) * normal(1)) * mat_n2n2;
            
            // calculate the component of Jacobian due to w-velocity
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_w);
            local_jac_xdot += (JxW[qp] * cp(2) * normal(2)) * mat_n2n2;
            
            
            // now use calculate the component of Jacobian due to deformation
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_dvdx);  // v: B^T dB/dx
            local_jac += (JxW[qp] * cp(1) * a_inf * mach * normal(1)) * mat_n2n2;
            
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_dwdx);  // w: B^T dB/dx
            local_jac += (JxW[qp] * cp(2) * a_inf * mach * normal(2)) * mat_n2n2;
        }
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= q_dyn * vec_n2;
    
    // if the Jacobian was requested, then transform it and add to the
    // global Jacobian
    if (request_jacobian) {
        transform_matrix_to_global_system(local_jac_xdot, mat_n2n2);
        jac_xdot -= q_dyn * mat_n2n2;

        transform_matrix_to_global_system(local_jac, mat_n2n2);
        jac      -= q_dyn * mat_n2n2;
    }
    
    return (request_jacobian);
}










bool
MAST::StructuralElement1D::
piston_theory_residual_sensitivity(bool request_jacobian,
                                   RealVectorX &f,
                                   RealMatrixX& jac_xdot,
                                   RealMatrixX& jac,
                                   MAST::BoundaryConditionBase& bc) {

    
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces

    // no parametric sensitivity is calculated for piston theory at this point.
    return (request_jacobian);
}



