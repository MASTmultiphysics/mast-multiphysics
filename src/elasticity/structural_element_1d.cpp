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
#include "elasticity/structural_element_1d.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/bending_operator.h"
#include "numerics/fem_operator_matrix.h"
#include "property_cards/element_property_card_1D.h"
#include "property_cards/material_property_card_base.h"
#include "base/system_initialization.h"
#include "base/boundary_condition_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/assembly_base.h"
#include "mesh/fe_base.h"
#include "mesh/geom_elem.h"


MAST::StructuralElement1D::StructuralElement1D(MAST::SystemInitialization& sys,
                                               MAST::AssemblyBase& assembly,
                                               const MAST::GeomElem& elem,
                                               const MAST::ElementPropertyCardBase& p):
MAST::BendingStructuralElem(sys, assembly, elem, p)  {
    
}



MAST::StructuralElement1D::~StructuralElement1D() {
    
}




void
MAST::StructuralElement1D::
initialize_direct_strain_operator(const unsigned int qp,
                                  const MAST::FEBase& fe,
                                  MAST::FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    
    unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi   = RealVectorX::Zero(n_phi);
    
    libmesh_assert_equal_to(Bmat.m(), 2);
    libmesh_assert_equal_to(Bmat.n(), 6*n_phi);
    libmesh_assert_less    (qp, dphi[0].size());
    
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
                                      const MAST::FEBase& fe,
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
    libmesh_assert_less    (qp, dphi[0].size());

    Real dv=0., dw=0.;
    vk_strain.setZero();
    vk_dvdxi_mat.setZero();
    vk_dwdxi_mat.setZero();
    
    RealVectorX phi_vec   = RealVectorX::Zero(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ ) {
        phi_vec(i_nd) = dphi[i_nd][qp](0);            // dphi/dx
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
                                                  const MAST::FEBase& fe,
                                                  RealMatrixX& vk_dvdxi_mat_sens,
                                                  RealMatrixX& vk_dwdxi_mat_sens) {
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    const unsigned int n_phi = (unsigned int)dphi.size();
    
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.rows(), 2);
    libmesh_assert_equal_to(vk_dvdxi_mat_sens.cols(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.rows(), 2);
    libmesh_assert_equal_to(vk_dwdxi_mat_sens.cols(), 2);
    libmesh_assert_less    (qp, dphi[0].size());

    Real dv=0., dw=0.;
    vk_dvdxi_mat_sens.setZero();
    vk_dwdxi_mat_sens.setZero();
    
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
                                            const MAST::FunctionBase* p,
                                            MAST::StressStrainOutputBase& output) {
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true, false));

    const unsigned int
    qp_loc_fe_size = (unsigned int)fe->get_qpoints().size(),
    n_added_qp     = 4;

    std::vector<libMesh::Point>
    qp_loc_fe = fe->get_qpoints(),
    qp_loc(qp_loc_fe_size*n_added_qp);
    
    
    // we will evaluate the stress at upper and lower layers of element,
    // so we will add two new points for each qp_loc
    // TODO: move this to element section property class for composite materials
    for (unsigned int i=0; i<qp_loc_fe.size(); i++) {
        
        qp_loc[i*4]        = qp_loc_fe[i];
        qp_loc[i*4](1)     = +1.;
        qp_loc[i*4](2)     = +1.; // upper right
        
        qp_loc[i*4+1]      = qp_loc_fe[i];
        qp_loc[i*4+1](1)   = -1.;
        qp_loc[i*4+1](2)   = +1.; // upper left

        qp_loc[i*4+2]      = qp_loc_fe[i];
        qp_loc[i*4+2](1)   = +1.;
        qp_loc[i*4+2](2)   = -1.; // lower right

        qp_loc[i*4+3]      = qp_loc_fe[i];
        qp_loc[i*4+3](1)   = -1.;
        qp_loc[i*4+3](2)   = -1.; // lower left
    }

    
    MAST::BendingOperatorType bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D>
    bend(MAST::build_bending_operator_1D(bending_model,
                                         *this,
                                         qp_loc_fe));

    
    // now that the FE object has been initialized, evaluate the stress values

    
    const std::vector<Real> &JxW              = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz    = fe->get_xyz();
    const unsigned int
    n_phi    = (unsigned int)fe->n_shape_functions(),
    n1       = this->n_direct_strain_components(),
    n2       = 6*n_phi,
    n3       = this->n_von_karman_strain_components();

    Real
    y     =  0.,
    z     =  0.,
    y_off =  0.,
    z_off =  0.,
    temp  =  0.,
    ref_t =  0.,
    alpha =  0.,
    dtemp =  0.,
    dref_t=  0.,
    dalpha=  0.;
    
    RealMatrixX
    material_mat,
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
    dstress_dp  = RealVectorX::Zero(n1),
    vec1        = RealVectorX::Zero(n2),
    vec2        = RealVectorX::Zero(n2);

    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1,  _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy

    // TODO: remove this const-cast, which may need change in API of
    // material card
    const MAST::FieldFunction<RealMatrixX >&
    mat_stiff  =
    const_cast<MAST::MaterialPropertyCardBase&>(_property.get_material()).stiffness_matrix(1);

    // get the thickness values for the bending strain calculation
    const MAST::FieldFunction<Real>
    &hy     =  _property.get<MAST::FieldFunction<Real> >("hy"),
    &hz     =  _property.get<MAST::FieldFunction<Real> >("hz"),
    &hy_off =  _property.get<MAST::FieldFunction<Real> >("hy_off"),
    &hz_off =  _property.get<MAST::FieldFunction<Real> >("hz_off");

    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN),
    if_bending = (_property.bending_model(_elem) != MAST::NO_BENDING);
    
    // a reference to the stress output data structure
    MAST::StressStrainOutputBase& stress_output =
    dynamic_cast<MAST::StressStrainOutputBase&>(output);

    // check to see if the element has any thermal loads specified
    // The object returns null
    MAST::BoundaryConditionBase *thermal_load =
    stress_output.get_thermal_load_for_elem(_elem);

    const MAST::FieldFunction<Real>
    *temp_func     = nullptr,
    *ref_temp_func = nullptr,
    *alpha_func    = nullptr;

    // get pointers to the temperature, if thermal load is specified
    if (thermal_load) {
        temp_func     =
        &(thermal_load->get<MAST::FieldFunction<Real> >("temperature"));
        ref_temp_func =
        &(thermal_load->get<MAST::FieldFunction<Real> >("ref_temperature"));
        alpha_func    =
        &(_property.get_material().get<MAST::FieldFunction<Real> >("alpha_expansion"));
    }
    
    // TODO: improve the stress calculation by including shear stress due
    // to torsion and shear flow due to beam bending.
    
    ///////////////////////////////////////////////////////////////////////
    // second for loop to calculate the residual and stiffness contributions
    ///////////////////////////////////////////////////////////////////////
    unsigned int
    qp = 0;
    for (unsigned int qp_loc_index=0; qp_loc_index<qp_loc_fe_size; qp_loc_index++)
        for (unsigned int section_qp_index=0; section_qp_index<n_added_qp; section_qp_index++)
        {
            qp = qp_loc_index*n_added_qp + section_qp_index;
            
            // get the material matrix
            mat_stiff(xyz[qp_loc_index], _time, material_mat);
            
            this->initialize_direct_strain_operator(qp_loc_index, *fe, Bmat_mem);
            
            // first handle constant throught the thickness stresses: membrane and vonKarman
            Bmat_mem.vector_mult(strain, _local_sol);
            
            // if thermal load was specified, then set the thermal strain
            // component of the total strain
            if (thermal_load) {
                (*temp_func)    (xyz[qp_loc_index], _time, temp);
                (*ref_temp_func)(xyz[qp_loc_index], _time, ref_t);
                (*alpha_func)   (xyz[qp_loc_index], _time, alpha);
                strain(0)  -=  alpha*(temp-ref_t);
            }
            
            
            if (if_bending) {
                // von Karman strain
                if (if_vk) {  // get the vonKarman strain operator if needed
                    
                    this->initialize_von_karman_strain_operator(qp_loc_index,
                                                                *fe,
                                                                strain_vk,
                                                                vk_dvdxi_mat,
                                                                vk_dwdxi_mat,
                                                                Bmat_v_vk,
                                                                Bmat_w_vk);
                    strain += strain_vk;
                }
                
                // add to this the bending strain
                hy    (xyz[qp_loc_index], _time,     y);
                hz    (xyz[qp_loc_index], _time,     z);
                hy_off(xyz[qp_loc_index], _time, y_off);
                hz_off(xyz[qp_loc_index], _time, z_off);
                
                // TODO: this assumes isotropic section. Multilayered sections need
                // special considerations
                // This operator depends on the y and z thickness values. Sensitivity
                // analysis should include the sensitivity of this operator on
                // these thickness values
                bend->initialize_bending_strain_operator_for_yz(*fe,
                                                                qp_loc_index,
                                                                qp_loc[qp](1) * y/2.+y_off,
                                                                qp_loc[qp](2) * z/2.+z_off,
                                                                Bmat_bend_v,
                                                                Bmat_bend_w);
                Bmat_bend_v.vector_mult(strain_bend, _local_sol);
                strain += strain_bend;
                
                Bmat_bend_w.vector_mult(strain_bend, _local_sol);
                strain += strain_bend;
            }
            
            
            // note that this assumes linear material laws
            stress = material_mat * strain;
            
            // now set the data for the 3D stress-strain vector
            // this is using only the direct strain/stress.
            // this can be improved by estimating the shear stresses from
            // torsion and shear flow from bending.
            strain_3D(0)  =   strain(0);
            stress_3D(0)  =   stress(0);
            
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
                                                                        xyz[qp_loc_index],
                                                                        stress_3D,
                                                                        strain_3D,
                                                                        JxW[qp_loc_index]));
            else
                data = &(stress_output.get_stress_strain_data_for_elem_at_qp(_elem,
                                                                             qp));
            
            // calculate the derivative if requested
            if (request_derivative || p) {
                
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
                    Bmat_bend_v.left_multiply(mat_n1n2, eye);
                    dstrain_dX  +=   mat_n1n2;

                    Bmat_bend_w.left_multiply(mat_n1n2, eye);
                    dstrain_dX  +=   mat_n1n2;
                }
                
                // note: this assumes linear material laws
                dstress_dX  = material_mat * dstrain_dX;
                
                // copy to the 3D structure
                vec1 = dstress_dX.row(0);
                this->transform_vector_to_global_system(vec1, vec2);
                dstress_dX_3D.row(0)  = vec2;
                vec1 = dstrain_dX.row(0);
                this->transform_vector_to_global_system(vec1, vec2);
                dstrain_dX_3D.row(0)  =  vec2;
                
                if (request_derivative)
                    data->set_derivatives(dstress_dX_3D, dstrain_dX_3D);
                
                
                if (p) {
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
                    
                    dstrain_dp  =  RealVectorX::Zero(n1);
                    
                    // if thermal load was specified, then set the thermal strain
                    // component of the total strain
                    if (thermal_load) {
                        temp_func->derivative(*p, xyz[qp_loc_index], _time, dtemp);
                        ref_temp_func->derivative(*p, xyz[qp_loc_index], _time, dref_t);
                        alpha_func->derivative(*p, xyz[qp_loc_index], _time, dalpha);
                        dstrain_dp(0)  -=  alpha*(dtemp-dref_t) + dalpha*(temp-ref_t);
                    }
                    
                    // include the dependence of strain on the thickness
                    if (if_bending) {
                        
                        // add to this the bending strain
                        hy.derivative(*p, xyz[qp_loc_index], _time, y);
                        hz.derivative(*p, xyz[qp_loc_index], _time, z);
                        hy_off.derivative(*p, xyz[qp_loc_index], _time, y_off);
                        hz_off.derivative(*p, xyz[qp_loc_index], _time, z_off);
                        
                        bend->initialize_bending_strain_operator_for_yz(*fe,
                                                                        qp_loc_index,
                                                                        qp_loc[qp](1) * y/2.+y_off,
                                                                        qp_loc[qp](2) * z/2.+z_off,
                                                                        Bmat_bend_v,
                                                                        Bmat_bend_w);
                        Bmat_bend_v.vector_mult(strain_bend, _local_sol);
                        dstrain_dp += strain_bend;
                        
                        Bmat_bend_w.vector_mult(strain_bend, _local_sol);
                        dstrain_dp += strain_bend;
                    }
                    
                    
                    // now use this to calculate the stress sensitivity.
                    dstress_dp  =  material_mat * dstrain_dp;
                    
                    // get the material matrix sensitivity
                    mat_stiff.derivative(*p, xyz[qp_loc_index], _time, material_mat);
                    
                    // partial sensitivity of strain is zero unless it is a
                    // shape parameter.
                    // TODO: shape sensitivity of strain operator
                    
                    // now use this to calculate the stress sensitivity.
                    dstress_dp +=  material_mat * strain;
                    
                    
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
                    data->set_sensitivity(*p,
                                          stress_3D,
                                          strain_3D);
                }
            }
        }
    
    // make sure that the number of data points for this element is
    // the same as the number of requested points
    libmesh_assert(qp_loc.size() ==
                   stress_output.n_stress_strain_data_for_elem(_elem));
    
    // if either derivative or sensitivity was requested, it was provided
    // by this routine
    return request_derivative || p;
}





bool
MAST::StructuralElement1D::internal_residual (bool request_jacobian,
                                              RealVectorX& f,
                                              RealMatrixX& jac)
{
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi    = (unsigned int)fe->get_phi().size(),
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
    
    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);

    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);

    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _property.stiffness_A_matrix(*this),
    mat_stiff_B  = _property.stiffness_B_matrix(*this),
    mat_stiff_D  = _property.stiffness_D_matrix(*this);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff_A)(xyz[qp], _time, material_A_mat);
        
        if (bend.get()) {
            (*mat_stiff_B)(xyz[qp], _time, material_B_mat);
            (*mat_stiff_D)(xyz[qp], _time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_residual_operation(if_vk, n2, qp, *fe, JxW,
                                     request_jacobian, 
                                     local_f, local_jac,
                                     bend.get(),
                                     Bmat_mem, Bmat_bend_v, Bmat_bend_w,
                                     Bmat_v_vk, Bmat_w_vk,
                                     stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                     material_A_mat,
                                     material_B_mat, material_D_mat, vec1_n1,
                                     vec2_n1, vec3_n2, vec4_n3,
                                     vec5_n3, mat1_n1n2, mat2_n2n2,
                                     mat3, mat4_n3n2);
        
    }
    
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (bend.get() && bend->include_transverse_shear_energy())
        bend->calculate_transverse_shear_residual(request_jacobian,
                                                  local_f,
                                                  local_jac);
    
    
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
MAST::StructuralElement1D::internal_residual_sensitivity (const MAST::FunctionBase& p,
                                                          bool request_jacobian,
                                                          RealVectorX& f,
                                                          RealMatrixX& jac)
{
    // this should be true if the function is called
    libmesh_assert(!p.is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(p);
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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

    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A = _property.stiffness_A_matrix(*this),
    mat_stiff_B = _property.stiffness_B_matrix(*this),
    mat_stiff_D = _property.stiffness_D_matrix(*this);
    
    // first calculate the sensitivity due to the parameter
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        mat_stiff_A->derivative(p, xyz[qp], _time, material_A_mat);
        
        if (bend.get()) {
            mat_stiff_B->derivative(p, xyz[qp], _time, material_B_mat);
            mat_stiff_D->derivative(p, xyz[qp], _time, material_D_mat);
        }
        
        // now calculte the quantity for these matrices
        _internal_residual_operation(if_vk, n2, qp, *fe, JxW,
                                     request_jacobian,
                                     local_f, local_jac,
                                     bend.get(),
                                     Bmat_mem, Bmat_bend_v, Bmat_bend_w,
                                     Bmat_v_vk, Bmat_w_vk,
                                     stress, stress_l, vk_dvdxi_mat, vk_dwdxi_mat,
                                     material_A_mat,
                                     material_B_mat, material_D_mat, vec1_n1,
                                     vec2_n1, vec3_n2, vec4_n3,
                                     vec5_n3, mat1_n1n2, mat2_n2n2,
                                     mat3, mat4_n3n2);
    }
    
    // now calculate the transverse shear contribution if appropriate for the
    // element
    if (bend.get() && bend->include_transverse_shear_energy())
        bend->calculate_transverse_shear_residual_sensitivity(p,
                                                              request_jacobian,
                                                              local_f,
                                                              local_jac);
    
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
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW            = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz  = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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
    local_jac         = RealMatrixX::Zero(n2,n2);
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    vec2_n1    = RealVectorX::Zero(n1),
    vec3_n2    = RealVectorX::Zero(n2),
    vec4_n3    = RealVectorX::Zero(n3),
    vec5_n3    = RealVectorX::Zero(n3);

    local_jac.setZero();


    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    // without the nonlinear strain, this matrix is zero.
    if (!if_vk)
        return false;
    
    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    mat_stiff_A  = _property.stiffness_A_matrix(*this),
    mat_stiff_B  = _property.stiffness_B_matrix(*this),
    mat_stiff_D  = _property.stiffness_D_matrix(*this);
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material matrix
        (*mat_stiff_A)(xyz[qp], _time, material_A_mat);
        
        (*mat_stiff_B)(xyz[qp], _time, material_B_mat);
        (*mat_stiff_D)(xyz[qp], _time, material_D_mat);
        
        // now calculte the quantity for these matrices
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // first handle constant throught the thickness stresses: membrane and vonKarman
        Bmat_mem.vector_mult(vec1_n1, _local_sol_sens);
        vec2_n1 = material_A_mat * vec1_n1; // linear direct stress
        
        // copy the stress values to a matrix
        stress(0,0)   = vec2_n1(0); // sigma_xx
        
        // get the bending strain operator
        vec2_n1.setZero(); // used to store vk strain, if applicable
        if (bend.get()) {
            bend->initialize_bending_strain_operator(*fe, qp,
                                                     Bmat_bend_v,
                                                     Bmat_bend_w);
            
            //  evaluate the bending stress and add that to the stress vector
            // for evaluation in the nonlinear stress term
            Bmat_bend_v.vector_mult(vec2_n1, _local_sol_sens);
            vec1_n1 = material_B_mat * vec2_n1;
            stress(0,0)   += vec1_n1(0);

            Bmat_bend_w.vector_mult(vec2_n1, _local_sol_sens);
            vec1_n1 = material_B_mat * vec2_n1;
            stress(0,0)   += vec1_n1(0);

            
            if (if_vk) {  // get the vonKarman strain operator if needed
                
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
                                                            vec2_n1,
                                                            vk_dvdxi_mat,
                                                            vk_dwdxi_mat,
                                                            Bmat_v_vk,
                                                            Bmat_w_vk);
                this->initialize_von_karman_strain_operator_sensitivity(qp,
                                                                        *fe,
                                                                        vk_dvdxi_mat_sens,
                                                                        vk_dwdxi_mat_sens);
                // sensitivity of von Karman strain
                vec2_n1.setZero();
                vec2_n1(0) = (vk_dvdxi_mat(0,0)*vk_dvdxi_mat_sens(0,0) +
                              vk_dwdxi_mat(0,0)*vk_dwdxi_mat_sens(0,0));
                vec1_n1 = material_A_mat * vec2_n1;
                stress(0,0) += vec1_n1(0);
            }
        }
        
        // copy the stress to use here.
        vec1_n1.setZero();
        
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
        Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;

        // bending - vk: w-displacement
        mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
        Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat_sens);
        mat3 = material_B_mat.transpose() * mat3;
        Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;

        // vk - bending: v-displacement
        Bmat_bend_v.left_multiply(mat1_n1n2, material_B_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;

        Bmat_bend_v.left_multiply(mat1_n1n2, material_B_mat);
        mat3 = vk_dwdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;

        // vk - bending: w-displacement
        Bmat_bend_w.left_multiply(mat1_n1n2, material_B_mat);
        mat3 = vk_dvdxi_mat_sens.transpose() * mat1_n1n2;
        Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
        local_jac += JxW[qp] * mat2_n2n2;
        
        Bmat_bend_w.left_multiply(mat1_n1n2, material_B_mat);
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
_internal_residual_operation(bool if_vk,
                             const unsigned int n2,
                             const unsigned int qp,
                             const MAST::FEBase& fe,
                             const std::vector<Real>& JxW,
                             bool request_jacobian,
                             RealVectorX& local_f,
                             RealMatrixX& local_jac,
                             MAST::BendingOperator1D* bend,
                             MAST::FEMOperatorMatrix& Bmat_mem,
                             MAST::FEMOperatorMatrix& Bmat_bend_v,
                             MAST::FEMOperatorMatrix& Bmat_bend_w,
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
    vec2_n1.setZero(); // used to store vk strain, if applicable
    if (bend) {
        bend->initialize_bending_strain_operator(fe, qp,
                                                 Bmat_bend_v,
                                                 Bmat_bend_w);
        
        //  evaluate the bending stress and add that to the stress vector
        // for evaluation in the nonlinear stress term
        Bmat_bend_v.vector_mult(vec2_n1, _local_sol);
        vec1_n1 = material_B_mat * vec2_n1;
        stress_l(0,0) += vec1_n1(0);
        stress(0,0)   += vec1_n1(0);

        Bmat_bend_w.vector_mult(vec2_n1, _local_sol);
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
    vec1_n1.setZero();
    vec1_n1(0)   = stress(0,0);
    vec1_n1(1)   = stress(0,1); // use the torsion strain from the temporary location
    stress(0, 1) = 0.;   // zero out the temporary value storing the torsion strain
    
    // now the internal force vector
    // this includes the membrane strain operator with all A and B material operators
    Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
    local_f += JxW[qp] * vec3_n2;
    
    if (bend) {
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
        vec1_n1 = material_B_mat.transpose() * vec2_n1;
        Bmat_bend_v.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        Bmat_bend_w.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;

        // now bending stress
        Bmat_bend_v.vector_mult(vec2_n1, _local_sol);
        vec1_n1 = material_D_mat * vec2_n1;
        Bmat_bend_v.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;

        Bmat_bend_w.vector_mult(vec2_n1, _local_sol);
        vec1_n1 = material_D_mat * vec2_n1;
        Bmat_bend_w.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
    }
    
    if (request_jacobian) {
        // membrane - membrane
        Bmat_mem.left_multiply(mat1_n1n2, material_A_mat);
        Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
        local_jac += JxW[qp] * mat2_n2n2;
                
        if (bend) {
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
                /*if (if_ignore_ho_jac) {
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
                else*/ {
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
                Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;

                // bending - vk: w-displacement
                mat3 = RealMatrixX::Zero(vk_dwdxi_mat.rows(), n2);
                Bmat_w_vk.left_multiply(mat3, vk_dwdxi_mat);
                mat3 = material_B_mat.transpose() * mat3;
                Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
                
                // vk - bending: v-displacement
                Bmat_bend_v.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;

                Bmat_bend_v.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;

                // vk - bending: w-displacement
                Bmat_bend_w.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dvdxi_mat.transpose() * mat1_n1n2;
                Bmat_v_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;

                Bmat_bend_w.left_multiply(mat1_n1n2, material_B_mat);
                mat3 = vk_dwdxi_mat.transpose() * mat1_n1n2;
                Bmat_w_vk.right_multiply_transpose(mat2_n2n2, mat3);
                local_jac += JxW[qp] * mat2_n2n2;
            }
            
            // bending - membrane
            mat3 = material_B_mat.transpose();
            Bmat_mem.left_multiply(mat1_n1n2, mat3);
            Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;
            Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;

            // membrane - bending
            Bmat_bend_v.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;

            Bmat_bend_w.left_multiply(mat1_n1n2, material_B_mat);
            Bmat_mem.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;

            // bending - bending
            Bmat_bend_v.left_multiply(mat1_n1n2, material_D_mat);
            Bmat_bend_v.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;

            Bmat_bend_w.left_multiply(mat1_n1n2, material_D_mat);
            Bmat_bend_w.right_multiply_transpose(mat2_n2n2, mat1_n1n2);
            local_jac += JxW[qp] * mat2_n2n2;
        }
    }
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

    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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
    
    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    prestress_A = _property.prestress_A_matrix(*this),
    prestress_B = _property.prestress_B_matrix(*this);
    
    // now calculate the quantity
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        (*prestress_A)(xyz[qp], _time, prestress_mat_A);
        (*prestress_B)(xyz[qp], _time, prestress_mat_B);
        _convert_prestress_A_mat_to_vector(prestress_mat_A, prestress_vec_A);
        _convert_prestress_B_mat_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.setZero(); // used to store vk strain, if applicable
        if (bend.get()) {
            bend->initialize_bending_strain_operator(*fe, qp,
                                                     Bmat_bend_v,
                                                     Bmat_bend_w);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
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
        
        if (bend.get()) {
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
            Bmat_bend_v.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
            Bmat_bend_w.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (bend.get() && if_vk) {
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
    return (request_jacobian);
}





bool
MAST::StructuralElement1D::prestress_residual_sensitivity (const MAST::FunctionBase& p,
                                                           bool request_jacobian,
                                                           RealVectorX& f,
                                                           RealMatrixX& jac)
{
    if (!_property.if_prestressed())
        return false;
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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

    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    std::unique_ptr<MAST::FieldFunction<RealMatrixX> >
    prestress_A = _property.prestress_A_matrix(*this),
    prestress_B = _property.prestress_B_matrix(*this);
    
    // transform to the local coordinate system
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        prestress_A->derivative(p, xyz[qp], _time, prestress_mat_A);
        prestress_B->derivative(p, xyz[qp], _time, prestress_mat_B);
        _convert_prestress_A_mat_to_vector(prestress_mat_A, prestress_vec_A);
        _convert_prestress_B_mat_to_vector(prestress_mat_B, prestress_vec_B);
        
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // get the bending strain operator if needed
        vec2_n1.setZero(); // used to store vk strain, if applicable
        if (bend.get()) {
            bend->initialize_bending_strain_operator(*fe, qp,
                                                     Bmat_bend_v,
                                                     Bmat_bend_w);
            
            if (if_vk)  // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
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
        
        if (bend.get()) {
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
            Bmat_bend_v.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
            Bmat_bend_w.vector_mult_transpose(vec3_n2, prestress_vec_B);
            local_f += JxW[qp] * vec3_n2; // epsilon_bend * sigma_0
        }
        
        if (request_jacobian) {
            if (bend.get() && if_vk) {
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
    return (request_jacobian);
}




bool
MAST::StructuralElement1D::
surface_pressure_residual(bool request_jacobian,
                          RealVectorX &f,
                          RealMatrixX &jac,
                          const unsigned int side,
                          MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(side, false));

    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const unsigned int
    n_phi  = (unsigned int)phi.size(),
    n1     = 3,
    n2     = 6*n_phi;
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& p_func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    
    const MAST::FieldFunction<Real>& A_func =
    dynamic_cast<const MAST::ElementPropertyCard1D&>(_property).A();
    
    FEMOperatorMatrix Bmat;
    Real
    press   = 0.,
    A_val   = 0.;
    
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
        
        // get pressure cross-sectional area value
        p_func(qpoint[qp], _time, press);
        A_func(qpoint[qp], _time, A_val);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = (press*A_val) * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    if (_elem.dim() < 3) {
        transform_vector_to_global_system(local_f, vec_n2);
        f -= vec_n2;
    }
    else
        f -= local_f;
    
    
    return (request_jacobian);
}





bool
MAST::StructuralElement1D::
surface_pressure_residual_sensitivity(const MAST::FunctionBase& p,
                                      bool request_jacobian,
                                      RealVectorX &f,
                                      RealMatrixX &jac,
                                      const unsigned int side,
                                      MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_elem.init_side_fe(side, false));

    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    const unsigned int
    n_phi  = (unsigned int)phi.size(),
    n1     = 3,
    n2     = 6*n_phi;
    
    
    // get the function from this boundary condition
    const MAST::FieldFunction<Real>& p_func =
    bc.get<MAST::FieldFunction<Real> >("pressure");
    const MAST::FieldFunction<Real>& A_func =
    dynamic_cast<const MAST::ElementPropertyCard1D&>(_property).A();

    
    FEMOperatorMatrix Bmat;
    Real
    press   =  0.,
    dpress  =  0.,
    A_val   =  0.,
    dA_val  =  0.;
    
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
        
        // get pressure and area values and their sensitivities
        p_func(qpoint[qp], _time, press);
        p_func.derivative(p, qpoint[qp], _time, dpress);
        A_func(qpoint[qp], _time, A_val);
        A_func.derivative(p, qpoint[qp], _time, dA_val);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = (press*dA_val + dpress*A_val) * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f += JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    if (_elem.dim() < 3) {
        transform_vector_to_global_system(local_f, vec_n2);
        f -= vec_n2;
    }
    else
        f -= local_f;
    
    
    return (request_jacobian);
}






bool
MAST::StructuralElement1D::thermal_residual (bool request_jacobian,
                                             RealVectorX& f,
                                             RealMatrixX& jac,
                                             MAST::BoundaryConditionBase& bc)
{
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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

    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    expansion_A = _property.thermal_expansion_A_matrix(*this),
    expansion_B = _property.thermal_expansion_B_matrix(*this);
    
    // temperature function
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real
    t       = 0.,
    t0      = 0.,
    scaling = 1.;
    
    if (bc.contains("thermal_jacobian_scaling"))
        bc.get<MAST::FieldFunction<Real>>("thermal_jacobian_scaling")(scaling);

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material property
        (*expansion_A)(xyz[qp], _time, material_exp_A_mat);
        (*expansion_B)(xyz[qp], _time, material_exp_B_mat);
        
        // get the temperature function
        temp_func    (xyz[qp], _time, t);
        ref_temp_func(xyz[qp], _time, t0);
        delta_t(0) = t-t0;
        
        vec1_n1 = material_exp_A_mat * delta_t; // [C]{alpha (T - T0)} (with membrane strain)
        stress(0,0) = vec1_n1(0); // sigma_xx
        vec2_n1 = material_exp_B_mat * delta_t; // [C]{alpha (T - T0)} (with bending strain)
        
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        
        if (bend.get()) {
            // bending strain
            bend->initialize_bending_strain_operator(*fe, qp,
                                                     Bmat_bend_v,
                                                     Bmat_bend_w);
            Bmat_bend_v.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;
            Bmat_bend_w.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;

            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
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
        jac -= scaling * mat2_n2n2;
    }
    
    // Jacobian contribution from von Karman strain
    return request_jacobian;
}




bool
MAST::StructuralElement1D::thermal_residual_sensitivity (const MAST::FunctionBase& p,
                                                         bool request_jacobian,
                                                         RealVectorX& f,
                                                         RealMatrixX& jac,
                                                         MAST::BoundaryConditionBase& bc)
{
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true,
                                                     false,
                                                     _property.extra_quadrature_order(_elem)));

    const std::vector<Real>& JxW           = fe->get_JxW();
    const std::vector<libMesh::Point>& xyz = fe->get_xyz();
    const unsigned int
    n_phi = (unsigned int)fe->get_phi().size(),
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
    
    MAST::FEMOperatorMatrix
    Bmat_mem,
    Bmat_bend_v,
    Bmat_bend_w,
    Bmat_v_vk,
    Bmat_w_vk;
    
    Bmat_mem.reinit(n1, _system.n_vars(), n_phi); // three stress-strain components
    Bmat_bend_v.reinit(n1, _system.n_vars(), n_phi);
    Bmat_bend_w.reinit(n1, _system.n_vars(), n_phi);
    Bmat_v_vk.reinit(n3, _system.n_vars(), n_phi); // only dv/dx and dv/dy
    Bmat_w_vk.reinit(n3, _system.n_vars(), n_phi); // only dw/dx and dw/dy
    
    bool if_vk = (_property.strain_type() == MAST::NONLINEAR_STRAIN);
    
    MAST::BendingOperatorType
    bending_model = _property.bending_model(_elem);
    
    std::unique_ptr<MAST::BendingOperator1D> bend;
    
    if (bending_model != MAST::NO_BENDING)
        bend.reset(MAST::build_bending_operator_1D(bending_model,
                                                   *this,
                                                   fe->get_qpoints()).release());

    std::unique_ptr<MAST::FieldFunction<RealMatrixX > >
    expansion_A = _property.thermal_expansion_A_matrix(*this),
    expansion_B = _property.thermal_expansion_B_matrix(*this);
    
    // temperature function
    const MAST::FieldFunction<Real>
    &temp_func     = bc.get<MAST::FieldFunction<Real> >("temperature"),
    &ref_temp_func = bc.get<MAST::FieldFunction<Real> >("ref_temperature");
    
    Real t, t0, t_sens, t0_sens;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // get the material property
        (*expansion_A)(xyz[qp], _time, material_exp_A_mat);
        (*expansion_B)(xyz[qp], _time, material_exp_B_mat);
        expansion_A->derivative(p, xyz[qp], _time, material_exp_A_mat_sens);
        expansion_B->derivative(p, xyz[qp], _time, material_exp_B_mat_sens);
        
        // get the temperature function
        temp_func(xyz[qp], _time, t);
        temp_func.derivative(p, xyz[qp], _time, t_sens);
        ref_temp_func(xyz[qp], _time, t0);
        ref_temp_func.derivative(p, xyz[qp], _time, t0_sens);
        delta_t(0)      = t-t0;
        delta_t_sens(0) = t_sens-t0_sens;
        
        // now prepare the membrane force sensitivity
        vec1_n1 = material_exp_A_mat * delta_t_sens; // [C]{alpha (dT/dp)} (with membrane strain)
        vec2_n1 = material_exp_A_mat_sens * delta_t; // d([C].{alpha})/dp (T - T0)} (with membrane
        vec1_n1 += vec2_n1;  // sensitivity of the thermal membrane force
        stress(0,0) = vec1_n1(0); // sigma_xx
        
        // now prepare the membrane-bending coupling force sensitivity
        vec2_n1 = material_exp_B_mat * delta_t_sens; // [C]{alpha dT/dp} (with bending strain)
        vec5_n1 = material_exp_B_mat_sens * delta_t; // d([C].{alpha})/dp (T - T0)} (with bending strain)
        vec2_n1 += vec5_n1;  // sensitivity of the thermal membrane force
        
        
        this->initialize_direct_strain_operator(qp, *fe, Bmat_mem);
        
        // membrane strain
        Bmat_mem.vector_mult_transpose(vec3_n2, vec1_n1);
        local_f += JxW[qp] * vec3_n2;
        
        if (bend.get()) {
            // bending strain
            bend->initialize_bending_strain_operator(*fe, qp,
                                                     Bmat_bend_v,
                                                     Bmat_bend_w);
            Bmat_bend_v.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;
            Bmat_bend_w.vector_mult_transpose(vec3_n2, vec2_n1);
            local_f += JxW[qp] * vec3_n2;

            // von Karman strain
            if (if_vk) {
                // get the vonKarman strain operator if needed
                this->initialize_von_karman_strain_operator(qp,
                                                            *fe,
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
    return request_jacobian;
}







bool
MAST::StructuralElement1D::
piston_theory_residual(bool request_jacobian,
                       RealVectorX &f,
                       RealMatrixX& jac_xdot,
                       RealMatrixX& jac,
                       MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true, false));

    const std::vector<Real> &JxW                = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi;
    
    
    
    // convert to piston theory boundary condition so that the necessary
    // flow properties can be obtained
    const MAST::PistonTheoryBoundaryCondition& piston_bc =
    dynamic_cast<MAST::PistonTheoryBoundaryCondition&>(bc);
    
    // create the constant field functions to pass the dwdx and dwdt values
    // to the piston theory pressure functions
    MAST::Parameter
    dwdx_p  ("dwdx", 0.),
    dwdt_p  ("dwdt", 0.);
    
    MAST::ConstantFieldFunction
    dwdx_f  ("dwdx", dwdx_p),
    dwdt_f  ("dwdx", dwdt_p);
    
    std::unique_ptr<MAST::FieldFunction<Real> >
    pressure        (piston_bc.get_pressure_function(dwdx_f, dwdt_f).release()),
    dpressure_dx    (piston_bc.get_dpdx_function    (dwdx_f, dwdt_f).release()),
    dpressure_dxdot (piston_bc.get_dpdxdot_function (dwdx_f, dwdt_f).release());

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
    p_val    = RealVectorX::Zero(3),
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
                                              *fe,
                                              dummy,
                                              dvdx,
                                              dwdx,
                                              Bmat_dvdx,
                                              Bmat_dwdx);
        
        
        // now use this information to evaluate the normal cp due to
        // both the v and w motions
        dwdx_p = dvdx(0,0);
        dwdt_p = vel(1);
        (*pressure)(qpoint[qp], _time, p_val(1));

        
        dwdx_p = dwdx(0,0);
        dwdt_p = vel(2);
        (*pressure)(qpoint[qp], _time, p_val(2));

        
        // calculate force and discrete vector for v-disp
        force.setZero();
        force(0) = p_val(1) * normal(1);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;

        // now the w-disp
        force.setZero();
        force(1) = p_val(2) * normal(2);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;
        
        
        // calculate the Jacobian if requested
        if (request_jacobian) {
            
            // we need the derivative of cp wrt normal velocity
            dwdx_p = dvdx(0,0);
            dwdt_p = vel(1);
            (*dpressure_dxdot)(qpoint[qp], _time, p_val(1));
            
            // calculate the component of Jacobian due to v-velocity
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_v);
            local_jac_xdot += JxW[qp] * p_val(1) * normal(1) * mat_n2n2;

            // now wrt v
            (*dpressure_dx)(qpoint[qp], _time, p_val(1));
            // now use calculate the component of Jacobian due to deformation
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_dvdx);  // v: B^T dB/dx
            local_jac += (JxW[qp] * p_val(1) * normal(1)) * mat_n2n2;
            

            dwdx_p = dwdx(0,0);
            dwdt_p = vel(2);
            (*dpressure_dxdot)(qpoint[qp], _time, p_val(2));

            // calculate the component of Jacobian due to w-velocity
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_w);
            local_jac_xdot += JxW[qp] * p_val(2) * normal(2) * mat_n2n2;
            
            // now wrt w
            (*dpressure_dx)(qpoint[qp], _time, p_val(2));
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_dwdx);  // w: B^T dB/dx
            local_jac += (JxW[qp] * p_val(2) * normal(2)) * mat_n2n2;
        }
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
    
    // if the Jacobian was requested, then transform it and add to the
    // global Jacobian
    if (request_jacobian) {
        transform_matrix_to_global_system(local_jac_xdot, mat_n2n2);
        jac_xdot -= mat_n2n2;

        transform_matrix_to_global_system(local_jac, mat_n2n2);
        jac      -= mat_n2n2;
    }
    
    return (request_jacobian);
}










bool
MAST::StructuralElement1D::
piston_theory_residual_sensitivity(const MAST::FunctionBase& p,
                                   bool request_jacobian,
                                   RealVectorX &f,
                                   RealMatrixX& jac_xdot,
                                   RealMatrixX& jac,
                                   MAST::BoundaryConditionBase& bc) {

    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    std::unique_ptr<MAST::FEBase>   fe(_elem.init_fe(true, false));

    const std::vector<Real> &JxW                = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint   = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi  = fe->get_phi();
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = this->n_direct_strain_components(),
    n2    = 6*n_phi;
    
    
    
    // convert to piston theory boundary condition so that the necessary
    // flow properties can be obtained
    const MAST::PistonTheoryBoundaryCondition& piston_bc =
    dynamic_cast<MAST::PistonTheoryBoundaryCondition&>(bc);
    
    // create the constant field functions to pass the dwdx and dwdt values
    // to the piston theory pressure functions
    MAST::Parameter
    dwdx_p  ("dwdx", 0.),
    dwdt_p  ("dwdt", 0.);
    
    MAST::ConstantFieldFunction
    dwdx_f  ("dwdx", dwdx_p),
    dwdt_f  ("dwdx", dwdt_p);
    
    std::unique_ptr<MAST::FieldFunction<Real> >
    pressure        (piston_bc.get_pressure_function(dwdx_f, dwdt_f).release()),
    dpressure_dx    (piston_bc.get_dpdx_function    (dwdx_f, dwdt_f).release()),
    dpressure_dxdot (piston_bc.get_dpdxdot_function (dwdx_f, dwdt_f).release());
    
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
    p_val    = RealVectorX::Zero(3),
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
                                              *fe,
                                              dummy,
                                              dvdx,
                                              dwdx,
                                              Bmat_dvdx,
                                              Bmat_dwdx);
        
        
        // now use this information to evaluate the normal cp due to
        // both the v and w motions
        dwdx_p = dvdx(0,0);
        dwdt_p = vel(1);
        pressure->derivative(p, qpoint[qp], _time, p_val(1));
        
        
        dwdx_p = dwdx(0,0);
        dwdt_p = vel(2);
        pressure->derivative(p, qpoint[qp], _time, p_val(2));
        
        
        // calculate force and discrete vector for v-disp
        force.setZero();
        force(0) = p_val(1) * normal(1);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;
        
        // now the w-disp
        force.setZero();
        force(1) = p_val(2) * normal(2);
        Bmat_v.vector_mult_transpose(vec_n2, force);
        local_f += JxW[qp] * vec_n2;
        
        
        // calculate the Jacobian if requested
        if (request_jacobian) {
            
            // we need the derivative of cp wrt normal velocity
            dwdx_p = dvdx(0,0);
            dwdt_p = vel(1);
            dpressure_dxdot->derivative(p, qpoint[qp], _time, p_val(1));
            
            // calculate the component of Jacobian due to v-velocity
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_v);
            local_jac_xdot += JxW[qp] * p_val(1) * normal(1) * mat_n2n2;
            
            // now wrt v
            dpressure_dx->derivative(p, qpoint[qp], _time, p_val(1));
            
            // now use calculate the component of Jacobian due to deformation
            Bmat_v.right_multiply_transpose(mat_n2n2, Bmat_dvdx);  // v: B^T dB/dx
            local_jac += (JxW[qp] * p_val(1) * normal(1)) * mat_n2n2;
            
            
            dwdx_p = dwdx(0,0);
            dwdt_p = vel(2);
            dpressure_dxdot->derivative(p, qpoint[qp], _time, p_val(2));
            
            // calculate the component of Jacobian due to w-velocity
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_w);
            local_jac_xdot += JxW[qp] * p_val(2) * normal(2) * mat_n2n2;
            
            // now wrt w
            dpressure_dx->derivative(p, qpoint[qp], _time, p_val(2));
            Bmat_w.right_multiply_transpose(mat_n2n2, Bmat_dwdx);  // w: B^T dB/dx
            local_jac += (JxW[qp] * p_val(2) * normal(2)) * mat_n2n2;
        }
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f -= vec_n2;
    
    // if the Jacobian was requested, then transform it and add to the
    // global Jacobian
    if (request_jacobian) {
        transform_matrix_to_global_system(local_jac_xdot, mat_n2n2);
        jac_xdot -= mat_n2n2;
        
        transform_matrix_to_global_system(local_jac, mat_n2n2);
        jac      -= mat_n2n2;
    }
    
    
    // no parametric sensitivity is calculated for piston theory at this point.
    return (request_jacobian);
}





/*bool
MAST::StructuralElement1D::
linearized_frequency_domain_surface_pressure_residual
 (bool request_jacobian,
 ComplexVectorX &f,
 ComplexMatrixX &jac,
 const unsigned int side,
 MAST::BoundaryConditionBase& bc) {
    
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(bc.type(), MAST::SMALL_DISTURBANCE_MOTION);
    
    
    MAST::FieldFunction<Real>&
    press_fn  = bc.get<MAST::FieldFunction<Real> >          ("pressure");
    MAST::FieldFunction<Complex>&
    dpress_fn = bc.get<MAST::FieldFunction<Complex> >       ("dpressure");
    MAST::FieldFunction<ComplexVectorX>&
    dn_rot_fn = bc.get<MAST::FieldFunction<ComplexVectorX> >("dnormal");
    
    
    libMesh::FEBase *fe_ptr    = nullptr;
    libMesh::QBase  *qrule_ptr = nullptr;
    _get_side_fe_and_qrule(get_elem_for_quadrature(),
                           side,
                           &fe_ptr,
                           &qrule_ptr,
                           false);
    std::unique_ptr<libMesh::FEBase> fe(fe_ptr);
    std::unique_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    
    // Physical location of the quadrature points
    const std::vector<Real> &JxW                    = fe->get_JxW();
    const std::vector<libMesh::Point>& qpoint       = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi      = fe->get_phi();
    const std::vector<libMesh::Point>& face_normals = fe->get_normals_for_local_coordinate();
    
    const unsigned int
    n_phi = (unsigned int)phi.size(),
    n1    = 3,
    n2    = 6*n_phi;
    
    RealVectorX phi_vec   = RealVectorX::Zero(n_phi);
    ComplexVectorX
    dn_rot  = ComplexVectorX::Zero(3),
    force   = ComplexVectorX::Zero(2*n1),
    local_f = ComplexVectorX::Zero(n2),
    vec_n2  = ComplexVectorX::Zero(n2);
    
    
    FEMOperatorMatrix Bmat;
    Real              press;
    Complex           dpress;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++) {
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        press_fn (qpoint[qp], _time,  press);
        dpress_fn(qpoint[qp], _time, dpress);
        dn_rot_fn(qpoint[qp], _time, dn_rot);
        
        //            press = 0.;
        //            dpress = Complex(2./4.*std::real(dn_rot(0)),  2./4./.1*std::imag(utrans(1)));
        //            libMesh::out << q_point[qp](0)
        //            << "  " << std::real(utrans(1))
        //            << "  " << std::imag(utrans(1))
        //            << "  " << std::real(dn_rot(0))
        //            << "  " << std::imag(dn_rot(0))
        //            << "  " << std::real(press)
        //            << "  " << std::imag(press)
        //            << "  " << std::real(dpress)
        //            << "  " << std::imag(dpress) << std::endl;
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) =  ( press * dn_rot(i_dim) + // steady pressure
                             dpress * face_normals[qp](i_dim)); // unsteady pressure
        
        
        Bmat.vector_mult_transpose(vec_n2, force);
        
        local_f -= JxW[qp] * vec_n2;
    }
    
    // now transform to the global system and add
    transform_vector_to_global_system(local_f, vec_n2);
    f  += vec_n2;
    
    
    return (request_jacobian);
}
*/


