/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "fluid/conservative_fluid_element_base.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/flight_condition.h"
#include "base/boundary_condition_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"


MAST::ConservativeFluidElementBase::
ConservativeFluidElementBase(MAST::SystemInitialization& sys,
                             const libMesh::Elem& elem,
                             const MAST::FlightCondition& f):
MAST::FluidElemBase(elem.dim(), f),
MAST::ElementBase(sys, elem) {
    
    // initialize the finite element data structures
    _init_fe_and_qrule(elem, &_fe, &_qrule);
}



MAST::ConservativeFluidElementBase::~ConservativeFluidElementBase() {
    
}




bool
MAST::ConservativeFluidElementBase::internal_residual (bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac) {
    const std::vector<Real>& JxW                  = _fe->get_JxW();
    const std::vector<std::vector<Real> >& phi    = _fe->get_phi();
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1,
    nphi   = _fe->n_shape_functions();
    
    RealMatrixX
    mat1_n1n1       = RealMatrixX::Zero(   n1,    n1),
    mat2_n1n1       = RealMatrixX::Zero(   n1,    n1),
    mat3_n1n2       = RealMatrixX::Zero(   n1,    n2),
    mat4_n2n2       = RealMatrixX::Zero(   n2,    n2),
    AiBi_adv        = RealMatrixX::Zero(   n1,    n2),
    A_sens          = RealMatrixX::Zero(   n1,    n2),
    LS              = RealMatrixX::Zero(   n1,    n2),
    LS_sens         = RealMatrixX::Zero(   n2,    n2);
    
    RealVectorX
    vec1_n1  = RealVectorX::Zero(n1),
    vec2_n1  = RealVectorX::Zero(n1),
    vec3_n2  = RealVectorX::Zero(n2),
    dc       = RealVectorX::Zero(dim);
    
    std::vector<RealMatrixX>
    Ai_adv  (dim);
    
    std::vector<std::vector<RealMatrixX> >
    Ai_sens  (dim);
    
    
    for (unsigned int i=0; i<dim; i++) {
        Ai_sens [i].resize(n1);
        Ai_adv  [i].setZero(n1, n1);
        for (unsigned int j=0; j<n1; j++)
            Ai_sens[i][j].setZero(n1, n1);
    }
    
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix Bmat;
    MAST::PrimitiveSolution      primitive_sol;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *_fe, Bmat);
        
        // calculate the local element solution
        Bmat.right_multiply(vec1_n1, _sol);
        
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        // initialize the FEM derivative operator
        _initialize_fem_gradient_operator(qp, dim, *_fe, dBmat);
        
        AiBi_adv.setZero();
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            calculate_advection_flux_jacobian(i_dim, primitive_sol, Ai_adv[i_dim]);
            calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
            (i_dim, primitive_sol, Ai_sens[i_dim]);
            
            dBmat[i_dim].left_multiply(mat3_n1n2, Ai_adv[i_dim]);
            AiBi_adv += mat3_n1n2;
        }
        
        // intrinsic time operator for this quadrature point
        calculate_differential_operator_matrix(qp,
                                               *_fe,
                                               _sol,
                                               primitive_sol,
                                               Bmat,
                                               dBmat,
                                               Ai_adv,
                                               AiBi_adv,
                                               Ai_sens,
                                               LS,
                                               LS_sens);
        
        // discontinuity capturing operator for this quadrature point
        calculate_aliabadi_discontinuity_operator(qp,
                                                  *_fe,
                                                  primitive_sol,
                                                  _sol,
                                                  dBmat,
                                                  AiBi_adv,
                                                  dc);
        
        // assemble the residual due to flux operator
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            
            // first the flux
            calculate_advection_flux(i_dim, primitive_sol, vec1_n1);
            dBmat[i_dim].vector_mult_transpose(vec3_n2, vec1_n1);
            f -= JxW[qp] * vec3_n2;
            
            // solution derivative in i^th direction
            // use this to calculate the discontinuity capturing term
            dBmat[i_dim].vector_mult(vec1_n1, _sol);
            dBmat[i_dim].vector_mult_transpose(vec3_n2, vec1_n1);
            f += JxW[qp] * dc(i_dim) * vec3_n2;
        }
        
        // stabilization term
        f += JxW[qp] * LS.transpose() * (AiBi_adv * _sol);
        
        
        if (request_jacobian) {
            
            A_sens.setZero();
            
            // contribution from flux Jacobian
            for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
                // flux term
                Bmat.left_multiply(mat3_n1n2, Ai_adv[i_dim]);                        // A_i B
                dBmat[i_dim].right_multiply_transpose(mat4_n2n2, mat3_n1n2);          // dB_i^T A_i B
                jac -= JxW[qp]*mat4_n2n2;
                
                // sensitivity of Ai_Bi with respect to U:   [dAi/dUj.Bi.U  ...  dAi/dUn.Bi.U]
                dBmat[i_dim].vector_mult(vec1_n1, _sol);
                for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++) {
                    
                    vec2_n1 = Ai_sens[i_dim][i_cvar] * vec1_n1;
                    for (unsigned int i_phi=0; i_phi<nphi; i_phi++)
                        A_sens.col(nphi*i_cvar+i_phi) += phi[i_phi][qp] *vec2_n1; // assuming that all variables have same n_phi
                }
                
                // discontinuity capturing term
                dBmat[i_dim].right_multiply_transpose(mat4_n2n2, dBmat[i_dim]);   // dB_i^T dc dB_i
                jac += JxW[qp] * dc(i_dim) * mat4_n2n2;
            }
            
            // stabilization term
            jac  += JxW[qp] * LS.transpose() * AiBi_adv;                          // A_i dB_i

            // linearization of the Jacobian terms
            jac += JxW[qp] * LS.transpose() * A_sens; // LS^T tau d^2F^adv_i / dx dU  (Ai sensitivity)
                                          // linearization of the LS terms
            jac += JxW[qp] * LS_sens;
            
        }
    }
    
    return request_jacobian;
}





bool
MAST::ConservativeFluidElementBase::velocity_residual (bool request_jacobian,
                                                       RealVectorX& f,
                                                       RealMatrixX& jac_xdot,
                                                       RealMatrixX& jac) {
    const std::vector<Real>& JxW           = _fe->get_JxW();
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealMatrixX
    tau              = RealMatrixX::Zero(n1, n1),
    mat1_n1n1        = RealMatrixX::Zero(n1, n1),
    mat2_n1n2        = RealMatrixX::Zero(n1, n2),
    mat3_n2n2        = RealMatrixX::Zero(n2, n2),
    mat4_n2n1        = RealMatrixX::Zero(n2, n1),
    AiBi_adv         = RealMatrixX::Zero(n1, n2),
    A_sens           = RealMatrixX::Zero(n1, n2),
    LS               = RealMatrixX::Zero(n1, n2),
    LS_sens          = RealMatrixX::Zero(n2, n2);
    RealVectorX
    vec1_n1          = RealVectorX::Zero(n1),
    vec2_n1          = RealVectorX::Zero(n1),
    vec3_n2          = RealVectorX::Zero(n2);
    
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    MAST::FEMOperatorMatrix      Bmat;
    MAST::PrimitiveSolution      primitive_sol;
    
    std::vector<RealMatrixX>
    Ai_adv  (dim);
    
    std::vector<std::vector<RealMatrixX> >
    Ai_sens  (dim);
    
    
    for (unsigned int i=0; i<dim; i++) {
        Ai_sens [i].resize(n1);
        Ai_adv  [i].setZero(n1, n1);
        for (unsigned int j=0; j<n1; j++)
            Ai_sens[i][j].setZero(n1, n1);
    }
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // first need to set the solution of the conservative operator
        _initialize_fem_interpolation_operator(qp, dim, *_fe, Bmat);
        Bmat.right_multiply(vec1_n1, _sol);                                     //  B * U
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        // initialize the FEM derivative operator
        _initialize_fem_gradient_operator(qp, dim, *_fe, dBmat);
        
        AiBi_adv.setZero();
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            calculate_advection_flux_jacobian(i_dim, primitive_sol, Ai_adv[i_dim]);
            
            dBmat[i_dim].left_multiply(mat2_n1n2, Ai_adv[i_dim]);
            AiBi_adv += mat2_n1n2;
        }
        
        // intrinsic time operator for this quadrature point
        calculate_differential_operator_matrix(qp,
                                               *_fe,
                                               _sol,
                                               primitive_sol,
                                               Bmat,
                                               dBmat,
                                               Ai_adv,
                                               AiBi_adv,
                                               Ai_sens,
                                               LS,
                                               LS_sens);
        
        // now evaluate the Jacobian due to the velocity term
        Bmat.right_multiply(vec1_n1, _vel);                                     //  B * U_dot
        Bmat.vector_mult_transpose(vec3_n2, vec1_n1);                           //  B^T * B * U_dot
        f += JxW[qp] * vec3_n2;
        
        // next, evaluate the contribution from the stabilization term
        f += JxW[qp] * LS.transpose() * vec1_n1;
        
        if (request_jacobian) {
            
            // contribution from the velocity term
            Bmat.right_multiply_transpose(mat3_n2n2, Bmat);
            jac_xdot += JxW[qp] * mat3_n2n2;                 // B^T B
            
            // contribution from the stabilization term
            // next, evaluate the contribution from the stabilization term
            mat4_n2n1   = LS.transpose();
            Bmat.left_multiply(mat3_n2n2, mat4_n2n1);     // LS^T B
            jac_xdot += JxW[qp]*mat3_n2n2;
        }
    }
    
    
    return request_jacobian;
}




bool
MAST::ConservativeFluidElementBase::
side_external_residual (bool request_jacobian,
                        RealVectorX& f,
                        RealMatrixX& jac,
                        std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SYMMETRY_WALL:
                        calculate_jac = (calculate_jac ||
                                         symmetry_surface_residual(request_jacobian,
                                                                   f, jac,
                                                                   n,
                                                                   *it.first->second));
                        break;
                        
                    case MAST::SLIP_WALL:
                        calculate_jac = (calculate_jac ||
                                         slip_wall_surface_residual(request_jacobian,
                                                                    f, jac,
                                                                    n,
                                                                    *it.first->second));
                        break;
                        
                    case MAST::FAR_FIELD:
                        calculate_jac = (calculate_jac ||
                                         far_field_surface_residual(request_jacobian,
                                                                    f, jac,
                                                                    n,
                                                                    *it.first->second));
                        break;
                        
                    case MAST::DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}





bool
MAST::ConservativeFluidElementBase::
side_external_residual_sensitivity (bool request_jacobian,
                                    RealVectorX& f,
                                    RealMatrixX& jac,
                                    std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SYMMETRY_WALL:
                        calculate_jac = (calculate_jac ||
                                         symmetry_surface_residual_sensitivity(request_jacobian,
                                                                               f, jac,
                                                                               n,
                                                                               *it.first->second));
                        break;
                        
                    case MAST::SLIP_WALL:
                        calculate_jac = (calculate_jac ||
                                         slip_wall_surface_residual_sensitivity(request_jacobian,
                                                                                f, jac,
                                                                                n,
                                                                                *it.first->second));
                        break;
                        
                    case MAST::FAR_FIELD:
                        calculate_jac = (calculate_jac ||
                                         far_field_surface_residual_sensitivity(request_jacobian,
                                                                                f, jac,
                                                                                n,
                                                                                *it.first->second));
                        break;
                        
                    case MAST::DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}






bool
MAST::ConservativeFluidElementBase::
internal_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}



bool
MAST::ConservativeFluidElementBase::
velocity_residual_sensitivity (bool request_jacobian,
                               RealVectorX& f,
                               RealMatrixX& jac) {
    
    return request_jacobian;
}








bool
MAST::ConservativeFluidElementBase::
symmetry_surface_residual(bool request_jacobian,
                          RealVectorX& f,
                          RealMatrixX& jac,
                          const unsigned int s,
                          MAST::BoundaryConditionBase& p) {
    
    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = NULL;
    libMesh::QBase  *qrule_ptr = NULL;
    _get_side_fe_and_qrule(_elem, s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    dnormal   = RealVectorX::Zero(dim);
    
    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1, _sol);
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        
        vec1_n1.setZero();
        // since vi=0, vi ni = 0, so the advection flux gets evaluated
        // using the slip wall condition
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            vec1_n1(i_dim+1) += primitive_sol.p * normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
        f += JxW[qp] * vec2_n2;
        
        if (request_jacobian) {
            
            calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             0.,
             normals[qp],
             dnormal,
             mat1_n1n1);
            
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2);
            jac += JxW[qp] * mat3_n2n2;
        }
    }
    
    // calculation of the load vector is independent of solution
    return false;
}








bool
MAST::ConservativeFluidElementBase::
symmetry_surface_residual_sensitivity(bool request_jacobian,
                                      RealVectorX& f,
                                      RealMatrixX& jac,
                                      const unsigned int s,
                                      MAST::BoundaryConditionBase& p) {
    
    return false;
}





bool
MAST::ConservativeFluidElementBase::
slip_wall_surface_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           const unsigned int s,
                           MAST::BoundaryConditionBase& p) {
    
    // inviscid boundary condition without any diffusive component
    // conditions enforced are
    // vi ni = wi_dot (ni + dni) - ui dni   (moving slip wall with deformation)
    // tau_ij nj = 0   (because velocity gradient at wall = 0)
    // qi ni = 0       (since heat flux occurs only on no-slip wall and far-field bc)
    
    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = NULL;
    libMesh::QBase  *qrule_ptr = NULL;
    _get_side_fe_and_qrule(_elem, s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    flux      = RealVectorX::Zero(n1),
    dnormal   = RealVectorX::Zero(dim);
    
    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    
    
    Real ui_ni = 0.;
    
    for (unsigned int qp=0; qp<JxW.size(); qp++)
    {
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1, _sol);
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        //        // copy the surface normal
        //        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        //            local_normal(i_dim) = face_normals[qp](i_dim);
        //
        //        // now check if the surface deformation is defined and
        //        // needs to be applied through transpiration boundary
        //        // condition
        //        ui_ni = 0.;
        //        primitive_sol.get_uvec(uvec);
        //
        //        if (surface_motion) // get the surface motion data
        //        {
        //            surface_motion->surface_velocity(this->time,
        //                                             qpoint[qp],
        //                                             face_normals[qp],
        //                                             surface_def,
        //                                             surface_vel,
        //                                             dnormal);
        //
        //            // update the normal with the deformation
        //            // this defines the normal of the surface that has been
        //            // deformed, although the geometry of the flow mesh does
        //            // not conform to that deformation
        //            local_normal.add(1., dnormal);
        //
        //            //    ui * (ni + dni) = wi_dot * (ni + dni)
        //            // => ui * ni = wi_dot * (ni + dni) - ui * dni
        //
        //            // now initialize the surface velocity
        //            // note that the perturbed local normal is used here
        //            // since it resembles the normal of the surface that
        //            // has undergone a static deformation, even though the
        //            // surface_vel might be identically zero for a static
        //            // body.
        //            ui_ni  = surface_vel.dot(local_normal);
        //            ui_ni -= uvec.dot(dnormal);
        //        }
        
        flux.setZero();
        flux += ui_ni * vec1_n1;
        flux(n1-1) += primitive_sol.p*ui_ni;
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            flux(i_dim+1) += primitive_sol.p * normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;
        
        if ( request_jacobian ) {
            
            this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             ui_ni,
             normals[qp],
             dnormal,
             mat1_n1n1);
            
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2);
            jac += JxW[qp] * mat3_n2n2;
        }
    }

    
    return false;
}




bool
MAST::ConservativeFluidElementBase::
slip_wall_surface_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       const unsigned int s,
                                       MAST::BoundaryConditionBase& p) {
    
    return false;
}



bool
MAST::ConservativeFluidElementBase::
far_field_surface_residual(bool request_jacobian,
                           RealVectorX& f,
                           RealMatrixX& jac,
                           const unsigned int s,
                           MAST::BoundaryConditionBase& p) {
    
    // conditions enforced are:
    // -- f_adv_i ni =  f_adv = f_adv(+) + f_adv(-)     (flux vector splitting for advection)
    // -- f_diff_i ni  = f_diff                         (evaluation of diffusion flux based on domain solution)

    // prepare the side finite element
    libMesh::FEBase *fe_ptr    = NULL;
    libMesh::QBase  *qrule_ptr = NULL;
    _get_side_fe_and_qrule(_elem, s, &fe_ptr, &qrule_ptr, false);
    std::auto_ptr<libMesh::FEBase> fe(fe_ptr);
    std::auto_ptr<libMesh::QBase>  qrule(qrule_ptr);
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n1   = RealVectorX::Zero(n1),
    vec3_n2   = RealVectorX::Zero(n2),
    flux      = RealVectorX::Zero(n1),
    eig_val   = RealVectorX::Zero(n1),
    dnormal   = RealVectorX::Zero(dim);
    
    RealMatrixX
    mat1_n1n1        = RealMatrixX::Zero( n1, n1),
    mat2_n1n2        = RealMatrixX::Zero( n1, n2),
    mat3_n2n2        = RealMatrixX::Zero( n2, n2),
    leig_vec         = RealMatrixX::Zero( n1, n1),
    leig_vec_inv_tr  = RealMatrixX::Zero( n1, n1);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        
        // first update the variables at the current quadrature point
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1, _sol);
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());

        this->calculate_advection_left_eigenvector_and_inverse_for_normal
        (primitive_sol,
         normals[qp],
         eig_val,
         leig_vec,
         leig_vec_inv_tr);
        
        // for all eigenalues that are less than 0, the characteristics are coming into the domain, hence,
        // evaluate them using the given solution.
        mat1_n1n1 = leig_vec_inv_tr;
        for (unsigned int j=0; j<n1; j++)
            if (eig_val(j) < 0)
                mat1_n1n1.col(j) *= eig_val(j); // L^-T [omaga]_{-}
            else
                mat1_n1n1.col(j) *= 0.0;
        
        mat1_n1n1 *= leig_vec.transpose(); // A_{-} = L^-T [omaga]_{-} L^T
        
        this->get_infinity_vars( vec2_n1 );
        flux = mat1_n1n1 * vec2_n1;  // f_{-} = A_{-} B U
        
        Bmat.vector_mult_transpose(vec3_n2, flux); // B^T f_{-}   (this is flux coming into the solution domain)
        f += JxW[qp] * vec3_n2;
        
        // now calculate the flux for eigenvalues greater than 0,
        // the characteristics go out of the domain, so that
        // the flux is evaluated using the local solution
        mat1_n1n1 = leig_vec_inv_tr;
        for (unsigned int j=0; j<n1; j++)
            if (eig_val(j) > 0)
                mat1_n1n1.col(j) *= eig_val(j); // L^-T [omaga]_{+}
            else
                mat1_n1n1.col(j) *= 0.0;
        
        mat1_n1n1 *= leig_vec.transpose(); // A_{+} = L^-T [omaga]_{+} L^T
        flux       = mat1_n1n1 * vec1_n1; // f_{+} = A_{+} B U
        
        Bmat.vector_mult_transpose(vec3_n2, flux); // B^T f_{+}   (this is flux going out of the solution domain)
        f += JxW[qp] * vec3_n2;
        
        
        if ( request_jacobian)
        {
            // terms with negative eigenvalues do not contribute to the Jacobian
            
            // now calculate the Jacobian for eigenvalues greater than 0,
            // the characteristics go out of the domain, so that
            // the flux is evaluated using the local solution
            mat1_n1n1 = leig_vec_inv_tr;
            for (unsigned int j=0; j<n1; j++)
                if (eig_val(j) > 0)
                    mat1_n1n1.col(j) *= eig_val(j); // L^-T [omaga]_{+}
                else
                    mat1_n1n1.col(j) *= 0.0;
            mat1_n1n1 *= leig_vec.transpose(); // A_{+} = L^-T [omaga]_{+} L^T
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2); // B^T A_{+} B   (this is flux going out of the solution domain)
            
            jac += JxW[qp] * mat3_n2n2;
        }
    }
    
    return false;
}




bool
MAST::ConservativeFluidElementBase::
far_field_surface_residual_sensitivity(bool request_jacobian,
                                       RealVectorX& f,
                                       RealMatrixX& jac,
                                       const unsigned int s,
                                       MAST::BoundaryConditionBase& p) {
    
    return false;
}





void
MAST::ConservativeFluidElementBase::
_initialize_fem_interpolation_operator(const unsigned int qp,
                                       const unsigned int dim,
                                       const libMesh::FEBase& fe,
                                       MAST::FEMOperatorMatrix& Bmat) {
    
    const std::vector<std::vector<Real> >& phi_fe = fe.get_phi();
    
    const unsigned int n_phi = (unsigned int)phi_fe.size();
    
    RealVectorX phi = RealVectorX::Zero(n_phi);
    
    // shape function values
    // N
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi(i_nd) = phi_fe[i_nd][qp];
    
    Bmat.reinit(dim+2, phi);
}




void
MAST::ConservativeFluidElementBase::
_initialize_fem_gradient_operator(const unsigned int qp,
                                  const unsigned int dim,
                                  const libMesh::FEBase& fe,
                                  std::vector<MAST::FEMOperatorMatrix>& dBmat) {
    
    libmesh_assert(dBmat.size() == dim);
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi = fe.get_dphi();
    
    const unsigned int n_phi = (unsigned int)dphi.size();
    RealVectorX phi = RealVectorX::Zero(n_phi);
    
    // now set the shape function values
    for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
        
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi(i_nd) = dphi[i_nd][qp](i_dim);
        dBmat[i_dim].reinit(dim+2, phi); //  dU/dx_i
    }
}


