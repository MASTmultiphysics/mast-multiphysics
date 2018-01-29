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
#include "fluid/conservative_fluid_element_base.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/flight_condition.h"
#include "base/boundary_condition_base.h"
#include "numerics/fem_operator_matrix.h"
#include "base/system_initialization.h"
#include "elasticity/normal_rotation_function_base.h"
#include "fluid/surface_integrated_pressure_output.h"
#include "base/nonlinear_system.h"
#include "mesh/fe_base.h"
#include "base/assembly_base.h"


MAST::ConservativeFluidElementBase::
ConservativeFluidElementBase(MAST::SystemInitialization&    sys,
                             MAST::AssemblyBase&            assembly,
                             const libMesh::Elem&           elem,
                             const MAST::FlightCondition&   f):
MAST::FluidElemBase(elem.dim(), f),
MAST::ElementBase(sys, assembly, elem) {
    
    // initialize the finite element data structures
    _fe = assembly.build_fe(_elem).release();
    _fe->init(elem);
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
    LS_sens         = RealMatrixX::Zero(   n2,    n2),
    stress          = RealMatrixX::Zero(  dim,   dim),
    dprim_dcons     = RealMatrixX::Zero(   n1,    n1),
    dcons_dprim     = RealMatrixX::Zero(   n1,    n1);
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n1   = RealVectorX::Zero(n1),
    vec3_n2   = RealVectorX::Zero(n2),
    dc        = RealVectorX::Zero(dim),
    temp_grad = RealVectorX::Zero(dim);

    
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
        
        if (if_viscous()) {
            
            calculate_conservative_variable_jacobian(primitive_sol,
                                                     dcons_dprim,
                                                     dprim_dcons);
            calculate_diffusion_tensors(_sol,
                                        dBmat,
                                        dprim_dcons,
                                        primitive_sol,
                                        stress,
                                        temp_grad);
        }
        
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
            
            // diffusive flux
            if (if_viscous()) {
                
                calculate_diffusion_flux(i_dim,
                                         primitive_sol,
                                         stress,
                                         temp_grad,
                                         vec1_n1);
                dBmat[i_dim].vector_mult_transpose(vec3_n2, vec1_n1);
                f += JxW[qp] * vec3_n2;
            }

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
                
                // viscous flux Jacobian
                if (if_viscous()) {
                    
                    for (unsigned int j_dim=0; j_dim<dim; j_dim++) {
                        
                        calculate_diffusion_flux_jacobian(i_dim,
                                                          j_dim,
                                                          primitive_sol,
                                                          mat1_n1n1);
                        
                        dBmat[j_dim].left_multiply(mat3_n1n2, mat1_n1n1);                     // Kij dB_j
                        dBmat[i_dim].right_multiply_transpose(mat4_n2n2, mat3_n1n2);          // dB_i^T Kij dB_j
                        jac += JxW[qp]*mat4_n2n2;
                    }
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
MAST::ConservativeFluidElementBase::
linearized_internal_residual (bool request_jacobian,
                              RealVectorX& f,
                              RealMatrixX& jac) {
    
    // first get the internal residual and Jacobian from the
    // conservative elems, and then add to it an extra dissipation
    // to control solution discontinuities emanating from the boundary.
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealMatrixX
    f_jac_x         = RealMatrixX::Zero(   n2,    n2),
    fm_jac_xdot     = RealMatrixX::Zero(   n2,    n2);
    
    RealVectorX
    local_f    = RealVectorX::Zero(n2);
    
    // df/dx. We always need the Jacobian, since it is used to calculate
    // the residual
    this->internal_residual(true, local_f, f_jac_x);
    
    if (request_jacobian)
        jac      +=  f_jac_x;
    
    f        +=  f_jac_x * _delta_sol;
    
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
linearized_velocity_residual (bool request_jacobian,
                              RealVectorX& f,
                              RealMatrixX& jac_xdot,
                              RealMatrixX& jac) {
    
    // first get the internal residual and Jacobian from the
    // conservative elems, and then add to it an extra dissipation
    // to control solution discontinuities emanating from the boundary.
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealMatrixX
    fm_jac_x         = RealMatrixX::Zero(   n2,    n2),
    fm_jac_xdot     = RealMatrixX::Zero(   n2,    n2);
    
    RealVectorX
    local_f    = RealVectorX::Zero(n2);
    
    // df/dx. We always need the Jacobian, since it is used to calculate
    // the residual
    this->velocity_residual(true, local_f, fm_jac_xdot, fm_jac_x);
    
    if (request_jacobian) {
        
        jac_xdot +=  fm_jac_xdot;
        jac      +=  fm_jac_x;
    }
    
    
    f        +=  (fm_jac_x    * _delta_sol +
                  fm_jac_xdot * _delta_vel);
    
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
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids;
        binfo.boundary_ids(&_elem, n, bc_ids);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SYMMETRY_WALL:
                        symmetry_surface_residual(request_jacobian,
                                                  f, jac,
                                                  n,
                                                  *it.first->second);
                        break;
                        
                    case MAST::SLIP_WALL:
                        slip_wall_surface_residual(request_jacobian,
                                                   f, jac,
                                                   n,
                                                   *it.first->second);
                        break;

                    case MAST::NO_SLIP_WALL:
                        noslip_wall_surface_residual(request_jacobian,
                                                     f, jac,
                                                     n,
                                                     *it.first->second);
                        break;

                    case MAST::FAR_FIELD:
                        far_field_surface_residual(request_jacobian,
                                                   f, jac,
                                                   n,
                                                   *it.first->second);
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
    return request_jacobian;
}




bool
MAST::ConservativeFluidElementBase::
linearized_side_external_residual (bool request_jacobian,
                                   RealVectorX& f,
                                   RealMatrixX& jac,
                                   std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*>& bc) {
    
    // the far-field and symmetry wall, which are assumed to be stationary will
    // only contribute to the real part of the complex Jacobian. Hence,
    // we will use the parent class's methods for these two. The
    // slip wall, which may be oscillating, is implemented for this element.
    
    const unsigned int
    dim     = _elem.dim(),
    n1      = dim+2,
    n2      = _fe->n_shape_functions()*n1;
    
    RealMatrixX
    f_jac_x = RealMatrixX::Zero(n2, n2);
    
    RealVectorX
    local_f = RealVectorX::Zero(n2);
    
    
    typedef std::multimap<libMesh::boundary_id_type, MAST::BoundaryConditionBase*> maptype;
    
    // iterate over the boundary ids given in the provided force map
    std::pair<maptype::const_iterator, maptype::const_iterator> it;
    
    const libMesh::BoundaryInfo& binfo = *_system.system().get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids;
        binfo.boundary_ids(&_elem, n, bc_ids);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SYMMETRY_WALL: {
                        
                        f_jac_x.setZero();
                        local_f.setZero();
                        
                        // We always need the Jacobian, since it is used to
                        // calculate the residual
                        
                        this->symmetry_surface_residual(true,
                                                        local_f,
                                                        f_jac_x,
                                                        n,
                                                        *it.first->second);
                        
                        // multiply jacobian with the nondimensionalizing
                        // factor (V/b for flutter analysis)
                        if (request_jacobian)
                            jac  += f_jac_x;
                        
                        f    += f_jac_x * _delta_sol;
                    }
                        break;
                        
                    case MAST::SLIP_WALL: {
                        
                        // this calculates the Jacobian and residual contribution
                        // including the nondimensionalizing factor.
                        
                        this->linearized_slip_wall_surface_residual(request_jacobian,
                                                                    f,
                                                                    jac,
                                                                    n,
                                                                    *it.first->second);
                    }
                        break;
                        
                    case MAST::FAR_FIELD: {
                        
                        f_jac_x.setZero();
                        local_f.setZero();
                        
                        // We always need the Jacobian, since it is used to
                        // calculate the residual
                        
                        this->far_field_surface_residual(true,
                                                         local_f,
                                                         f_jac_x,
                                                         n,
                                                         *it.first->second);
                        
                        // multiply jacobian with the nondimensionalizing
                        // factor (V/b for flutter analysis)
                        if (request_jacobian)
                            jac  += f_jac_x;
                        
                        f    += f_jac_x * _delta_sol;
                    }
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
    
    
    return request_jacobian;
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
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        
        // if no boundary ids have been specified for the side, then
        // move to the next side.
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        // check to see if any of the specified boundary ids has a boundary
        // condition associated with them
        std::vector<libMesh::boundary_id_type> bc_ids;
        binfo.boundary_ids(&_elem, n, bc_ids);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            
            // find the loads on this boundary and evaluate the f and jac
            it = bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SYMMETRY_WALL:
                        symmetry_surface_residual_sensitivity(request_jacobian,
                                                              f, jac,
                                                              n,
                                                              *it.first->second);
                        break;
                        
                    case MAST::SLIP_WALL:
                        slip_wall_surface_residual_sensitivity(request_jacobian,
                                                               f, jac,
                                                               n,
                                                               *it.first->second);
                        break;
                        
                    case MAST::FAR_FIELD:
                        far_field_surface_residual_sensitivity(request_jacobian,
                                                               f, jac,
                                                               n,
                                                               *it.first->second);
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
    return request_jacobian;
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
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, false);

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
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, false);

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();

    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    flux      = RealVectorX::Zero(n1),
    dnormal   = RealVectorX::Zero(dim),
    uvec      = RealVectorX::Zero(3),
    ni        = RealVectorX::Zero(3),
    vel_fe    = RealVectorX::Zero(6),
    dwdot_i   = RealVectorX::Zero(3),
    dni       = RealVectorX::Zero(3);
    
    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    Real
    ui_ni = 0.;

    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    
    
    // get the surface motion object from the boundary condition object
    MAST::FieldFunction<RealVectorX>
    *vel   = nullptr;
    MAST::NormalRotationFunctionBase<RealVectorX>
    *n_rot = nullptr;
    
    if (p.contains("velocity"))
        vel = &p.get<MAST::FieldFunction<RealVectorX> >("velocity");
    if (p.contains("normal_rotation")) {
        
        MAST::FieldFunction<RealVectorX>&
        tmp = p.get<MAST::FieldFunction<RealVectorX> >("normal_rotation");
        n_rot = dynamic_cast<MAST::NormalRotationFunctionBase<RealVectorX>*>(&tmp);
    }

    // if displ is provided then n_rot must also be provided
    if (vel) libmesh_assert(n_rot);
    
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
        
        ////////////////////////////////////////////////////////////
        //   Calculation of the surface velocity term.
        //        vi (ni + dni)  =  wdot_i (ni + dni)
        //  or    vi ni = wdot_i (ni + dni) - vi dni
        ////////////////////////////////////////////////////////////
        
        
        // copy the surface normal
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            ni(i_dim) = normals[qp](i_dim);
        
        
        // now check if the surface deformation is defined and
        // needs to be applied through transpiration boundary
        // condition

        primitive_sol.get_uvec(uvec);

        if (vel) { // get the surface motion data
            
            (*vel)(qpoint[qp], _time, vel_fe);
            dwdot_i = vel_fe.topRows(3);
        }
        
        if (n_rot)
            (*n_rot)(qpoint[qp], normals[qp], _time, dni);
            
        ui_ni  = dwdot_i.dot(ni+dni) - uvec.dot(dni);
        
        
        flux.setZero();
        flux                 += ui_ni * vec1_n1;           // vi ni cons_flux
        flux(n1-1)           += ui_ni * primitive_sol.p;   // vi ni {0, 0, 0, 0, p}
        flux.segment(1,dim)  += primitive_sol.p * ni.segment(0,dim);      // p * ni for the momentum eqs.
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;
        
        if ( request_jacobian ) {
            
            this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             ui_ni,
             normals[qp],
             dni,
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
linearized_slip_wall_surface_residual(bool request_jacobian,
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
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, false);

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    
    const unsigned int
    dim        = _elem.dim(),
    n1         = dim+2,
    n2         = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1    = RealVectorX::Zero(n1),
    uvec       = RealVectorX::Zero(3),
    dwdot_i    = RealVectorX::Zero(3),
    ni         = RealVectorX::Zero(3),
    dni        = RealVectorX::Zero(3),
    Dw_i       = RealVectorX::Zero(3),
    Dni        = RealVectorX::Zero(3),
    Duvec      = RealVectorX::Zero(3),
    vec2_n1    = RealVectorX::Zero(n1),
    vec2_n2    = RealVectorX::Zero(n2),
    tmp        = RealVectorX::Zero(6),
    flux       = RealVectorX::Zero(n1);
    
    RealMatrixX
    mat1_n1n1  = RealMatrixX::Zero( n1, n1),
    mat2_n1n2  = RealMatrixX::Zero( n1, n2),
    mat3_n2n2  = RealMatrixX::Zero( n2, n2);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution                            primitive_sol;
    MAST::SmallPerturbationPrimitiveSolution<Real>  sd_primitive_sol;
    
    Real
    Dvi_ni              = 0.,
    ui_ni_steady        = 0.;
    
    
    // get the surface motion object from the boundary condition object
    MAST::FieldFunction<RealVectorX>
    *vel   = nullptr;
    MAST::NormalRotationFunctionBase<RealVectorX>
    *n_rot = nullptr;
    
    if (p.contains("velocity"))
        vel = &p.get<MAST::FieldFunction<RealVectorX> >("velocity");
    if (p.contains("normal_rotation")) {
        
        MAST::FieldFunction<RealVectorX>&
        tmp = p.get<MAST::FieldFunction<RealVectorX> >("normal_rotation");
        n_rot = dynamic_cast<MAST::NormalRotationFunctionBase<RealVectorX>*>(&tmp);
    }
    
    // if displ is provided then n_rot must also be provided
    if (vel) libmesh_assert(n_rot);

    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        Bmat.right_multiply(vec1_n1,         _sol);  // conservative sol
        Bmat.right_multiply(vec2_n1,   _delta_sol);  // perturbation sol
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());
        
        // initialize the small-disturbance primitive sol
        sd_primitive_sol.zero();
        sd_primitive_sol.init(primitive_sol, vec2_n1);
        
        ////////////////////////////////////////////////////////////
        //   Calculation of the surface velocity term. For a
        //   steady-state system,
        //        vi (ni + dni)  =  wdot_i (ni + dni)
        //  or    vi ni = wdot_i (ni + dni) - vi dni
        //
        //    The first order perturbation, D(.), of this is
        //        (vi + Dvi) (ni + dni + Dni) =
        //          (wdot_i + Dwdot_i) (ni + dni + Dni)
        //  or    Dvi ni = Dwdot_i (ni + dni + Dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi (dni + Dni)
        //     Neglecting HOT, this becomes
        //        Dvi ni = Dwdot_i (ni + dni) +
        //                  wdot_i (ni + dni + Dni) -
        //                  vi (ni + dni + Dni) -
        //                  Dvi dni
        //               = Dwdot_i (ni + dni) +
        //                 wdot_i Dni - vi Dni - Dvi dni
        //
        ////////////////////////////////////////////////////////////
        
        // copy the surface normal
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            ni(i_dim) = normals[qp](i_dim);
        
        // now check if the surface deformation is defined and
        // needs to be applied through transpiration boundary
        // condition
        primitive_sol.get_uvec(uvec);
        sd_primitive_sol.get_duvec(Duvec);
        
        flux.setZero();
        
        if (vel) {
            
            // base flow contribution
            (*vel)(qpoint[qp], _time, tmp);
            dwdot_i = tmp.topRows(3);

            // small-disturbance contribution
            vel->perturbation(qpoint[qp], _time, tmp);
            Dw_i = tmp.topRows(3);
        }
        
        if (n_rot) {
            
            // base flow contribution
            (*n_rot)(qpoint[qp], normals[qp], _time, dni);
        
            // small-disturbance contribution
            n_rot->perturbation(qpoint[qp], normals[qp], _time, Dni);
        }
        
        // base flow contribution to the flux
        ui_ni_steady  =  dwdot_i.dot(ni+dni) - uvec.dot(dni);
        flux         += ui_ni_steady * vec2_n1;               // vi_ni  dcons_flux
        flux(n1-1)   += ui_ni_steady * sd_primitive_sol.dp;   // vi_ni {0,0,0,0,Dp}
        
        // perturbed quantity contribution to the flux
        Dvi_ni            = (Dw_i.dot(ni+dni) +
                             dwdot_i.dot(Dni) -
                             uvec.dot(Dni) -
                             Duvec.dot(dni));
        flux             += Dvi_ni * vec1_n1;              // Dvi_ni cons_flux
        flux(n1-1)       += Dvi_ni * primitive_sol.p;      // Dvi_ni {0,0,0,0,p}

        // ni Dp  (only for the momentun eq)
        flux.segment(1, dim) += (sd_primitive_sol.dp * ni.segment(0,dim));
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;
        
        if ( request_jacobian ) {
            
            // the Jacobian contribution comes only from the dp term in the
            // residual
            this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             ui_ni_steady,
             normals[qp],
             dni,
             mat1_n1n1);
            
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2);
            jac += JxW[qp] * mat3_n2n2;
        }
    }
    
    return request_jacobian;
}





bool
MAST::ConservativeFluidElementBase::
noslip_wall_surface_residual(bool request_jacobian,
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
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, true);

    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    const std::vector<libMesh::Point>& qpoint    = fe->get_xyz();
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    flux      = RealVectorX::Zero(n1),
    dnormal   = RealVectorX::Zero(dim),
    uvec      = RealVectorX::Zero(3),
    ni        = RealVectorX::Zero(3),
    dwdot_i   = RealVectorX::Zero(3),
    dni       = RealVectorX::Zero(3),
    temp_grad = RealVectorX::Zero(dim);
    
    RealMatrixX
    mat1_n1n1       = RealMatrixX::Zero( n1, n1),
    mat2_n1n2       = RealMatrixX::Zero( n1, n2),
    mat3_n2n2       = RealMatrixX::Zero( n2, n2),
    stress          = RealMatrixX::Zero(  dim,   dim),
    dprim_dcons     = RealMatrixX::Zero(   n1,    n1),
    dcons_dprim     = RealMatrixX::Zero(   n1,    n1);
    
    Real
    ui_ni = 0.;
    
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    std::vector<MAST::FEMOperatorMatrix> dBmat(dim);
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    
    
    // get the surface motion object from the boundary condition object
    MAST::FieldFunction<RealVectorX>
    *vel   = nullptr;
    MAST::NormalRotationFunctionBase<RealVectorX>
    *n_rot = nullptr;
    
    if (p.contains("velocity"))
        vel = &p.get<MAST::FieldFunction<RealVectorX> >("velocity");
    if (p.contains("normal_rotation")) {
        
        MAST::FieldFunction<RealVectorX>&
        tmp = p.get<MAST::FieldFunction<RealVectorX> >("normal_rotation");
        n_rot = dynamic_cast<MAST::NormalRotationFunctionBase<RealVectorX>*>(&tmp);
    }

    
    for (unsigned int qp=0; qp<JxW.size(); qp++)
    {
        // initialize the Bmat operator for this term
        _initialize_fem_interpolation_operator(qp, dim, *fe, Bmat);
        _initialize_fem_gradient_operator(qp, dim, *fe, dBmat);
        
        Bmat.right_multiply(vec1_n1, _sol);
        
        // initialize the primitive solution
        primitive_sol.zero();
        primitive_sol.init(dim,
                           vec1_n1,
                           flight_condition->gas_property.cp,
                           flight_condition->gas_property.cv,
                           if_viscous());

        // copy the surface normal
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            ni(i_dim) = normals[qp](i_dim);
        

        /*////////////////////////////////////////////////////////////
        //   Calculation of the surface velocity term.
        //        vi (ni + dni)  =  wdot_i (ni + dni)
        //  or    vi ni = wdot_i (ni + dni) - vi dni
        ////////////////////////////////////////////////////////////
        
        // now check if the surface deformation is defined and
        // needs to be applied through transpiration boundary
        // condition
        
        primitive_sol.get_uvec(uvec);
        
        if (motion) { // get the surface motion data
            
            motion->time_domain_motion(_time,
                                       qpoint[qp],
                                       normals[qp],
                                       dwdot_i,
                                       dni);
            
            ui_ni  = dwdot_i.dot(ni+dni) - uvec.dot(dni);
        }*/
        
        
        flux.setZero();
        
        // convective flux terms
        flux                 += ui_ni * vec1_n1;           // vi ni cons_flux
        flux(n1-1)           += ui_ni * primitive_sol.p;   // vi ni {0, 0, 0, 0, p}
        flux.segment(1,dim)  += primitive_sol.p * ni.segment(0,dim);      // p * ni for the momentum eqs.
        
        // now, add the viscous flux terms
        calculate_conservative_variable_jacobian(primitive_sol,
                                                 dcons_dprim,
                                                 dprim_dcons);
        calculate_diffusion_tensors(_sol,
                                    dBmat,
                                    dprim_dcons,
                                    primitive_sol,
                                    stress,
                                    temp_grad);

        temp_grad.setZero(); // assuming adiabatic for now
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            
            calculate_diffusion_flux(i_dim,
                                     primitive_sol,
                                     stress,
                                     temp_grad,
                                     vec1_n1);
            //flux -= ni(i_dim) * vec1_n1;
        }
        
        
        Bmat.vector_mult_transpose(vec2_n2, flux);
        f += JxW[qp] * vec2_n2;

        if ( true) { //request_jacobian ) {
            
            this->calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
            (primitive_sol,
             ui_ni,
             normals[qp],
             dni,
             mat1_n1n1);
            
            Bmat.left_multiply(mat2_n1n2, mat1_n1n1);
            Bmat.right_multiply_transpose(mat3_n2n2, mat2_n1n2);
            jac += JxW[qp] * mat3_n2n2;
            
            for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
                
                for (unsigned int j_dim=0; j_dim<dim; j_dim++) {
                    
                    calculate_diffusion_flux_jacobian(i_dim,
                                                      j_dim,
                                                      primitive_sol,
                                                      mat1_n1n1);
                    
                    dBmat[j_dim].left_multiply(mat2_n1n2, mat1_n1n1);                     // Kij dB_j
                    dBmat[i_dim].right_multiply_transpose(mat3_n2n2, mat2_n1n2);          // dB_i^T Kij dB_j
                                                                                          //jac -= JxW[qp]*mat3_n2n2;
                }
            }
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
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, false);

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
_calculate_surface_integrated_load(bool request_derivative,
                                   bool request_sensitivity,
                                   const unsigned int s,
                                   MAST::OutputAssemblyElemOperations& output) {
    
    // prepare the side finite element
    std::unique_ptr<MAST::FEBase> fe(_assembly.build_fe(_elem));
    fe->init_for_side(_elem, s, false);
    
    const std::vector<Real> &JxW                 = fe->get_JxW();
    const std::vector<libMesh::Point>& normals   = fe->get_normals();
    
    const unsigned int
    dim    = _elem.dim(),
    n1     = dim+2,
    n2     = _fe->n_shape_functions()*n1;
    
    RealVectorX
    vec1_n1   = RealVectorX::Zero(n1),
    vec2_n2   = RealVectorX::Zero(n2),
    load      = RealVectorX::Zero(3),
    load_sens = RealVectorX::Zero(3);
    
    RealMatrixX
    dload_dX   = RealMatrixX::Zero( 3, n1);
    
    libMesh::Point pt;
    MAST::FEMOperatorMatrix Bmat;
    
    // create objects to calculate the primitive solution, flux, and Jacobian
    MAST::PrimitiveSolution      primitive_sol;
    MAST::SmallPerturbationPrimitiveSolution<Real>  primitive_sol_sens;
    
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
        
        
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            load(i_dim) += JxW[qp] * primitive_sol.p * normals[qp](i_dim);
        

        // calculate the derivative, if requested
        if (request_derivative) {
            
            calculate_pressure_derivative_wrt_conservative_variables
            (primitive_sol, vec1_n1);
            Bmat.vector_mult_transpose(vec2_n2, vec1_n1);
            
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                dload_dX.row(i_dim) += JxW[qp] * vec2_n2 * normals[qp](i_dim);
        }
        
        

        // calculate the sensitivity, if requested
        if (request_sensitivity) {
            
            // calculate the solution sensitivity
            Bmat.right_multiply(vec1_n1, _sol_sens);
            
            // initialize the perturbation in primite solution, which will
            // give us sensitivity of pressure
            primitive_sol_sens.zero();
            primitive_sol_sens.init(primitive_sol, vec1_n1);
            
            for (unsigned int i_dim=0; i_dim<dim; i_dim++)
                load_sens(i_dim) += JxW[qp] *
                (primitive_sol_sens.dp * normals[qp](i_dim) +
                 primitive_sol.p       * 0.);
        }
    }
}




void
MAST::ConservativeFluidElementBase::
_initialize_fem_interpolation_operator(const unsigned int qp,
                                       const unsigned int dim,
                                       const MAST::FEBase& fe,
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
                                  const MAST::FEBase& fe,
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


