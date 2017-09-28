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

// C++ includes
#include <iomanip>

// MAST includes
#include "fluid/fluid_elem_base.h"
#include "fluid/primitive_fluid_solution.h"
#include "fluid/small_disturbance_primitive_fluid_solution.h"
#include "fluid/flight_condition.h"

// Basic include files
#include "libmesh/mesh.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"



MAST::FluidElemBase::FluidElemBase(const unsigned int d,
                                   const MAST::FlightCondition& f):
_if_viscous(f.gas_property.if_viscous),
_include_pressure_switch(false),
flight_condition(&f),
dim(d),
_dissipation_scaling(1.) {
    
    
    // prepare the variable vector
    _active_primitive_vars.push_back(RHO_PRIM);
    _active_primitive_vars.push_back(VEL1);
    _active_conservative_vars.push_back(RHO_CONS);
    _active_conservative_vars.push_back(RHOVEL1);
    
    if (dim > 1)
    {
        _active_primitive_vars.push_back(VEL2);
        _active_conservative_vars.push_back(RHOVEL2);
    }
    
    if (dim > 2)
    {
        _active_primitive_vars.push_back(VEL3);
        _active_conservative_vars.push_back(RHOVEL3);
    }
    _active_primitive_vars.push_back(TEMP);
    _active_conservative_vars.push_back(ETOT);
    
}




MAST::FluidElemBase::~FluidElemBase() {
    
    
}





void MAST::FluidElemBase::
get_infinity_vars( RealVectorX& vars_inf ) const {
    
    vars_inf(0)       = flight_condition->rho();
    
    vars_inf(1)       = flight_condition->rho_u1();
    
    if (dim > 1)
        vars_inf(2)   = flight_condition->rho_u2();
    
    if (dim > 2)
        vars_inf(3)   = flight_condition->rho_u3();
    
    vars_inf(dim+2-1) = flight_condition->rho_e();
}




void
MAST::FluidElemBase::
update_solution_at_quadrature_point(const unsigned int qp,
                                    const libMesh::FEBase& fe,
                                    const RealVectorX& elem_solution,
                                    RealVectorX& conservative_sol,
                                    PrimitiveSolution& primitive_sol,
                                    MAST::FEMOperatorMatrix& B_mat,
                                    std::vector<MAST::FEMOperatorMatrix>& dB_mat) {
    
    conservative_sol.setZero();
    
    
    const std::vector<std::vector<Real> >& phi = fe.get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();

    RealVectorX
    phi_vals               = RealVectorX::Zero(n_phi);
    
    for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
        phi_vals(i_nd) = phi[i_nd][qp];
    
    B_mat.reinit(dim+2, phi_vals); // initialize the operator matrix
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi =
    fe.get_dphi();
    
    
    for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
    {
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vals(i_nd) = dphi[i_nd][qp](i_dim);
        dB_mat[i_dim].reinit(dim+2, phi_vals);
    }
    
    B_mat.vector_mult( conservative_sol, elem_solution );
    
    primitive_sol.zero();
    primitive_sol.init(dim,
                       conservative_sol,
                       flight_condition->gas_property.cp,
                       flight_condition->gas_property.cv,
                       _if_viscous);
}




void
MAST::FluidElemBase::calculate_advection_flux(const unsigned int calculate_dim,
                                              const MAST::PrimitiveSolution& sol,
                                              RealVectorX& flux) {
    
    const unsigned int n1 = 2 + dim;
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    p = sol.p,
    e_tot = sol.e_tot;
    
    flux.setZero();
    
    // calculate the flux using given flow parameters at this point
    switch (calculate_dim)
    {
        case 0:
        {
            flux(0) =  rho * u1;
            flux(n1-1) =  u1 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u1 * u3;
                case 2:
                    flux(2) =  rho * u1 * u2;
                case 1:
                    flux(1) =  rho * u1 * u1 + p;
            }
        }
            break;
            
        case 1:
        {
            flux(0) =  rho * u2;
            flux(n1-1) =  u2 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u2 * u3;
                case 2:
                    flux(2) =  rho * u2 * u2 + p;
                case 1:
                    flux(1) =  rho * u2 * u1;
            }
        }
            break;
            
        case 2:
        {
            flux(0) =  rho * u3;
            flux(n1-1) =  u3 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u3 * u3 + p;
                case 2:
                    flux(2) =  rho * u3 * u2;
                case 1:
                    flux(1) =  rho * u3 * u1;
            }
        }
            break;
    }
}




void
MAST::FluidElemBase::calculate_diffusion_flux(const unsigned int calculate_dim,
                                              const MAST::PrimitiveSolution& sol,
                                              const RealMatrixX& stress_tensor,
                                              const RealVectorX& temp_gradient,
                                              RealVectorX& flux) {
    
    const unsigned int n1 = 2 + dim;
    
    RealVectorX
    uvec               = RealVectorX::Zero(3);

    uvec(0) = sol.u1;
    uvec(1) = sol.u2;
    uvec(2) = sol.u3;
    
    flux.setZero();
    
    // the momentum flux is obtained from the stress tensor
    for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
        
        flux(1+i_dim) += stress_tensor(calculate_dim, i_dim); // tau_ij
        
        flux(n1-1) += uvec(i_dim) * stress_tensor(calculate_dim, i_dim); // u_j tau_ij
    }

    flux(n1-1) += sol.k_thermal * temp_gradient(calculate_dim);
}




void
MAST::FluidElemBase::
calculate_diffusion_tensors(const RealVectorX& elem_sol,
                            const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                            const RealMatrixX& dprim_dcons,
                            const MAST::PrimitiveSolution& psol,
                            RealMatrixX& stress_tensor,
                            RealVectorX& temp_gradient) {
    
    const unsigned int n1 = dim+2;
    
    RealVectorX
    dprim_dx               = RealVectorX::Zero(n1),
    dcons_dx               = RealVectorX::Zero(n1);
    
    stress_tensor.setZero();
    temp_gradient.setZero();
    
    for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
        
        dB_mat[i_dim].vector_mult(dcons_dx, elem_sol); // dUcons/dx_i
        dprim_dx = dprim_dcons * dcons_dx; // dUprim/dx_i
        
        for (unsigned int j_dim=0; j_dim<dim; j_dim++) {
            
            stress_tensor(i_dim, j_dim) += psol.mu * dprim_dx(j_dim+1); // mu * duj/dxi
            stress_tensor(j_dim, i_dim) += psol.mu * dprim_dx(j_dim+1); // mu * dui/dxj
            
            stress_tensor(j_dim, j_dim) += psol.lambda * dprim_dx(i_dim+1);  // + delta_ij lambda duk/dxk
        }
        
        temp_gradient(i_dim) = dprim_dx(n1-1); // dT/dx_i
    }
}




void
MAST::FluidElemBase::
calculate_conservative_variable_jacobian(const MAST::PrimitiveSolution& sol,
                                         RealMatrixX& dcons_dprim,
                                         RealMatrixX& dprim_dcons) {
    
    
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    dcons_dprim.setZero();
    dprim_dcons.setZero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    rho = sol.rho,
    k = sol.k,
    e_tot = sol.e_tot,
    cv = flight_condition->gas_property.cv;
    
    switch (dim)
    {
        case 3:
        {
            dcons_dprim(3, 0) = u3;
            dcons_dprim(3, 3) = rho;
            
            dcons_dprim(n1-1, 3) = rho*u3;
        }
            
        case 2:
        {
            dcons_dprim(2, 0) = u2;
            dcons_dprim(2, 2) = rho;
            
            dcons_dprim(n1-1, 2) = rho*u2;
        }
            
        case 1:
        {
            dcons_dprim(0, 0) = 1.;
            
            dcons_dprim(1, 0) = u1;
            dcons_dprim(1, 1) = rho;
            
            dcons_dprim(n1-1, 0) = e_tot;
            dcons_dprim(n1-1, 1) = rho*u1;
            dcons_dprim(n1-1, n1-1) = rho*cv;
        }
    }
    
    switch (dim)
    {
        case 3:
        {
            dprim_dcons(3, 0) = -u3/rho;
            dprim_dcons(3, 3) = 1./rho;
            
            dprim_dcons(n1-1, 3) = -u3/cv/rho;
        }
            
        case 2:
        {
            dprim_dcons(2, 0) = -u2/rho;
            dprim_dcons(2, 2) = 1./rho;
            
            dprim_dcons(n1-1, 2) = -u2/cv/rho;
        }
            
        case 1:
        {
            dprim_dcons(0, 0) = 1.;
            
            dprim_dcons(1, 0) = -u1/rho;
            dprim_dcons(1, 1) = 1./rho;
            
            dprim_dcons(n1-1, 0) = (-e_tot+2*k)/cv/rho;
            dprim_dcons(n1-1, 1) = -u1/cv/rho;
            dprim_dcons(n1-1, n1-1) = 1./cv/rho;
        }
    }
    
}



void
MAST::FluidElemBase::
calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                  const MAST::PrimitiveSolution& sol,
                                  RealMatrixX& mat) {
    
    
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.setZero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    e_tot = sol.e_tot,
    T = sol.T,
    gamma = flight_condition->gas_property.gamma,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(1, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u1*u3;
                    mat(3, 1) = u3;
                    mat(3, 3) = u1;
                    
                    mat(n1-1, 3) = -u1*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(1, 2) = -u2*R/cv;
                    
                    mat(2, 0) = -u1*u2;
                    mat(2, 1) = u2;
                    mat(2, 2) = u1;
                    
                    mat(n1-1, 2) = -u1*u2*R/cv;
                }
                    
                case 1:
                {
                    mat(0, 1) = 1.0; // d U / d (rho u1)
                    
                    mat(1, 0) = -u1*u1+R*k/cv;
                    mat(1, 1) = u1*(2.0-R/cv);
                    mat(1, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u1*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = e_tot+R*(T-u1*u1/cv);
                    mat(n1-1, n1-1) = u1*gamma;
                }
                    break;
            }
        }
            break;
            
        case 1:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(2, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u2*u3;
                    mat(3, 2) = u3;
                    mat(3, 3) = u2;
                    
                    mat(n1-1, 3) = -u2*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(0, 2) = 1.0; // d U / d (rho u2)
                    
                    mat(1, 0) = -u1*u2;
                    mat(1, 1) = u2;
                    mat(1, 2) = u1;
                    
                    mat(2, 0) = -u2*u2+R*k/cv;
                    mat(2, 1) = -u1*R/cv;
                    mat(2, 2) = u2*(2.0-R/cv);
                    mat(2, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u2*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = -u1*u2*R/cv;
                    mat(n1-1, 2) = e_tot+R*(T-u2*u2/cv);
                    mat(n1-1, n1-1) = u2*gamma;
                }
                    break;
                    
                case 1:
                    // if second coordinate divergence is being asked for, then the element is atleast 2D
                    libmesh_assert_msg(false, "invalid dim");
                    break;
            }
        }
            break;
            
        case 2:
        {
            mat(0, 3) = 1.0; // d U / d (rho u3)
            
            mat(1, 0) = -u1*u3;
            mat(1, 1) = u3;
            mat(1, 3) = u1;
            
            mat(2, 0) = -u2*u3;
            mat(2, 2) = u3;
            mat(2, 3) = u2;
            
            mat(3, 0) = -u3*u3+R*k/cv;
            mat(3, 1) = -u1*R/cv;
            mat(3, 2) = -u2*R/cv;
            mat(3, 3) = u3*(2.0-R/cv);
            mat(3, n1-1) = R/cv;
            
            mat(n1-1, 0) = u3*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
            mat(n1-1, 1) = -u1*u3*R/cv;
            mat(n1-1, 2) = -u2*u3*R/cv;
            mat(n1-1, 3) = e_tot+R*(T-u3*u3/cv);
            mat(n1-1, n1-1) = u3*gamma;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}





void
MAST::FluidElemBase::
calculate_diffusion_flux_jacobian (const unsigned int flux_dim,
                                   const unsigned int deriv_dim,
                                   const MAST::PrimitiveSolution& sol,
                                   RealMatrixX& mat) {
    
    const unsigned int n1 = 2 + dim;
    
    mat.setZero();
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    e_tot = sol.e_tot,
    mu = sol.mu,
    lambda = sol.lambda,
    kth = sol.k_thermal,
    cv = flight_condition->gas_property.cv;
    
    switch (flux_dim)
    {
        case 0:
        {
            switch (deriv_dim)
            {
                case 0: // K11
                {
                    switch (dim)
                    {
                        case 3:
                        {
                            mat(3,0) = -u3*mu/rho;
                            mat(3,3) = mu/rho;
                            
                            mat(n1-1,3) = u3*(-kth+cv*mu)/cv/rho;
                        }
                            
                        case 2:
                        {
                            mat(2,0) = -u2*mu/rho;
                            mat(2,2) = mu/rho;
                            
                            mat(n1-1,2) = u2*(-kth+cv*mu)/cv/rho;
                        }
                            
                        case 1:
                        {
                            mat(1,0) = -u1*(lambda+2.*mu)/rho;
                            mat(1,1) = (lambda+2.*mu)/rho;
                            
                            mat(n1-1,0) = (kth*(2.*k-e_tot)-cv*(2.*k*mu+u1*u1*(mu+lambda)))/cv/rho;
                            mat(n1-1,1) = u1*(-kth+cv*(lambda+2.*mu))/cv/rho;
                            mat(n1-1,n1-1) = kth/cv/rho;
                        }
                    }
                }
                    break;
                    
                case 1: // K12
                {
                    mat(1,0) = -u2*lambda/rho;
                    mat(1,2) = lambda/rho;
                    
                    mat(2,0) = -u1*mu/rho;
                    mat(2,1) = mu/rho;
                    
                    mat(n1-1,0) = -u1*u2*(lambda+mu)/rho;
                    mat(n1-1,1) = u2*mu/rho;
                    mat(n1-1,2) = u1*lambda/rho;
                }
                    break;
                    
                case 2: // K13
                {
                    mat(1,0) = -u3*lambda/rho;
                    mat(1,3) = lambda/rho;
                    
                    mat(3,0) = -u1*mu/rho;
                    mat(3,1) = mu/rho;
                    
                    mat(n1-1,0) = -u1*u3*(lambda+mu)/rho;
                    mat(n1-1,1) = u3*mu/rho;
                    mat(n1-1,3) = u1*lambda/rho;
                }
                    break;
            }
        }
            break;
            
        case 1:
        {
            switch (deriv_dim)
            {
                case 0: // K21
                {
                    mat(1,0) = -u2*mu/rho;
                    mat(1,2) = mu/rho;
                    
                    mat(2,0) = -u1*lambda/rho;
                    mat(2,1) = lambda/rho;
                    
                    mat(n1-1,0) = -u1*u2*(lambda+mu)/rho;
                    mat(n1-1,1) = u2*lambda/rho;
                    mat(n1-1,2) = u1*mu/rho;
                }
                    break;
                    
                case 1: // K22
                {
                    switch (dim)
                    {
                        case 3:
                        {
                            mat(3,0) = -u3*mu/rho;
                            mat(3,3) = mu/rho;
                            
                            mat(n1-1,3) = u3*(-kth+cv*mu)/cv/rho;
                        }
                            
                        case 2:
                        case 1:
                        {
                            mat(1,0) = -u1*mu/rho;
                            mat(1,1) = mu/rho;
                            
                            mat(2,0) = -u2*(lambda+2.*mu)/rho;
                            mat(2,2) = (lambda+2.*mu)/rho;
                            
                            mat(n1-1,0) = (kth*(2.*k-e_tot)-cv*(2.*k*mu+u2*u2*(mu+lambda)))/cv/rho;
                            mat(n1-1,1) = u1*(-kth+cv*mu)/cv/rho;
                            mat(n1-1,2) = u2*(-kth+cv*(lambda+2.*mu))/cv/rho;
                            mat(n1-1,n1-1) = kth/cv/rho;
                        }
                    }
                }
                    break;
                    
                case 2: // K23
                {
                    mat(2,0) = -u3*lambda/rho;
                    mat(2,3) = lambda/rho;
                    
                    mat(3,0) = -u2*mu/rho;
                    mat(3,2) = mu/rho;
                    
                    mat(n1-1,0) = -u2*u3*(lambda+mu)/rho;
                    mat(n1-1,2) = u3*mu/rho;
                    mat(n1-1,3) = u2*lambda/rho;
                }
                    break;
            }
        }
            break;
            
        case 2:
        {
            switch (deriv_dim)
            {
                case 0: // K31
                {
                    mat(1,0) = -u3*mu/rho;
                    mat(1,3) = mu/rho;
                    
                    mat(3,0) = -u1*lambda/rho;
                    mat(3,1) = lambda/rho;
                    
                    mat(n1-1,0) = -u1*u3*(lambda+mu)/rho;
                    mat(n1-1,1) = u3*lambda/rho;
                    mat(n1-1,3) = u1*mu/rho;
                }
                    break;
                    
                case 1: // K32
                {
                    mat(2,0) = -u3*mu/rho;
                    mat(2,3) = mu/rho;
                    
                    mat(3,0) = -u2*lambda/rho;
                    mat(3,2) = lambda/rho;
                    
                    mat(n1-1,0) = -u2*u3*(lambda+mu)/rho;
                    mat(n1-1,2) = u3*lambda/rho;
                    mat(n1-1,3) = u2*mu/rho;
                }
                    break;
                    
                case 2: // K33
                {
                    mat(1,0) = -u1*mu/rho;
                    mat(1,1) = mu/rho;
                    
                    mat(2,0) = -u2*mu/rho;
                    mat(2,2) = mu/rho;
                    
                    mat(3,0) = -u3*(lambda+2.*mu)/rho;
                    mat(3,3) = (lambda+2.*mu)/rho;
                    
                    mat(n1-1,0) = (kth*(2.*k-e_tot)-cv*(2.*k*mu+u3*u3*(mu+lambda)))/cv/rho;
                    mat(n1-1,1) = u1*(-kth+cv*mu)/cv/rho;
                    mat(n1-1,2) = u2*(-kth+cv*mu)/cv/rho;
                    mat(n1-1,3) = u3*(-kth+cv*(lambda+2.*mu))/cv/rho;
                    mat(n1-1,n1-1) = kth/cv/rho;
                }
                    break;
            }
        }
            break;
    }
}




void
MAST::FluidElemBase::
calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
(const unsigned int calculate_dim,
 const MAST::PrimitiveSolution& sol,
 std::vector<RealMatrixX >& jac) {
    
    const unsigned int n1 = 2 + dim;
    
    RealMatrixX
    dprim_dcons      = RealMatrixX::Zero(n1, n1),
    mat              = RealMatrixX::Zero(n1, n1);
    
    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        jac[i_cvar].setZero();
    
    this->calculate_conservative_variable_jacobian(sol, mat, dprim_dcons);
    
    // calculate based on chain rule of the primary variables
    for (unsigned int i_pvar=0; i_pvar<n1; i_pvar++) // iterate over the primitive variables for chain rule
    {
        this->calculate_advection_flux_jacobian_sensitivity_for_primitive_variable
        (calculate_dim, i_pvar, sol, mat);
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        {
            if (fabs(dprim_dcons(i_pvar, i_cvar)) > 0.0)
                jac[i_cvar] += dprim_dcons(i_pvar, i_cvar) * mat;
        }
    }
}



void
MAST::FluidElemBase::
calculate_advection_flux_jacobian_sensitivity_for_primitive_variable
(const unsigned int calculate_dim,
 const unsigned int primitive_var,
 const MAST::PrimitiveSolution& sol,
 RealMatrixX& mat) {
    
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.setZero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    e_tot = sol.e_tot,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:  // Ax
        {
            switch (_active_primitive_vars[primitive_var])
            {
                case RHO_PRIM:
                {
                    // nothing to be done for density; matrix is zero
                }
                    break;
                    
                case VEL1:
                {
                    switch (dim)
                    {
                        case 3:
                        {
                            mat(3, 0) = -u3;
                            mat(3, 3) = 1.0;
                            
                            mat(n1-1, 3) = -u3*R/cv;
                        }
                            
                        case 2:
                        {
                            mat(2, 0) = -u2;
                            mat(2, 2) = 1.0;
                            
                            mat(n1-1, 2) = -u2*R/cv;
                        }
                            
                        case 1:
                        {
                            mat(1, 0) = -u1*(2.0-R/cv);
                            mat(1, 1) =     (2.0-R/cv);
                            
                            mat(n1-1, 0) = (-e_tot*(cv+R)+2.0*R*k + u1*u1*(-cv+R))/cv;
                            mat(n1-1, 1) = u1*(1.0-2.0*R/cv);
                            mat(n1-1, n1-1) = (cv+R)/cv;
                        }
                            break;
                    }
                }
                    break;
                    
                case VEL2:
                {
                    mat(1, 0) =  u2*R/cv;
                    mat(1, 2) =    -R/cv;
                    
                    mat(2, 0) = -u1;
                    mat(2, 1) = 1.0;
                    
                    mat(n1-1, 0) = (-cv+R)*u1*u2/cv;
                    mat(n1-1, 1) = u2;
                    mat(n1-1, 2) = -u1*R/cv;
                }
                    break;
                    
                case VEL3:
                {
                    mat(1, 0) =  u3*R/cv;
                    mat(1, 3) =    -R/cv;
                    
                    mat(3, 0) = -u1;
                    mat(3, 1) = 1.0;
                    
                    mat(n1-1, 0) = (-cv+R)*u1*u3/cv;
                    mat(n1-1, 1) = u3;
                    mat(n1-1, 3) = -u1*R/cv;
                }
                    break;
                    
                case TEMP:
                {
                    // all dimensions have the same format
                    mat(n1-1, 0) = -u1*(cv+R);
                    mat(n1-1, 1) = cv+R;
                }
                    break;
                    
                default:
                    libmesh_assert_msg(false, "Invalid primitive variable number");
                    break;
            }
        }
            break;
            
        case 1:  // Ay
        {
            switch (_active_primitive_vars[primitive_var])
            {
                case RHO_PRIM:
                {
                    // nothing to be done for density; matrix is zero
                }
                    break;
                    
                case VEL1:
                {
                    mat(1, 0) = -u2;
                    mat(1, 2) = 1.0;
                    
                    mat(2, 0) =  u1*R/cv;
                    mat(2, 1) =    -R/cv;
                    
                    
                    mat(n1-1, 0) = (-cv+R)*u1*u2/cv;
                    mat(n1-1, 1) = -u2*R/cv;
                    mat(n1-1, 2) = u1;
                }
                    break;
                    
                case VEL2:
                {
                    switch (dim)
                    {
                        case 3:
                        {
                            mat(3, 0) = -u3;
                            mat(3, 3) = 1.0;
                            
                            mat(n1-1, 3) = -u3*R/cv;
                        }
                            
                        case 1:
                        case 2:
                        {
                            mat(1, 0) = -u1;
                            mat(1, 1) = 1.0;
                            
                            mat(2, 0) = -u2*(2.0-R/cv);
                            mat(2, 2) =     (2.0-R/cv);
                            
                            mat(n1-1, 0) = (-e_tot*(cv+R)+2.0*R*k + u2*u2*(-cv+R))/cv;
                            mat(n1-1, 1) = -u1*R/cv;
                            mat(n1-1, 2) = u2*(1.0-2.0*R/cv);
                            mat(n1-1, n1-1) = (cv+R)/cv;
                        }
                            break;
                    }
                }
                    break;
                    
                case VEL3:
                {
                    mat(2, 0) =  u3*R/cv;
                    mat(2, 3) =    -R/cv;
                    
                    mat(3, 0) = -u2;
                    mat(3, 2) = 1.0;
                    
                    mat(n1-1, 0) = (-cv+R)*u2*u3/cv;
                    mat(n1-1, 2) = u3;
                    mat(n1-1, 3) = -u2*R/cv;
                }
                    break;
                    
                case TEMP:
                {
                    // all dimensions have the same format
                    mat(n1-1, 0) = -u2*(cv+R);
                    mat(n1-1, 2) = cv+R;
                }
                    break;
                    
                default:
                    libmesh_assert_msg(false, "Invalid primitive variable number");
                    break;
            }
        }
            break;
            
        case 2:  // Az
        {
            switch (_active_primitive_vars[primitive_var])
            {
                case RHO_PRIM:
                {
                    // nothing to be done for density; matrix is zero
                }
                    break;
                    
                case VEL1:
                {
                    mat(1, 0) = -u3;
                    mat(1, 3) = 1.0;
                    
                    mat(3, 0) = u1*R/cv;
                    mat(3, 1) =   -R/cv;
                    
                    mat(n1-1, 0) = (-cv+R)*u1*u3/cv;
                    mat(n1-1, 1) = -u3*R/cv;
                    mat(n1-1, 3) =  u1;
                }
                    break;
                    
                case VEL2:
                {
                    mat(2, 0) = -u3;
                    mat(2, 3) = 1.0;
                    
                    mat(3, 0) =  u2*R/cv;
                    mat(3, 2) =    -R/cv;
                    
                    mat(n1-1, 0) = (-cv+R)*u2*u3/cv;
                    mat(n1-1, 2) = -u3*R/cv;
                    mat(n1-1, 3) = u2;
                }
                    break;
                    
                case VEL3:
                {
                    mat(1, 0) = -u1;
                    mat(1, 1) = 1.0;
                    
                    mat(2, 0) = -u2;
                    mat(2, 2) = 1.0;
                    
                    mat(3, 0) =  -u3*(2.0-R/cv);
                    mat(3, 3) =      (2.0-R/cv);
                    
                    mat(n1-1, 0) = (-e_tot*(cv+R)+2.0*R*k + u3*u3*(-cv+R))/cv;
                    mat(n1-1, 1) = -u1*R/cv;
                    mat(n1-1, 2) = -u2*R/cv;
                    mat(n1-1, 3) = u3*(1.0-2.0*R/cv);
                    mat(n1-1, n1-1) = (cv+R)/cv;
                }
                    break;
                    
                case TEMP:
                {
                    // all dimensions have the same format
                    mat(n1-1, 0) = -u3*(cv+R);
                    mat(n1-1, 3) = cv+R;
                }
                    break;
                    
                default:
                    libmesh_assert_msg(false, "Invalid primitive variable number");
                    break;
            }
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}



void
MAST::FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_for_normal
(const MAST::PrimitiveSolution& sol,
 const libMesh::Point& normal,
 RealVectorX& eig_vals,
 RealMatrixX& l_eig_mat,
 RealMatrixX& l_eig_mat_inv_tr) {
    
    
    const unsigned int n1 = 2 + dim;
    
    eig_vals.setZero(); l_eig_mat.setZero(); l_eig_mat_inv_tr.setZero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a,
    gamma = flight_condition->gas_property.gamma,
    cv = flight_condition->gas_property.cv;
    
    // initialize the values
    switch (dim)
    {
        case 3:
        {
            nz = normal(2);
            u += u3*nz;
        }
            
        case 2:
        {
            ny = normal(1);
            u += u2*ny;
        }
            
        case 1:
        {
            nx = normal(0);
            u += u1*nx;
        }
    }
    
    // select the largest value of surface normal component, so that the appropriate section can be chosen
    if ((fabs(nx)>=fabs(ny)) && (fabs(nx)>=fabs(nz)))
        dim_for_eig_vec = 1;
    else if ((fabs(ny)>=fabs(nx)) && (fabs(ny)>=fabs(nz)))
        dim_for_eig_vec = 2;
    else if ((fabs(nz)>=fabs(nx)) && (fabs(nz)>=fabs(ny)))
        dim_for_eig_vec = 3;
    
    // set eigenvalues
    switch (dim)
    {
        case 3:
            eig_vals(2) = u;
        case 2:
            eig_vals(1) = u;
        case 1:
        {
            eig_vals(0) = u;
            eig_vals(n1-2) = u-a;
            eig_vals(n1-1) = u+a;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u-a
            l_eig_mat(3, n1-2) = -u3-nz*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(3, n1-1) = -u3+nz*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(3, 0) = (gamma-1.0)*u1*u3/(a*a)-nx*nz;
            
            // for u
            l_eig_mat_inv_tr(3, 1) = (gamma-1.0)*u2*u3/(a*a)-ny*nz;
            
            // for u
            l_eig_mat_inv_tr(0, 2) = (gamma-1.0)*u3/(a*a);
            l_eig_mat_inv_tr(1, 2) = (gamma-1.0)*u3*u1/(a*a)-nz*nx;
            l_eig_mat_inv_tr(2, 2) = (gamma-1.0)*u3*u2/(a*a)-nz*ny;
            l_eig_mat_inv_tr(3, 2) = (gamma-1.0)*u3*u3/(a*a)+1.0-nz*nz;
            l_eig_mat_inv_tr(n1-1, 2) = (gamma-1.0)*k*u3/(a*a)+u3-nz*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = u3;
            
            // for u-a
            l_eig_mat_inv_tr(3, n1-2) = 0.5*(gamma-1.0)*(-nz+u3/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(3, n1-1) = 0.5*(gamma-1.0)*(nz+u3/a)/a;
            
        }
            
        case 2:
        {
            // for u-a
            l_eig_mat(2, n1-2) = -u2-ny*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(2, n1-1) = -u2+ny*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(2, 0) = (gamma-1.0)*u1*u2/(a*a)-nx*ny;
            
            // for u
            l_eig_mat_inv_tr(0, 1) = (gamma-1.0)*u2/(a*a);
            l_eig_mat_inv_tr(1, 1) = (gamma-1.0)*u2*u1/(a*a)-ny*nx;
            l_eig_mat_inv_tr(2, 1) = (gamma-1.0)*u2*u2/(a*a)+1.0-ny*ny;
            l_eig_mat_inv_tr(n1-1, 1) = (gamma-1.0)*k*u2/(a*a)+u2-ny*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = u2;
            
            // for u-a
            l_eig_mat_inv_tr(2, n1-2) = 0.5*(gamma-1.0)*(-ny+u2/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(2, n1-1) = 0.5*(gamma-1.0)*(ny+u2/a)/a;
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = k+u*a/(gamma-1.0);
            l_eig_mat(1, n1-2) = -u1-nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-2) = 1.0;
            
            // for u+a
            l_eig_mat(0, n1-1) = k-u*a/(gamma-1.0);
            l_eig_mat(1, n1-1) = -u1+nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-1) = 1.0;
            
            // for u
            l_eig_mat_inv_tr(0, 0) = (gamma-1.0)*u1/(a*a);
            l_eig_mat_inv_tr(1, 0) = (gamma-1.0)*u1*u1/(a*a)+1.0-nx*nx;
            l_eig_mat_inv_tr(n1-1, 0) = (gamma-1.0)*k*u1/(a*a)+u1-nx*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 1.0;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = u1;
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = k;
            l_eig_mat_inv_tr.col(dim_for_eig_vec-1) *= -1.0/(cv*T*gamma);
            
            // for u-a
            l_eig_mat_inv_tr(0, n1-2) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-2) = 0.5*(gamma-1.0)*(-nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-2) = 0.5*(1.0+(gamma-1.0)*(-u+k/a)/a);
            
            // for u+a
            l_eig_mat_inv_tr(0, n1-1) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-1) = 0.5*(gamma-1.0)*(nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-1) = 0.5*(1.0+(gamma-1.0)*(u+k/a)/a);
            
        }
    }
    
    
    switch (dim_for_eig_vec)
    {
        case 1:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u1/nx-u3;
                    l_eig_mat(1, 2) = -nz/nx;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = ny*u1/nx-u2;
                    l_eig_mat(1, 1) = -ny/nx;
                    l_eig_mat(2, 1) = 1.0;
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = u*u1/nx-(cv*T*gamma+k);
                    l_eig_mat(1, 0) = -u/nx;
                    l_eig_mat(n1-1, 0) = 1.0;
                }
            }
        }
            break;
            
        case 2:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u2/ny-u3;
                    l_eig_mat(2, 2) = -nz/ny;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = nx*u2/ny-u1;
                    l_eig_mat(1, 0) = 1.0;
                    l_eig_mat(2, 0) = -nx/ny;
                    
                    // for u
                    l_eig_mat(0, 1) = u*u2/ny-(cv*T*gamma+k);
                    l_eig_mat(2, 1) = -u/ny;
                    l_eig_mat(n1-1, 1) = 1.0;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 0) = nx*u3/nz-u1;
            l_eig_mat(1, 0) = 1.0;
            l_eig_mat(3, 0) = -nx/nz;
            
            // for u
            l_eig_mat(0, 1) = ny*u3/nz-u2;
            l_eig_mat(2, 1) = 1.0;
            l_eig_mat(3, 1) = -ny/nz;
            
            // for u
            l_eig_mat(0, 2) = u*u3/nz-(cv*T*gamma+k);
            l_eig_mat(3, 2) = -u/nz;
            l_eig_mat(n1-1, 2) = 1.0;
        }
            break;
    }
}



void
MAST::FluidElemBase::
calculate_pressure_derivative_wrt_conservative_variables(const MAST::PrimitiveSolution& sol,
                                                         RealVectorX& dpdX) {
    
    const unsigned int n1 = 2 + dim;
    
    dpdX.setZero();
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    R = flight_condition->gas_property.R,
    cv= flight_condition->gas_property.cv,
    e = sol.e_tot,
    p = sol.p;
    
    switch (dim)
    {
        case 3:
            dpdX(3)        = -R/cv*u3;
            
        case 2:
            dpdX(2)        = -R/cv*u2;
            
        case 1: {
            
            dpdX(0)        =  R/cv*k;
            dpdX(1)        = -R/cv*u1;
            dpdX(n1-1)     =  R/cv;
        }
    }
}




void
MAST::FluidElemBase::
calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
(const MAST::PrimitiveSolution& sol,
 const Real ui_ni,
 const libMesh::Point& nvec,
 const RealVectorX& dnormal,
 RealMatrixX& mat) {
    
    // the solid wall boundary flux is defined as
    //     fn     = ui_ni (U_c + {0,0,0,0,p}^T) + {0, p n_hat, 0}
    // then,
    //     dfn/dU = ui_ni (I + {0,0,0,0,dp/dU}^T) +
    //               dui_ni/dU (U_c + {0,0,0,0,p}^T) +
    //               {0, dp/dU n1, dp/dU n2, dp/dU n3, 0}
    //
    // and using  ui_ni = wdot (ni + dni) - ui dni
    //      dui_ni/dU = - dni (dui/dU)
    //
    // the first and third terms of dfn/dU are calculated first.
    // the last section calculates the contribution of the second term
    // dfn/dU
    
    
    const unsigned int n1 = 2 + dim;
    
    mat.setZero();
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    R = flight_condition->gas_property.R,
    cv= flight_condition->gas_property.cv,
    e = sol.e_tot,
    p = sol.p;
    
    switch (dim)
    {
        case 3:
        {
            mat(1, 3)        = -R/cv*nvec(0)*u3;         //               dp/dU n1
            
            mat(2, 3)        = -R/cv*nvec(1)*u3;         //               dp/dU n2
            
            mat(3, 0)        =  nvec(2)*R/cv*k;          //               dp/dU n3
            mat(3, 1)        = -R/cv*nvec(2)*u1;         //               dp/dU n3
            mat(3, 2)        = -R/cv*nvec(2)*u2;         //               dp/dU n3
            mat(3, 3)        =  ui_ni-R/cv*nvec(2)*u3;   // ui_ni dU/dU + dp/dU n3
            mat(3, n1-1)     =  R/cv*nvec(2);            //               dp/dU n3
            
            mat(n1-1, 3)     = -ui_ni*R/cv*u3;           // ui_ni dp/dU
        }
            
        case 2:
        {
            mat(1, 2)        = -R/cv*nvec(0)*u2;         //               dp/dU n1
            
            mat(2, 0)        =  nvec(1)*R/cv*k;          //               dp/dU n2
            mat(2, 1)        = -R/cv*nvec(1)*u1;         //               dp/dU n2
            mat(2, 2)        =  ui_ni-R/cv*nvec(1)*u2;   // ui_ni dU/dU + dp/dU n2
            mat(2, n1-1)     =  R/cv*nvec(1);            //               dp/dU n2
            
            mat(n1-1, 2)     = -ui_ni*R/cv*u2;           // ui_ni dp/dU
        }
            
        case 1:
        {
            mat(0, 0)        =  ui_ni;                   // ui_ni dU/dU
            
            mat(1, 0)        =  nvec(0)*R/cv*k;          //               dp/dU n1
            mat(1, 1)        =  ui_ni-R/cv*nvec(0)*u1;   // ui_ni dU/dU + dp/dU n1
            mat(1, n1-1)     =  R/cv*nvec(0);            //               dp/dU n1
            
            mat(n1-1, 0)     =  ui_ni*R/cv*k;            // ui_ni dp/dU
            mat(n1-1, 1)     = -ui_ni*R/cv*u1;           // ui_ni dp/dU
            mat(n1-1, n1-1)  =  ui_ni*(1.+R/cv);         // ui_ni dp/dU
        }
            break;
    }
    
    
    //
    // this calculates the Jacobian contribution from
    // -dni (dui/dU) (U_c + {0,0,0,0,p}^T)
    //
    //
    // note that
    //     ui * ni = wi_dot * (ni + dni) - (ui * dni)
    // hence,
    //     d (ui * ni)/d U = - d ui / dU * dni
    //
    // from primitive variable Jacobian,
    //     du1/dU  =  1/rho {-u1, 1, 0, 0, 0}
    //     du2/dU  =  1/rho {-u2, 0, 1, 0, 0}
    //     du3/dU  =  1/rho {-u3, 0, 0, 1, 0}
    //
    // Hence,
    //    dni dui/dU = 1/rho{-dn1 u1 - dn2 u2 - dn3 u3, dn1, dn2 dn3, 0}
    //
    //  and
    //    (U_c + {0,0,0,0,p}^T) dni dui/dU =
    // {rho, rho u1, rho u2, rho u3, rho E + p}^T . 1/rho{-dn1 u1 - dn2 u2 - dn3 u3, dn1, dn2 dn3, 0}
    //    = {1, u1, u2, u3, E+p/rho} {-dn1 u1 - dn2 u2 - dn3 u3, dn1, dn2 dn3, 0}
    //
    if (dnormal.norm() > 0.)
    {
        RealVectorX
        duini_dU          = RealVectorX::Zero(n1),
        tmp               = RealVectorX::Zero(n1);

        // initialze the Jacobian of ui_ni wrt sol variables
        tmp(0)    = 1.;
        tmp(n1-1) = e+p/rho;
        for (unsigned int i_dim=0; i_dim<dim; i_dim++) {
            
            switch (i_dim) {
                case 0: {
                    
                    duini_dU(0) -=  u1 * dnormal(0);
                    tmp(1)       =  u1;
                }
                    break;
                case 1: {
                    
                    duini_dU(0) -=  u2 * dnormal(1);
                    tmp(2)       =  u2;
                }
                    break;
                case 2: {
                    
                    duini_dU(0) -=  u3 * dnormal(2);
                    tmp(3)       =  u3;
                }
                    break;
            }
            duini_dU(i_dim+1) = dnormal(i_dim);
        }
        
        // now add to the Jacobian matrix
        for (unsigned int i=0; i<n1; i++)
            for (unsigned int j=0; j<n1; j++)
                // look at the note above to see reason for negative sign
                mat(i,j) -= tmp(i)*duini_dU(j);
    }
}




void
MAST::FluidElemBase::
calculate_entropy_variable_jacobian(const MAST::PrimitiveSolution& sol,
                                    RealMatrixX& dUdV,
                                    RealMatrixX& dVdU) {
    
    // calculates dU/dV where V is the Entropy variable vector
    
    // calculate A0 = d U / d Y , where U = conservation variable vector,
    // and Y = unknown variable vector
    // note that for conservation variables as the unknown, this is an
    // identity matrix
    
    const unsigned int n1 = 2 + dim;
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    e_tot = sol.e_tot,
    gamma = flight_condition->gas_property.gamma,
    cv = flight_condition->gas_property.cv;
    
    dUdV.setZero(); dVdU.setZero();
    
    // du/dv
    switch (dim)
    {
        case 3:
        {
            dUdV(0, 3) = u3;
            
            dUdV(1, 3) = u1*u3;
            
            dUdV(2, 3) = u2*u3;
            
            dUdV(3, 0) = dUdV(0, 3);
            dUdV(3, 1) = dUdV(1, 3);
            dUdV(3, 2) = dUdV(2, 3);
            dUdV(3, 3) = u3*u3+cv*T*(gamma-1.0);
            dUdV(3, n1-1) = u3*(cv*T*gamma+k);
            
            dUdV(n1-1, 3) = dUdV(3, n1-1);
        }
            
        case 2:
        {
            dUdV(0, 2) = u2;
            
            dUdV(1, 2) = u1*u2;
            
            dUdV(2, 0) = dUdV(0, 2);
            dUdV(2, 1) = dUdV(1, 2);
            dUdV(2, 2) = u2*u2+cv*T*(gamma-1.0);
            dUdV(2, n1-1) = u2*(cv*T*gamma+k);
            
            dUdV(n1-1, 2) = dUdV(2, n1-1);
        }
            
        case 1:
        {
            dUdV(0, 0) = 1.0;
            dUdV(0, 1) = u1;
            dUdV(0, n1-1) = e_tot;
            
            dUdV(1, 0) = dUdV(0, 1);
            dUdV(1, 1) = u1*u1+cv*T*(gamma-1.0);
            dUdV(1, n1-1) = u1*(cv*T*gamma+k);
            
            dUdV(n1-1, 0) = dUdV(0, n1-1);
            dUdV(n1-1, 1) = dUdV(1, n1-1);
            dUdV(n1-1, n1-1) = k*k+gamma*cv*T*(cv*T+2*k);
            
        }
            break;
            
        default:
            break;
    }
    
    dUdV *= rho/(gamma-1.0);
    
    
    // dv/du
    switch (dim)
    {
        case 3:
        {
            dVdU(0, 3) = -u3*k;
            
            dVdU(1, 3) = u1*u3;
            
            dVdU(2, 3) = u2*u3;
            
            dVdU(3, 0) = dVdU(0, 3);
            dVdU(3, 1) = dVdU(1, 3);
            dVdU(3, 2) = dVdU(2, 3);
            dVdU(3, 3) = u3*u3+cv*T;
            dVdU(3, n1-1) = -u3;
            
            dVdU(n1-1, 3) = dVdU(3, n1-1);
        }
            
        case 2:
        {
            dVdU(0, 2) = -u2*k;
            
            dVdU(1, 2) = u1*u2;
            
            dVdU(2, 0) = dVdU(0, 2);
            dVdU(2, 1) = dVdU(1, 2);
            dVdU(2, 2) = u2*u2+cv*T;
            dVdU(2, n1-1) = -u2;
            
            dVdU(n1-1, 2) = dVdU(2, n1-1);
        }
            
        case 1:
        {
            dVdU(0, 0) = k*k+cv*cv*T*T*gamma;
            dVdU(0, 1) = -u1*k;
            dVdU(0, n1-1) = -e_tot+2.0*k;
            
            dVdU(1, 0) = dVdU(0, 1);
            dVdU(1, 1) = u1*u1+cv*T;
            dVdU(1, n1-1) = -u1;
            
            dVdU(n1-1, 0) = dVdU(0, n1-1);
            dVdU(n1-1, 1) = dVdU(1, n1-1);
            dVdU(n1-1, n1-1) = 1.0;
        }
            break;
            
        default:
            break;
    }
    
    dVdU /= (rho*cv*cv*T*T);
}




bool
MAST::FluidElemBase::
calculate_barth_tau_matrix (const unsigned int qp,
                            const libMesh::FEBase& fe,
                            const MAST::PrimitiveSolution& sol,
                            RealMatrixX& tau,
                            std::vector<RealMatrixX >& tau_sens) {
    
    const unsigned int n1 = 2 + dim;
    
    libMesh::Point nvec;
    RealVectorX
    eig_val               = RealVectorX::Zero(n1);
    
    RealMatrixX
    l_eig_vec             = RealMatrixX::Zero(n1, n1),
    l_eig_vec_inv_tr      = RealMatrixX::Zero(n1, n1),
    tmp1                  = RealMatrixX::Zero(n1, n1);
    
    Real nval = 0.;
    
    const std::vector<std::vector<libMesh::RealVectorValue> >&
    dphi = fe.get_dphi(); // assuming that all variables have the same interpolation
    
    for (unsigned int i_node=0; i_node<dphi.size(); i_node++)
    {
        nvec = dphi[i_node][qp];
        nval = nvec.size();
        if (nval > 0.) {
            
            nvec /= nval;
            this->calculate_advection_left_eigenvector_and_inverse_for_normal
            (sol, nvec, eig_val, l_eig_vec, l_eig_vec_inv_tr);
            
            for (unsigned int i_var=0; i_var<n1; i_var++)
                l_eig_vec_inv_tr.col(i_var) *= fabs(eig_val(i_var)); // L^-T [omaga]
            
            l_eig_vec_inv_tr *= l_eig_vec.transpose(); // A = L^-T [omaga] L^T
            
            tmp1 += nval * l_eig_vec_inv_tr;  // sum_inode  | A_i |
        }
    }
    
    
    // now invert the tmp matrix to get the tau matrix
    RealVectorX
    x    =   RealVectorX::Zero(n1),
    b    =   RealVectorX::Zero(n1);

    for (unsigned int i_var=0; i_var<n1; i_var++) {
        
        tau_sens[i_var].setZero(); // zero the sensitivity matrix for now
        b.setZero(); b(i_var) = 1.;
        x = tmp1.lu().solve(b);
        tau.col(i_var) = x;
    }
    
    return false;
}



bool
MAST::FluidElemBase::
calculate_aliabadi_tau_matrix (const unsigned int qp,
                               const libMesh::FEBase& fe,
                               const MAST::PrimitiveSolution& sol,
                               RealMatrixX& tau,
                               std::vector<RealMatrixX >& tau_sens) {
    
    const unsigned int n1 = 2 + dim;
    
    const std::vector<std::vector<libMesh::RealVectorValue> >& dphi =
    fe.get_dphi(); // assuming that all variables have the same interpolation
    
    RealVectorX
    u                = RealVectorX::Zero(dim),
    dN               = RealVectorX::Zero(dim);

    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    a = sol.a,
    dt = 0., //c.get_deltat_value(),
    cp = flight_condition->gas_property.cp,
    R = flight_condition->gas_property.R,
    gamma = flight_condition->gas_property.gamma;
    
    tau.setZero();
    
    // calculate the gradients
    switch (dim)
    {
        case 3:
            u(2) = u3;
            
        case 2:
            u(1) = u2;
            
        case 1:
            u(0) = u1;
            break;
            
        default:
            break;
    }
    
    // calculate the dot product of velocity times gradient of shape function
    Real h = 0, u_val = u.norm(), tau_rho, tau_m, tau_e;
    u /= u_val;
    
    for (unsigned int i_nodes=0; i_nodes<dphi.size(); i_nodes++)
    {
        // set value of shape function gradient
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            dN(i_dim) = dphi[i_nodes][qp](i_dim);
        
        h += fabs(dN.dot(u));
    }
    
    h = 2.0/h;
    
    // now set the value of streamwise dissipation
    tau_rho = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    if (!_if_viscous)
    {
        tau_m   = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
        tau_e   = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    }
    else
    {
        tau_m   = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2) +
                           pow(4.0*sol.mu/sol.rho/h/h, 2));
        tau_e   = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2)+
                           pow(4.0*sol.k_thermal/sol.rho/cp/h/h, 2));
    }
    
    tau(0, 0) = tau_rho;
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        tau(1+i_dim, 1+i_dim) = tau_m;
    tau(n1-1, n1-1) = tau_e;
    
    // calculation sensitivity of the tau matrix for each conservative variable
    std::vector<Real> primitive_sens(n1), cons_sens(n1);
    
    // sensitivity wrt primitive variables
    for (unsigned int i_pvar=0; i_pvar<n1; i_pvar++)
    {
        switch (_active_primitive_vars[i_pvar])
        {
            case RHO_PRIM:
                primitive_sens[i_pvar] = 0.;
                break;
                
            case VEL1:
                primitive_sens[i_pvar] = 2.0*u1/sqrt(2.0*sol.k)/h;
                break;
                
            case VEL2:
                primitive_sens[i_pvar] = 2.0*u2/sqrt(2.0*sol.k)/h;
                break;
                
            case VEL3:
                primitive_sens[i_pvar] = 2.0*u3/sqrt(2.0*sol.k)/h;
                break;
                
            case TEMP:
                primitive_sens[i_pvar] = sqrt(gamma*R/sol.T)/h;
                break;
                
            default:
                break;
        }
        
        RealMatrixX
        dprim_dcons      = RealMatrixX::Zero(n1, n1),
        dcons_dprim      = RealMatrixX::Zero(n1, n1);

        std::fill(cons_sens.begin(), cons_sens.end(), 0.);
        
        // apply chain rule
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        {
            for (unsigned int i_pvar=0; i_pvar<n1; i_pvar++)
                cons_sens[i_cvar] += primitive_sens[i_pvar] *
                dprim_dcons(i_pvar, i_cvar);
        }
        
        // set the values in the sensitivity matrices
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        {
            tau_sens[i_cvar].setZero();
            cons_sens[i_cvar] *= -pow(2.0/h*(u_val+a),-2.0);
            for (unsigned int i=0; i<n1; i++)
                tau_sens[i_cvar](i,i) = cons_sens[i_cvar];
        }
    }
    
    return true;
}



void
MAST::FluidElemBase::
calculate_dxidX (const unsigned int qp, const libMesh::FEBase& fe,
                 RealMatrixX& dxi_dX,
                 RealMatrixX& dX_dxi) {
    
    
    // initialize dxi_dX and dX_dxi
    dxi_dX.setZero(); dX_dxi.setZero();
    Real val=0., val2=0.;
    
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        for (unsigned int j_dim=0; j_dim<dim; j_dim++)
        {
            switch (i_dim)
            {
                case 0:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val  = fe.get_dxidx()[qp];
                            val2 = fe.get_dxyzdxi()[qp](i_dim);
                            break;
                        case 1:
                            val = fe.get_dxidy()[qp];
                            val2 = fe.get_dxyzdeta()[qp](i_dim);
                            break;
                        case 2:
                            val = fe.get_dxidz()[qp];
                            val2 = fe.get_dxyzdzeta()[qp](i_dim);
                            break;
                    }
                }
                    break;
                case 1:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe.get_detadx()[qp];
                            val2 = fe.get_dxyzdxi()[qp](i_dim);
                            break;
                        case 1:
                            val = fe.get_detady()[qp];
                            val2 = fe.get_dxyzdeta()[qp](i_dim);
                            break;
                        case 2:
                            val = fe.get_detadz()[qp];
                            val2 = fe.get_dxyzdzeta()[qp](i_dim);
                            break;
                    }
                }
                    break;
                case 2:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe.get_dzetadx()[qp];
                            val2 = fe.get_dxyzdxi()[qp](i_dim);
                            break;
                        case 1:
                            val = fe.get_dzetady()[qp];
                            val2 = fe.get_dxyzdeta()[qp](i_dim);
                            break;
                        case 2:
                            val = fe.get_dzetadz()[qp];
                            val2 = fe.get_dxyzdzeta()[qp](i_dim);
                            break;
                    }
                }
                    break;
            }
            dxi_dX(i_dim, j_dim) = val;
            dX_dxi(i_dim, j_dim) = val2;
        }
}




void MAST::FluidElemBase::
calculate_hartmann_discontinuity_operator (const unsigned int qp,
                                           const libMesh::FEBase& fe,
                                           const MAST::PrimitiveSolution& sol,
                                           const RealVectorX& elem_solution,
                                           const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                                           const RealMatrixX& Ai_Bi_advection,
                                           RealVectorX& discontinuity_val) {
    
    discontinuity_val.setZero();
    const unsigned int n1 = 2 + dim;
    
    RealMatrixX
    tmpmat1_n1n1       = RealMatrixX::Zero(dim+2, dim+2),
    dxi_dX             = RealMatrixX::Zero(dim, dim),
    dX_dxi             = RealMatrixX::Zero(dim, dim),
    dpdc               = RealMatrixX::Zero(n1, n1),
    dcdp               = RealMatrixX::Zero(n1, n1);
    RealVectorX
    vec1               = RealVectorX::Zero(n1),
    vec2               = RealVectorX::Zero(n1),
    dpress_dp          = RealVectorX::Zero(n1),
    dp                 = RealVectorX::Zero(dim),
    hk                 = RealVectorX::Zero(dim);

    
    // residual of the strong form of the equation: assuming steady flow
    vec2            = Ai_Bi_advection * elem_solution;// Ai dU/dxi
    
    // calculate dp/dU
    calculate_conservative_variable_jacobian(sol, dcdp, dpdc);
    dpress_dp(0)    = (sol.cp - sol.cv)*sol.T; // R T
    dpress_dp(n1-1) = (sol.cp - sol.cv)*sol.rho; // R rho
    vec1            = dpdc.transpose() * dpress_dp; // dpress/dprimitive * dprimitive/dconservative
    
    Real rp = fabs(vec2.dot(vec1)) / sol.p;  // sum_i=1..m  (dp/dU_i)*R_i / p
    
    // calculate the approximate element size
    this->calculate_dxidX (qp, fe, dxi_dX, dX_dxi);
    for (unsigned int i=0; i<dim; i++) {
        
        for (unsigned int j=0; j<dim; j++)
            hk(i) = fmax(hk(i), fabs(dX_dxi(i, j))); // get the mqximum value out of each xi
        
        // calculate gradient of p = dp/dprim * dprim/dcons * dcons/dx_i
        dB_mat[i].vector_mult(vec1, elem_solution);
        vec2 = dpdc * vec1;
        dp(i) = vec2.dot(dpress_dp);
    }
    
    const unsigned int fe_order = fe.get_fe_type().order;
    
    // now set the value of the dissipation coefficient
    for (unsigned int i=0; i<dim; i++)
        discontinuity_val(i) =
        (dp.norm() * hk(i) / sol.p / (fe_order + 1)) * // pressure switch
        pow(hk(i), 2) * rp;
    discontinuity_val *= _dissipation_scaling;
}




void MAST::FluidElemBase::
calculate_aliabadi_discontinuity_operator(const unsigned int qp,
                                          const libMesh::FEBase& fe,
                                          const MAST::PrimitiveSolution& sol,
                                          const RealVectorX& elem_solution,
                                          const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                                          const RealMatrixX& Ai_Bi_advection,
                                          RealVectorX& discontinuity_val) {
    
    
    discontinuity_val.setZero();
    const unsigned int n1 = 2 + dim;
    
    std::vector<RealVectorX >
    diff_vec(3);
    RealMatrixX
    A_inv_entropy      = RealMatrixX::Zero(dim+2, dim+2),
    A_entropy          = RealMatrixX::Zero(dim+2, dim+2),
    dxi_dX             = RealMatrixX::Zero(dim, dim),
    dX_dxi             = RealMatrixX::Zero(dim, dim),
    tmpmat1_n1n1       = RealMatrixX::Zero(n1, n1);
    RealVectorX
    vec1               = RealVectorX::Zero(n1),
    vec2               = RealVectorX::Zero(n1);

    for (unsigned int i=0; i<dim; i++) diff_vec[i].setZero(n1);
    
    Real dval;
    
    this->calculate_dxidX (qp, fe, dxi_dX, dX_dxi);
    this->calculate_entropy_variable_jacobian ( sol, A_entropy, A_inv_entropy );
    
    
    for (unsigned int i=0; i<dim; i++)
        dB_mat[i].vector_mult(diff_vec[i], elem_solution); // dU/dxi
    vec1 = Ai_Bi_advection * elem_solution; // Ai dU/dxi
    
    // TODO: divergence of diffusive flux
    
    // add the velocity and calculate the numerator of the discontinuity
    // capturing term coefficient
    //vec2 += c.elem_solution; // add velocity TODO: how to get the
    vec2 = A_inv_entropy * vec1;
    dval = vec1.dot(vec2);  // this is the numerator term
    
    // now evaluate the dissipation factor for the discontinuity capturing term
    // this is the denominator term
    
    // add a small number to avoid division of zero by zero
    Real val1 = 1.0e-6;
    for (unsigned int i=0; i<dim; i++) {
        vec1.setZero();
        
        for (unsigned int j=0; j<dim; j++)
            vec1 += dxi_dX(i, j) * diff_vec[j];
        
        // calculate the value of denominator
        vec2 = A_inv_entropy * vec1;
        val1 += vec1.dot(vec2);
    }
    
    /*    // now calculate the Ducros shock sensor
     //    Real d_ducros = 0., div = 0.;
     //    RealVectorX u, dudx, dudy, dudz, dN, curl;
     //    u.resize(dim); dudx.resize(3); dudy.resize(3); dudz.resize(3); curl.resize(3);
     //    sol.get_uvec(u);
     //
     //    for (unsigned int i=0; i<dim; i++) {
     //        dudx(i) = diff_vec[0](i+1); // du/dx
     //        dudx(i) -= diff_vec[0](0)*u(i); // -drho/dx * u
     //        dudx(i) /= sol.rho;
     //    }
     //    div += dudx(0);
     //    if (dim > 1) {
     //        for (unsigned int i=0; i<dim; i++) {
     //            dudy(i) = diff_vec[1](i+2); // du/dy
     //            dudy(i) -= diff_vec[1](0)*u(i); // -drho/dy * u
     //            dudy(i) /= sol.rho;
     //        }
     //        div += dudy(1);
     //    }
     //    if (dim > 2) {
     //        for (unsigned int i=0; i<dim; i++) {
     //            dudz(i) = diff_vec[2](i+3); // du/dz
     //            dudz(i) -= diff_vec[2](0)*u(i); // -drho/dz * u
     //            dudz(i) /= sol.rho;
     //        }
     //        div += dudz(2);
     //    }
     //    curl(0) =   dudy(2)-dudz(1);
     //    curl(1) = -(dudx(2)-dudz(0));
     //    curl(2) =   dudx(1)-dudy(0);
     //    d_ducros = (div*div) / (div*div + pow(curl.l2_norm(),2) + 1.0e-6); */
    
    dval = sqrt(dval/val1);
    
    if (_include_pressure_switch) {
        // also add a pressure switch q
        RealMatrixX
        dpdc      = RealMatrixX::Zero(n1, n1),
        dcdp      = RealMatrixX::Zero(n1, n1);
        
        RealVectorX
        dpress_dp          = RealVectorX::Zero(n1, n1),
        dp                 = RealVectorX::Zero(dim);

        Real p_sensor = 0., hk = 0.;
        calculate_conservative_variable_jacobian(sol, dcdp, dpdc);
        dpress_dp(0) = (sol.cp - sol.cv)*sol.T; // R T
        dpress_dp(n1-1) = (sol.cp - sol.cv)*sol.rho; // R rho
        for (unsigned int i=0; i<dim; i++) {
            dB_mat[i].vector_mult(vec1, elem_solution);
            vec2 = dpdc * vec1;
            dp(i) = vec2.dot(dpress_dp);
            for (unsigned int j=0; j<dim; j++)
                hk = fmax(hk, fabs(dX_dxi(i, j)));
        }
        p_sensor = dp.norm() * hk / sol.p;
        dval *= (p_sensor * _dissipation_scaling);
    }
    
    dval *= _dissipation_scaling;
    
    // set value in all three dimensions to be the same
    const unsigned int fe_order = fe.get_fe_type().order;
    for (unsigned int i=0; i<dim; i++)
        discontinuity_val(i) = dval/fe_order;
}




template <typename ValType>
void
MAST::FluidElemBase::
calculate_small_disturbance_aliabadi_discontinuity_operator
(const unsigned int qp,
 const libMesh::FEBase& fe,
 const MAST::PrimitiveSolution& sol,
 const SmallPerturbationPrimitiveSolution<ValType>& dsol,
 const RealVectorX& elem_solution,
 const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
 const RealMatrixX& Ai_Bi_advection,
 RealVectorX& discontinuity_val) {
    
    
    discontinuity_val.setZero();
    const unsigned int n1 = 2 + dim;
    
    std::vector<RealVectorX >
    diff_vec(3);
    
    RealMatrixX
    A_inv_entropy      = RealMatrixX::Zero(dim+2, dim+2),
    A_entropy          = RealMatrixX::Zero(dim+2, dim+2),
    dxi_dX             = RealMatrixX::Zero(dim, dim),
    dX_dxi             = RealMatrixX::Zero(dim, dim),
    tmpmat1_n1n1       = RealMatrixX::Zero(dim+2, dim+2);
    RealVectorX
    vec1               = RealVectorX::Zero(n1),
    vec2               = RealVectorX::Zero(n1);

    for (unsigned int i=0; i<dim; i++) diff_vec[i].setZero(n1);
    
    Real dval;
    
    this->calculate_dxidX (qp, fe, dxi_dX, dX_dxi);
    this->calculate_entropy_variable_jacobian ( sol, A_entropy, A_inv_entropy );
    
    
    for (unsigned int i=0; i<dim; i++)
        dB_mat[i].vector_mult(diff_vec[i], elem_solution); // dU/dxi
    vec1 = Ai_Bi_advection * elem_solution; // Ai dU/dxi
    
    // TODO: divergence of diffusive flux
    
    // add the velocity and calculate the numerator of the discontinuity
    // capturing term coefficient
    //vec2 += c.elem_solution; // add velocity TODO: how to get the
    vec2 = A_inv_entropy * vec1;
    dval = vec1.dot(vec2);  // this is the numerator term
    
    // now evaluate the dissipation factor for the discontinuity capturing term
    // this is the denominator term
    
    // add a small number to avoid division of zero by zero
    Real val1 = 1.0e-6;
    for (unsigned int i=0; i<dim; i++) {
        vec1.setZero();
        
        for (unsigned int j=0; j<dim; j++)
            vec1 += dxi_dX(i, j) * diff_vec[j];
        
        // calculate the value of denominator
        vec2 = A_inv_entropy * vec1;
        val1 += vec1.dot(vec2);
    }
    
    dval = std::min(1., sqrt(dval/val1));
    
    if (_include_pressure_switch) {
        // also add a pressure switch q
        RealMatrixX
        dpdc                    = RealMatrixX::Zero(n1, n1),
        dcdp                    = RealMatrixX::Zero(n1, n1);
        RealVectorX
        dpress_dp               = RealVectorX::Zero(n1),
        dp                      = RealVectorX::Zero(dim);

        Real p_sensor = 0., hk = 0.;
        calculate_conservative_variable_jacobian(sol, dcdp, dpdc);
        dpress_dp(0) = (sol.cp - sol.cv)*sol.T; // R T
        dpress_dp(n1-1) = (sol.cp - sol.cv)*sol.rho; // R rho
        for (unsigned int i=0; i<dim; i++) {
            dB_mat[i].vector_mult(vec1, elem_solution);
            vec2 = dpdc * vec1;
            dp(i) = vec2.dot(dpress_dp);
            for (unsigned int j=0; j<dim; j++)
                hk = fmax(hk, fabs(dX_dxi(i, j)));
        }
        p_sensor = std::min(1.,
                            dp.norm() * hk /
                            (sol.p + abs(dsol.dp)));
        dval *= (p_sensor * _dissipation_scaling);
    }
    
    dval *= _dissipation_scaling;
    
    // set value in all three dimensions to be the same
    const unsigned int fe_order = fe.get_fe_type().order;
    for (unsigned int i=0; i<dim; i++)
        discontinuity_val(i) = dval/fe_order;
}



void MAST::FluidElemBase::
calculate_differential_operator_matrix (const unsigned int qp,
                                        const libMesh::FEBase& fe,
                                        const RealVectorX& elem_solution,
                                        const MAST::PrimitiveSolution& sol,
                                        const MAST::FEMOperatorMatrix& B_mat,
                                        const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
                                        const std::vector<RealMatrixX >& Ai_advection,
                                        const RealMatrixX& Ai_Bi_advection,
                                        const std::vector<std::vector<RealMatrixX > >& Ai_sens,
                                        RealMatrixX& LS_operator,
                                        RealMatrixX& LS_sens) {
    
    const unsigned int n1 = 2 + dim, n2 = B_mat.n();
    
    RealMatrixX
    mat                = RealMatrixX::Zero(n1, n2),
    mat2               = RealMatrixX::Zero(n1, n2),
    tau                = RealMatrixX::Zero(n1, n1);
    RealVectorX
    vec1               = RealVectorX::Zero(n1),
    vec2               = RealVectorX::Zero(n1),
    vec3               = RealVectorX::Zero(n1),
    vec4_n2            = RealVectorX::Zero(n2);

    
    const std::vector<std::vector<Real> >& phi =
    fe.get_phi(); // assuming that all variables have the same interpolation
    const unsigned int n_phi = phi.size();
    std::vector<RealMatrixX > tau_sens(n1);
    for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        tau_sens[i_cvar].setZero(n1, n1);
    
    // contribution of unsteady term
    LS_operator.setZero();
    LS_sens.setZero();
    
    bool if_diagonal_tau = false;
    
    vec2.setZero();
    vec2 = Ai_Bi_advection * elem_solution; // sum A_i dU/dx_i
    
    //if_diagonal_tau = this->calculate_aliabadi_tau_matrix
    //(qp, c, sol, tau, tau_sens);
    if_diagonal_tau = this->calculate_barth_tau_matrix
    (qp, fe, sol, tau, tau_sens);
    
    // contribution of advection flux term
    for (unsigned int i=0; i<dim; i++)
    {
        mat = Ai_advection[i].transpose();
        dB_mat[i].left_multiply(mat2, mat);
        LS_operator += mat2;  // A_i^T dB/dx_i
        
        
        // sensitivity of the LS operator times strong form of residual
        // Bi^T dAi/dalpha tau
        vec1 = tau * vec2;
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        {
            vec3 = Ai_sens[i][i_cvar] * vec1;
            dB_mat[i].vector_mult_transpose(vec4_n2, vec3);
            for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                LS_sens.col((n_phi*i_cvar)+i_phi) += phi[i_phi][qp] * vec4_n2;
        }
        
        // Bi^T Ai dtau/dalpha
        for (unsigned int i_cvar=0; i_cvar<n1; i_cvar++)
        {
            vec1 = tau_sens[i_cvar] * vec2;
            vec3 = Ai_advection[i] * vec1;
            dB_mat[i].vector_mult_transpose(vec4_n2, vec3);
            for (unsigned int i_phi=0; i_phi<n_phi; i_phi++)
                LS_sens.col((n_phi*i_cvar)+i_phi) += phi[i_phi][qp] * vec4_n2;
        }
    }
    
    // scale the LS matrix with the correct factor
    if (if_diagonal_tau) {
        
        for (unsigned int i=0; i<n1; i++)
            LS_operator.row(i) *= tau(i,i);
    }
    else
        LS_operator = tau.transpose() * LS_operator;
}



// template instantiations
template void
MAST::FluidElemBase::
calculate_small_disturbance_aliabadi_discontinuity_operator<Real>
(const unsigned int qp,
 const libMesh::FEBase& fe,
 const MAST::PrimitiveSolution& sol,
 const SmallPerturbationPrimitiveSolution<Real>& dsol,
 const RealVectorX& elem_solution,
 const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
 const RealMatrixX& Ai_Bi_advection,
 RealVectorX& discontinuity_val);



template void
MAST::FluidElemBase::
calculate_small_disturbance_aliabadi_discontinuity_operator<Complex>
(const unsigned int qp,
 const libMesh::FEBase& fe,
 const MAST::PrimitiveSolution& sol,
 const SmallPerturbationPrimitiveSolution<Complex>& dsol,
 const RealVectorX& elem_solution,
 const std::vector<MAST::FEMOperatorMatrix>& dB_mat,
 const RealMatrixX& Ai_Bi_advection,
 RealVectorX& discontinuity_val);



