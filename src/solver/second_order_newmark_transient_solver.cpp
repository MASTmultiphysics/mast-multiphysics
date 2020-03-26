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
#include "solver/second_order_newmark_transient_solver.h"
#include "base/transient_assembly_elem_operations.h"
#include "base/elem_base.h"
#include "base/system_initialization.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/numeric_vector.h"


MAST::SecondOrderNewmarkTransientSolver::SecondOrderNewmarkTransientSolver():
MAST::TransientSolverBase(2, 2),
beta(0.25),
gamma(0.5)
{ }


MAST::SecondOrderNewmarkTransientSolver::~SecondOrderNewmarkTransientSolver()
{ }




void
MAST::SecondOrderNewmarkTransientSolver::
set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                 const std::vector<libMesh::NumericVector<Real>*>& sols){
    
    libmesh_assert_equal_to(sols.size(), 3);

    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol          = RealVectorX::Zero(n_dofs),
    vel          = RealVectorX::Zero(n_dofs),
    accel        = RealVectorX::Zero(n_dofs);
    
    
    // get the references to current and previous sol and velocity
    const libMesh::NumericVector<Real>
    &sol_global     =   *sols[0],
    &vel_global     =   *sols[1],
    &acc_global     =   *sols[2];
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
        accel(i)        = acc_global(dof_indices[i]);
    }
    
    _assembly_ops->set_elem_solution(sol);
    _assembly_ops->set_elem_velocity(vel);
    _assembly_ops->set_elem_acceleration(accel);
}



void
MAST::SecondOrderNewmarkTransientSolver::
extract_element_sensitivity_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                                 const std::vector<libMesh::NumericVector<Real>*>& sols,
                                 std::vector<RealVectorX>& local_sols) {
    
    libmesh_assert_equal_to(sols.size(), 3);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    local_sols.resize(3);
    
    RealVectorX
    &sol         = local_sols[0],
    &vel         = local_sols[1],
    &accel       = local_sols[2];
    
    sol          = RealVectorX::Zero(n_dofs);
    vel          = RealVectorX::Zero(n_dofs);
    accel        = RealVectorX::Zero(n_dofs);
    
    
    // get the references to current and previous sol and velocity
    const libMesh::NumericVector<Real>
    &sol_global     =   *sols[0],
    &vel_global     =   *sols[1],
    &acc_global     =   *sols[2];
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
        accel(i)        = acc_global(dof_indices[i]);
    }
    
    
}



void
MAST::SecondOrderNewmarkTransientSolver::
set_element_perturbed_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                           const std::vector<libMesh::NumericVector<Real>*>& sols){
    
    libmesh_assert_equal_to(sols.size(), 3);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol          = RealVectorX::Zero(n_dofs),
    vel          = RealVectorX::Zero(n_dofs),
    accel        = RealVectorX::Zero(n_dofs);
    
    
    // get the references to current and previous sol and velocity
    const libMesh::NumericVector<Real>
    &sol_global     =   *sols[0],
    &vel_global     =   *sols[1],
    &acc_global     =   *sols[2];
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
        accel(i)        = acc_global(dof_indices[i]);
    }
    
    _assembly_ops->set_elem_perturbed_solution(sol);
    _assembly_ops->set_elem_perturbed_velocity(vel);
    _assembly_ops->set_elem_perturbed_acceleration(accel);
}




void
MAST::SecondOrderNewmarkTransientSolver::
update_velocity(libMesh::NumericVector<Real>& vec,
                const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution(1),
    &prev_vel = this->velocity(1),
    &prev_acc = this->acceleration(1);

    // first calculate the acceleration
    vec = sol;
    vec.add(                  -1., prev_sol);
    vec.scale(gamma/beta/dt);
    vec.add(        1.-gamma/beta, prev_vel);
    vec.add((1.-gamma/2./beta)*dt, prev_acc);
    vec.close();
}




void
MAST::SecondOrderNewmarkTransientSolver::
update_acceleration(libMesh::NumericVector<Real>& vec,
                    const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution(1),
    &prev_vel = this->velocity(1),
    &prev_acc = this->acceleration(1);
    
    // first calculate the acceleration
    vec = sol;
    vec.add(            -1., prev_sol);
    vec.scale(1./beta/dt/dt);
    vec.add(    -1./beta/dt, prev_vel);
    vec.add(-(.5-beta)/beta, prev_acc);
    vec.close();
}



void
MAST::SecondOrderNewmarkTransientSolver::
update_sensitivity_velocity(libMesh::NumericVector<Real>& vec,
                            const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution_sensitivity(1),
    &prev_vel = this->velocity_sensitivity(1),
    &prev_acc = this->acceleration_sensitivity(1);
    
    // first calculate the acceleration
    vec = sol;
    vec.add(                  -1., prev_sol);
    vec.scale(gamma/beta/dt);
    vec.add(        1.-gamma/beta, prev_vel);
    vec.add((1.-gamma/2./beta)*dt, prev_acc);
    vec.close();
}



void
MAST::SecondOrderNewmarkTransientSolver::
update_sensitivity_acceleration(libMesh::NumericVector<Real>& vec,
                                const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution_sensitivity(1),
    &prev_vel = this->velocity_sensitivity(1),
    &prev_acc = this->acceleration_sensitivity(1);
    
    // first calculate the acceleration
    vec = sol;
    vec.add(            -1., prev_sol);
    vec.scale(1./beta/dt/dt);
    vec.add(    -1./beta/dt, prev_vel);
    vec.add(-(.5-beta)/beta, prev_acc);
    vec.close();
}




void
MAST::SecondOrderNewmarkTransientSolver::
update_delta_velocity(libMesh::NumericVector<Real>& vec,
                      const libMesh::NumericVector<Real>& sol) {
    
    // first calculate the acceleration
    vec = sol;
    vec.scale( gamma/beta/dt);
    vec.close();
}




void
MAST::SecondOrderNewmarkTransientSolver::
update_delta_acceleration(libMesh::NumericVector<Real>& vec,
                          const libMesh::NumericVector<Real>& sol) {
    
    // first calculate the acceleration
    vec = sol;
    vec.scale(  1./beta/dt/dt);
    vec.close();
}





void
MAST::SecondOrderNewmarkTransientSolver::
elem_calculations(bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    // make sure that the assembly object is provided
    libmesh_assert(_assembly);
    unsigned int n_dofs = (unsigned int)vec.size();
    
    RealVectorX
    f_x     = RealVectorX::Zero(n_dofs),
    f_m     = RealVectorX::Zero(n_dofs);
    
    RealMatrixX
    f_m_jac_xddot    = RealMatrixX::Zero(n_dofs, n_dofs),
    f_m_jac_xdot     = RealMatrixX::Zero(n_dofs, n_dofs),
    f_m_jac          = RealMatrixX::Zero(n_dofs, n_dofs),
    f_x_jac_xdot     = RealMatrixX::Zero(n_dofs, n_dofs),
    f_x_jac          = RealMatrixX::Zero(n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_calculations(if_jac,
                                     f_m,           // mass vector
                                     f_x,           // forcing vector
                                     f_m_jac_xddot, // Jac of mass wrt x_dotdot
                                     f_m_jac_xdot,  // Jac of mass wrt x_dot
                                     f_m_jac,       // Jac of mass wrt x
                                     f_x_jac_xdot,  // Jac of forcing vector wrt x_dot
                                     f_x_jac);      // Jac of forcing vector wrt x

    
    if (_if_highest_derivative_solution) {
    
        //
        //  The residual here is modeled as
        // r(x, xdot, xddot) = f_m(x, xdot, xddot) + f_x(x, xdot)= 0
        //
        //  then, given x = x0, and xdot = xdot0, the residual can be used to
        //  evaluate xddot at the same time step as
        //
        //  r(x0, xdot0, xddot) + dr/dxddot dxddot = 0
        //
        
        // system residual
        vec  = (f_m + f_x);
        
        // system Jacobian
        if (if_jac)
            mat = f_m_jac_xddot;
    }
    else {
        
        // N-R solver: r(x) = 0, and provides estimates for x (current solution)
        // in an iterative fashion.
        // Newmark solver uses x and evaluates x_ddot and x_dot (acc and vel at
        // current x).
        //
        // N-R requires   r(x) and dr/dx   (residual and Jacobian)
        // N-R asks Newmark to provide these quantities
        //
        // Newmark asks elements (which provide physics) for contributions to
        // f_m and f_x (using the current estimates of x_ddot, x_dot, x)
        //
        //
        //  The residual here is modeled as
        // r = (f_m + f_x )= 0
        // where, (for example)
        // f_m = int_Omega phi u_dot   [typical mass vector in conduction, for example]
        // f_x = int_Omega phi_i u_i - int_Gamma phi q_n [typical conductance and heat flux combination, for example]
        //
        // This method assumes
        //     x      = x0 + dt x0_dot + (1/2-beta) dt^2 x0_ddot + beta dt^2 x_ddot
        // or, x_ddot = (x-x0)/beta/dt^2 - 1/beta/dt x0_dot - (1/2-beta)/beta x0_ddot
        //
        //     x_dot  = x0_dot + (1-gamma) dt x0_ddot + gamma dt x_ddot
        // or, x_dot  = x0_dot + (1-gamma) dt x0_ddot +
        //              gamma dt [(x-x0)/beta/dt^2 - 1/beta/dt x0_dot - (1/2-beta)/beta x0_ddot]
        //            = x0_dot + (1-gamma) dt x0_ddot +
        //              gamma/beta/dt (x-x0) - gamma/beta x0_dot - gamma (1/2-beta)/beta dt x0_ddot
        //            = gamma/beta/dt (x-x0) + (1 - gamma/beta) x0_dot + (1 - gamma/2/beta) dt x0_ddot
        //
        // Both f_m and f_x can be functions of x_dot and x. Then, the
        // Jacobian is
        // dr/dx = df_m/dx + df_x/dx +
        //         df_m/dx_dot dx_dot/dx + df_x/dx_dot dx_dot/dx
        //       = (df_m/dx + df_x/dx) +
        //         (df_m/dx_dot + df_x/dx_dot) (gamma/beta/dt) +
        //         df_m/dx_ddot (1/beta/dt/dt) +
        //
        
        
        // system residual
        vec  = (f_m + f_x);
        
        // system Jacobian
        if (if_jac)
            mat = ((1./beta/dt/dt) *  f_m_jac_xddot +
                   (gamma/beta/dt)* (f_m_jac_xdot+f_x_jac_xdot) +
                   (f_m_jac + f_x_jac));
    }
}



void
MAST::SecondOrderNewmarkTransientSolver::
elem_linearized_jacobian_solution_product(RealVectorX& vec) {
    
    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    
    // perform the element assembly
    _assembly_ops->linearized_jacobian_solution_product(vec);
}




void
MAST::SecondOrderNewmarkTransientSolver::
elem_sensitivity_calculations(const MAST::FunctionBase& f,
                              RealVectorX& vec) {
    
    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    unsigned int n_dofs = (unsigned int)vec.size();
    
    RealVectorX
    f_x     = RealVectorX::Zero(n_dofs),
    f_m     = RealVectorX::Zero(n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_sensitivity_calculations(f,
                                                 f_m,           // mass vector
                                                 f_x);          // forcing vector
    
    // system residual
    vec  = (f_m + f_x);
}


void
MAST::SecondOrderNewmarkTransientSolver::
elem_sensitivity_contribution_previous_timestep(const std::vector<RealVectorX>& prev_sols,
                                                RealVectorX& vec) {
    
    // nothing to be done for highest derivative term
    if (_if_highest_derivative_solution) return;

    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    libmesh_assert_equal_to(prev_sols.size(), 3);
    
    unsigned int n_dofs = (unsigned int)vec.size();
    
    const RealVectorX
    &sol          = prev_sols[0],
    &vel          = prev_sols[1],
    &accel        = prev_sols[2];
    RealVectorX
    dummy_vec     = RealVectorX::Zero(n_dofs);
    
    RealMatrixX
    f_m_jac_xdot  = RealMatrixX::Zero(n_dofs, n_dofs),
    dummy_mat     = RealMatrixX::Zero(n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_calculations(true,
                                     dummy_vec,           // mass vector
                                     dummy_vec,           // forcing vector
                                     f_m_jac_xdot,        // Jac of mass wrt x_dot
                                     dummy_mat,           // Jac of mass wrt x
                                     dummy_mat);          // Jac of forcing vector wrt x
    
    libmesh_assert(false); // this expression is in error and needs to be set
    vec  -= f_m_jac_xdot * ( (1./beta/dt)*sol + (1.-beta)/beta * vel);
}


void
MAST::SecondOrderNewmarkTransientSolver::
elem_shape_sensitivity_calculations(const MAST::FunctionBase& f,
                                    RealVectorX& vec) {
    
    libmesh_assert(false); // to be implemented
}



void
MAST::SecondOrderNewmarkTransientSolver::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       RealVectorX& vec) {
    libmesh_assert(false); // to be implemented
}



void
MAST::SecondOrderNewmarkTransientSolver::
elem_topology_sensitivity_calculations(const MAST::FunctionBase& f,
                                       const MAST::FieldFunction<RealVectorX>& vel,
                                       RealVectorX& vec) {
    libmesh_assert(false); // to be implemented
}



