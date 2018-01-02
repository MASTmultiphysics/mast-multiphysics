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
#include "solver/first_order_newmark_transient_solver.h"
#include "base/transient_assembly_elem_operations.h"
#include "base/elem_base.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/numeric_vector.h"


MAST::FirstOrderNewmarkTransientSolver::FirstOrderNewmarkTransientSolver():
MAST::TransientSolverBase(),
beta(1.)
{ }


MAST::FirstOrderNewmarkTransientSolver::~FirstOrderNewmarkTransientSolver()
{ }



void
MAST::FirstOrderNewmarkTransientSolver::solve() {
    
    // make sure that the system has been specified
    libmesh_assert_msg(_system, "System pointer is nullptr.");
    
    // ask the Newton solver to solve for the system solution
    _system->solve();
    
}



void
MAST::FirstOrderNewmarkTransientSolver::advance_time_step() {
    
    // use the parent class' advance time scheme
    MAST::TransientSolverBase::advance_time_step();
}



void
MAST::FirstOrderNewmarkTransientSolver::
set_element_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                 const std::vector<libMesh::NumericVector<Real>*>& sols,
                 MAST::ElementBase &elem){
    
    libmesh_assert_equal_to(sols.size(), 2);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol          = RealVectorX::Zero(n_dofs),
    vel          = RealVectorX::Zero(n_dofs);


    const libMesh::NumericVector<Real>
    &sol_global = *sols[0],
    &vel_global = *sols[1];
    
    // get the references to current and previous sol and velocity
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
    }
    
    elem.set_solution(sol);
    elem.set_velocity(vel);
}



void
MAST::FirstOrderNewmarkTransientSolver::
set_element_perturbed_data(const std::vector<libMesh::dof_id_type>& dof_indices,
                           const std::vector<libMesh::NumericVector<Real>*>& sols,
                           MAST::ElementBase &elem){
    
    libmesh_assert_equal_to(sols.size(), 2);
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol          = RealVectorX::Zero(n_dofs),
    vel          = RealVectorX::Zero(n_dofs);
    
    
    const libMesh::NumericVector<Real>
    &sol_global = *sols[0],
    &vel_global = *sols[1];
    
    // get the references to current and previous sol and velocity
    
    for (unsigned int i=0; i<n_dofs; i++) {
        
        sol(i)          = sol_global(dof_indices[i]);
        vel(i)          = vel_global(dof_indices[i]);
    }
    
    elem.set_perturbed_solution(sol);
    elem.set_perturbed_velocity(vel);
}





void
MAST::FirstOrderNewmarkTransientSolver::
update_velocity(libMesh::NumericVector<Real>&       vec,
                 const libMesh::NumericVector<Real>& sol) {
    
    const libMesh::NumericVector<Real>
    &prev_sol = this->solution(1),
    &prev_vel = this->velocity(1);
    
    vec.zero();
    vec.add( 1.,      sol);
    vec.add(-1., prev_sol);
    vec.scale(1./beta/dt);
    vec.close();
    vec.add(-(1.-beta)/beta, prev_vel);
    
    vec.close();
}



void
MAST::FirstOrderNewmarkTransientSolver::
update_delta_velocity(libMesh::NumericVector<Real>&       vec,
                       const libMesh::NumericVector<Real>& sol) {
    
    vec.zero();
    vec.add( 1./beta/dt,      sol);
    vec.close();
}





void
MAST::FirstOrderNewmarkTransientSolver::
elem_calculations(MAST::ElementBase& elem,
                  bool if_jac,
                  RealVectorX& vec,
                  RealMatrixX& mat) {
    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    unsigned int n_dofs = (unsigned int)vec.size();

    RealVectorX
    f_x     = RealVectorX::Zero(n_dofs),
    f_m     = RealVectorX::Zero(n_dofs);
    
    RealMatrixX
    f_m_jac_xdot  = RealMatrixX::Zero(n_dofs, n_dofs),
    f_m_jac       = RealMatrixX::Zero(n_dofs, n_dofs),
    f_x_jac       = RealMatrixX::Zero(n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly_ops->elem_calculations(elem,
                                     if_jac,
                                     f_m,           // mass vector
                                     f_x,           // forcing vector
                                     f_m_jac_xdot,  // Jac of mass wrt x_dot
                                     f_m_jac,       // Jac of mass wrt x
                                     f_x_jac);      // Jac of forcing vector wrt x

    if (_if_highest_derivative_solution) {
        
        //
        //  The residual here is modeled as
        // r(x, xdot) = f_m(x, xdot) + f_x(x, xdot)= 0
        //
        //  then, given x = x0, the residual can be used to evaluate xdot0 at
        //  the same time step as
        //
        //  r(x0, xdot) + dr/dxdot dxdot = 0
        //
        
        // system residual
        vec  = (f_m + f_x);
        
        // system Jacobian
        if (if_jac)
            mat = f_m_jac_xdot;
    }
    else {
        //
        //  The residual here is modeled as
        // r = (f_m + f_x )= 0
        // where, (for example)
        // f_m = int_Omega phi u_dot   [typical mass vector in conduction, for example]
        // f_x = int_Omega phi_i u_i - int_Gamma phi q_n [typical conductance and heat flux combination, for example]
        //
        // This method assumes
        //     x     = x0 + (1-beta) dt x0_dot + beta dt x_dot
        // or, x_dot = (x-x0)/beta/dt - (1-beta)/beta x0_dot
        //
        // Both f_m and f_x can be functions of x_dot and x. Then, the
        // Jacobian is
        // dr/dx =[df_m/dx + df_x/dx +
        //         df_m/dx_dot dx_dot/dx + df_x/dx_dot dx_dot/dx]
        //       = [(df_m/dx + df_x/dx) +
        //          (df_m/dx_dot + df_x/dx_dot) (1/beta/dt)]
        //       = (df_m/dx + df_x/dx) +
        //         1/(beta*dt)(df_m/dx_dot + df_x/dx_dot)
        // Note that this form of equations makes it a good candidate for
        // use as implicit solver, ie, for a nonzero beta.
        //
        
        
        // system residual
        vec  = (f_m + f_x);
        
        // system Jacobian
        if (if_jac)
            mat = (1./beta/dt)*f_m_jac_xdot + (f_m_jac + f_x_jac);
    }
}



void
MAST::FirstOrderNewmarkTransientSolver::
elem_linearized_jacobian_solution_product(MAST::ElementBase& elem,
                                          RealVectorX& vec) {

    // make sure that the assembly object is provided
    libmesh_assert(_assembly_ops);
    
    // perform the element assembly
    _assembly_ops->linearized_jacobian_solution_product(elem, vec);
}




void
MAST::FirstOrderNewmarkTransientSolver::
elem_sensitivity_calculations(MAST::ElementBase& elem,
                              RealVectorX& vec) {

    // to be implemented
    libmesh_error();
}
