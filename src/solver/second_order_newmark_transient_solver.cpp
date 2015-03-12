
// MAST includes
#include "solver/second_order_newmark_transient_solver.h"
#include "base/transient_assembly.h"
#include "base/elem_base.h"

// libMesh includes
#include "libmesh/numeric_vector.h"


MAST::SecondOrderNewmarkTransientSolver::SecondOrderNewmarkTransientSolver():
MAST::TransientSolverBase(),
beta(0.5),
gamma(0.25)
{ }


MAST::SecondOrderNewmarkTransientSolver::~SecondOrderNewmarkTransientSolver()
{ }



void
MAST::SecondOrderNewmarkTransientSolver::solve() {
    
    // make sure that the system has been specified
    libmesh_assert_msg(_system, "System pointer is NULL.");
    
    // ask the Newton solver to solve for the system solution
    _system->solve();
    
}



void
MAST::SecondOrderNewmarkTransientSolver::advance_time_step() {
    
    // use the parent class' advance time scheme
    MAST::TransientSolverBase::advance_time_step();
}



void
MAST::SecondOrderNewmarkTransientSolver::
_set_element_data(std::vector<libMesh::dof_id_type>& dof_indices,
                  MAST::ElementBase &elem){
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol         (n_dofs),
    vel         (n_dofs),
    accel       (n_dofs),
    prev_sol    (n_dofs),
    prev_vel    (n_dofs),
    prev_accel  (n_dofs);
    
    
    // get the references to current and previous sol and velocity
    const libMesh::NumericVector<Real>
    &sol_global     =   solution(0),
    &old_sol_global =   solution(1),
    &old_vel_global =   velocity(1),
    &old_acc_global =   acceleration(1);
    
    for (unsigned int i=0; i<n_dofs; i++) {
        sol(i)          = sol_global(dof_indices[i]);
        prev_sol(i)     = old_sol_global(dof_indices[i]);
        prev_vel(i)     = old_vel_global(dof_indices[i]);
        prev_accel(i)   = old_acc_global(dof_indices[i]);
    }
    
    
    
    accel = (sol-prev_sol - dt*prev_vel - (0.5-beta)*dt*dt*prev_accel)/(beta*dt*dt);
    vel   = (prev_vel + (1.-gamma)*dt*prev_accel + gamma*dt*accel);
    
    elem.set_solution(sol);
    elem.set_velocity(vel);
    elem.set_acceleration(accel);
}




void
MAST::SecondOrderNewmarkTransientSolver::
_update_velocity(libMesh::NumericVector<Real>& vec) {
    
    const libMesh::NumericVector<Real>
    &sol      = _system->get_vector("transient_solution_0"),
    &prev_sol = _system->get_vector("transient_solution_1"),
    &prev_vel = _system->get_vector("transient_velocity_1"),
    &prev_acc = _system->get_vector("transient_acceleration_1");

    // first calculate the acceleration
    vec = sol;
    vec.add(-1., prev_sol);
    vec.add( dt, prev_vel);
    vec.add(-(.5-beta)*dt*dt, prev_acc);
    vec.scale(gamma/(beta*dt));
    // now add the rest of the terms
    vec.add(1., prev_vel);
    vec.add((1.-gamma)*dt, prev_acc);

    vec.close();
}




void
MAST::SecondOrderNewmarkTransientSolver::
_update_acceleration(libMesh::NumericVector<Real>& vec) {
    
    const libMesh::NumericVector<Real>
    &sol      = _system->get_vector("transient_solution_0"),
    &prev_sol = _system->get_vector("transient_solution_1"),
    &prev_vel = _system->get_vector("transient_velocity_1"),
    &prev_acc = _system->get_vector("transient_acceleration_1");
    
    vec = sol;
    vec.add(-1., prev_sol);
    vec.add( dt, prev_vel);
    vec.add(-(.5-beta)*dt*dt, prev_acc);
    vec.scale(1./(beta*dt*dt));
    
    vec.close();
}





void
MAST::SecondOrderNewmarkTransientSolver::
_elem_calculations(MAST::ElementBase& elem,
                   const std::vector<libMesh::dof_id_type>& dof_indices,
                   bool if_jac,
                   RealVectorX& vec,
                   RealMatrixX& mat) {
    // make sure that the assembly object is provides
    libmesh_assert(_assembly);
    unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    RealVectorX
    f_x    (n_dofs),
    f_m    (n_dofs);
    
    RealMatrixX
    f_m_jac_xddot   (n_dofs, n_dofs),
    f_m_jac_xdot    (n_dofs, n_dofs),
    f_m_jac         (n_dofs, n_dofs),
    f_x_jac_xdot    (n_dofs, n_dofs),
    f_x_jac         (n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly->_elem_calculations(elem,
                                  if_jac,
                                  f_m,           // mass vector
                                  f_x,           // forcing vector
                                  f_m_jac_xddot, // Jac of mass wrt x_dotdot
                                  f_m_jac_xdot,  // Jac of mass wrt x_dot
                                  f_m_jac,       // Jac of mass wrt x
                                  f_x_jac_xdot,  // Jac of forcing vector wrt x_dot
                                  f_x_jac);      // Jac of forcing vector wrt x
    
    //
    //  The residual here is modeled as
    // r = beta*dt(f_m + f_x )= 0
    // where, (for example)
    // f_m = int_Omega phi u_dot   [typical mass vector in conduction, for example]
    // f_x = int_Omega phi_i u_i - int_Gamma phi q_n [typical conductance and heat flux combination, for example]
    //
    // This method assumes
    //     x     = x0 + (1-beta) dt x0_dot + beta dt x_dot
    // or, x_dot = (x-x0)/beta/dt - (1-beta)/beta x0_dot
    // Note that the residual expression multiplies the expression by beta*dt
    // for convenience in the following expressions
    //
    // Both f_m and f_x can be functions of x_dot and x. Then, the
    // Jacobian is
    // dr/dx = beta*dt [df_m/dx + df_x/dx +
    //                 df_m/dx_dot dx_dot/dx + df_x/dx_dot dx_dot/dx]
    //       = beta*dt [(df_m/dx + df_x/dx) +
    //                  (df_m/dx_dot + df_x/dx_dot) (1/beta/dt)]
    //       = beta*dt (df_m/dx + df_x/dx) +
    //                 (df_m/dx_dot + df_x/dx_dot)
    // Note that this form of equations makes it a good candidate for
    // use as implicit solver, ie, for a nonzero beta.
    //
    
    
    // system residual
    vec  = beta*dt*dt*(f_m + f_x);
    
    // system Jacobian
    if (if_jac)
        mat = f_m_jac_xddot +
        gamma*dt*    (f_m_jac_xdot+f_x_jac_xdot) +
        beta*dt*dt*  (f_m_jac + f_x_jac);
}



void
MAST::SecondOrderNewmarkTransientSolver::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               const std::vector<libMesh::dof_id_type>& dof_indices,
                               RealVectorX& vec) {
    
}
