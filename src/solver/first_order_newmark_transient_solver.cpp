
// MAST includes
#include "solver/first_order_newmark_transient_solver.h"
#include "base/transient_assembly.h"
#include "base/elem_base.h"

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
    libmesh_assert_msg(_system, "System pointer is NULL.");
    
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
_set_element_data(std::vector<libMesh::dof_id_type>& dof_indices,
                  MAST::ElementBase &elem){
    
    const unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    // get the current state and velocity estimates
    // also get the current discrete velocity replacement
    RealVectorX
    sol         (n_dofs),
    vel         (n_dofs),
    prev_sol    (n_dofs),
    prev_vel    (n_dofs);

    
    // get the references to current and previous sol and velocity
    const libMesh::NumericVector<Real>
    &sol_global     =   solution(0),
    &old_sol_global =   solution(1),
    &old_vel_global =   velocity(1);
    
    for (unsigned int i=0; i<n_dofs; i++) {
        sol(i)          = sol_global(dof_indices[i]);
        prev_sol(i)     = old_sol_global(dof_indices[i]);
        prev_vel(i)     = old_vel_global(dof_indices[i]);
    }
    
    
    
    vel = (sol-prev_sol - (1.-beta)*dt*prev_vel)/(beta*dt);
    
    elem.set_solution(sol);
    elem.set_velocity(vel);
}




void
MAST::FirstOrderNewmarkTransientSolver::
_update_velocity(libMesh::NumericVector<Real>& vec) {

    const libMesh::NumericVector<Real>
    &sol      = _system->get_vector("transient_solution_0"),
    &prev_sol = _system->get_vector("transient_solution_1"),
    &prev_vel = _system->get_vector("transient_velocity_1");
    
    vec = sol;
    vec.add(-1., prev_sol);
    vec.scale(1./beta/dt);
    vec.add(-(1.-beta)/beta, prev_vel);
    
    vec.close();
}





void
MAST::FirstOrderNewmarkTransientSolver::
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
    f_m_jac_xdot(n_dofs, n_dofs),
    f_m_jac     (n_dofs, n_dofs),
    f_x_jac     (n_dofs, n_dofs);
    
    // perform the element assembly
    _assembly->_elem_calculations(elem,
                                  if_jac,
                                  f_m,           // mass vector
                                  f_x,           // forcing vector
                                  f_m_jac_xdot, // Jac of mass wrt x_dot
                                  f_m_jac,       // Jac of mass wrt x
                                  f_x_jac);      // Jac of forcing vector wrt x
    
    // system residual
    vec  = beta*dt*(f_m + f_x);
    
    // system Jacobian
    if (if_jac)
        mat = f_m_jac_xdot + beta*dt*(f_m_jac + f_x_jac);
}



void
MAST::FirstOrderNewmarkTransientSolver::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               const std::vector<libMesh::dof_id_type>& dof_indices,
                               RealVectorX& vec) {
    
}
