
// MAST includes
#include "solver/first_order_newmark_transient_solver.h"
#include "base/transient_assembly.h"

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
    
    // ask the Newton solver to solve for the system solution
    _system->solve();
    
}



void
MAST::FirstOrderNewmarkTransientSolver::advance_time_step() {
    
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
    
    // get vectors from previous iterations
    const libMesh::NumericVector<Real>
    & curr_solution = solution(0),
    & old_solution  = solution(1),
    & old_velocity  = velocity(1);
    
    RealVectorX x_dot, vec2, old_vel, old_sol, sol;
    RealMatrixX mass, jac;
    
    unsigned int n_dofs = (unsigned int)dof_indices.size();
    
    x_dot.resize(n_dofs);
    vec2.resize(n_dofs);
    old_vel.resize(n_dofs);
    old_sol.resize(n_dofs);
    sol.resize(n_dofs);
    
    for (unsigned int i=0; i<n_dofs; i++) {
        old_vel(i) = old_velocity(dof_indices[i]);
        old_sol(i) = old_solution(dof_indices[i]);
        sol(i)     = curr_solution(dof_indices[i]);
    }
    
    // perform the element assembly
    _assembly->_elem_calculations(elem, if_jac, x_dot, mass, jac);
    
    vec2.resize(x_dot.size());
    
    vec = x_dot;
    vec *= -beta*dt;
    
    vec2 = mass*old_vel;
    vec -= dt*(1.-beta)*vec2;
    
    vec2 = mass*old_sol;
    vec -= vec2;
    
    vec2 = mass*sol;
    vec += vec2;
    
    // set the system residual
    if (if_jac)
        mat = mass - beta*dt*jac;
}



void
MAST::FirstOrderNewmarkTransientSolver::
_elem_sensitivity_calculations(MAST::ElementBase& elem,
                               const std::vector<libMesh::dof_id_type>& dof_indices,
                               RealVectorX& vec) {
    
}
