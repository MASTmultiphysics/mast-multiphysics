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
#include "solver/transient_solver_base.h"
#include "base/transient_assembly_elem_operations.h"
#include "base/assembly_base.h"
#include "base/transient_assembly.h"
#include "base/nonlinear_system.h"
#include "base/system_initialization.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"



MAST::TransientSolverBase::TransientSolverBase(unsigned int o,
                                               unsigned int n):
MAST::NonlinearImplicitAssemblyElemOperations(),
dt                              (0.),
_first_step                     (true),
_first_sensitivity_step         (true),
_ode_order                      (o),
_n_iters_to_store               (n),
_assembly_ops                   (nullptr),
_if_highest_derivative_solution (false) {

}




MAST::TransientSolverBase::~TransientSolverBase() {
    
    this->clear_elem_operation_object();
}



void
MAST::TransientSolverBase::set_assembly(MAST::AssemblyBase &assembly) {
    
    libmesh_assert(_assembly_ops);
    MAST::NonlinearImplicitAssemblyElemOperations::set_assembly(assembly);
    _assembly_ops->set_assembly(assembly);
}


void
MAST::TransientSolverBase::clear_assembly() {
    
    MAST::NonlinearImplicitAssemblyElemOperations::clear_assembly();
    if (_assembly_ops)
        _assembly_ops->clear_assembly();
}


MAST::TransientAssemblyElemOperations&
MAST::TransientSolverBase::get_elem_operation_object() {
    
    libmesh_assert(_assembly_ops);
    
    return *_assembly_ops;
}


void
MAST::TransientSolverBase::
set_elem_operation_object(MAST::TransientAssemblyElemOperations &assembly_ops) {
    
    libmesh_assert(_system);
    libmesh_assert(_discipline);

    MAST::NonlinearSystem
    &sys = _system->system();
    
    // make sure that the previous assembly association has been cleared
    this->clear_elem_operation_object();
    
    _assembly_ops = &assembly_ops;
    
    // now, add the vectors
    std::string nm;
    for (unsigned int i=0; i<_n_iters_to_store; i++) {
        std::ostringstream iter;
        iter << i;
        
        // add the solution only for previous iterations. The current
        // solution is obtained from system->solution
        if (i>0) {
            
            nm = "transient_solution_";
            nm += iter.str();
            sys.add_vector(nm, true, libMesh::GHOSTED);

            nm = "transient_solution_sensitivity_";
            nm += iter.str();
            sys.add_vector(nm, true, libMesh::GHOSTED);
        }

        // add the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        sys.add_vector(nm, true, libMesh::GHOSTED);
        nm = "transient_velocity_sensitivity_";
        nm += iter.str();
        sys.add_vector(nm, true, libMesh::GHOSTED);

        if (_ode_order > 1) {
            // add the acceleration
            nm = "transient_acceleration_";
            nm += iter.str();
            sys.add_vector(nm, true, libMesh::GHOSTED);
            nm = "transient_acceleration_sensitivity_";
            nm += iter.str();
            sys.add_vector(nm, true, libMesh::GHOSTED);
        }
    }
    
    _first_step = true;
}



void
MAST::TransientSolverBase::clear_elem_operation_object() {
    
    // if no system has been set so far, nothing needs to be done
    if (!_system)
        return;
    
    MAST::NonlinearSystem
    &sys = _system->system();
    
    // clear the transient solutions stored in system for solution
    // number of time steps to store
    std::string nm;
    for (unsigned int i=0; i<_n_iters_to_store; i++) {
        std::ostringstream iter;
        iter << i;
        
        // remove the solution
        nm = "transient_solution_";
        nm += iter.str();
        if (sys.have_vector(nm)) sys.remove_vector(nm);
        nm = "transient_solution_sensitivity_";
        nm += iter.str();
        if (sys.have_vector(nm)) sys.remove_vector(nm);

        // remove the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        if (sys.have_vector(nm)) sys.remove_vector(nm);
        nm = "transient_velocity_sensitivity_";
        nm += iter.str();
        if (sys.have_vector(nm)) sys.remove_vector(nm);

        if (_ode_order > 1) {
            // remove the acceleration
            nm = "transient_acceleration_";
            nm += iter.str();
            if (sys.have_vector(nm)) sys.remove_vector(nm);
            nm = "transient_acceleration_sensitivity_";
            nm += iter.str();
            if (sys.have_vector(nm)) sys.remove_vector(nm);

        }
    }
    
    _assembly_ops   = nullptr;
    _first_step     = true;
}




libMesh::NumericVector<Real>&
MAST::TransientSolverBase::solution(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    if (prev_iter) {
        // get references to the solution
        std::string
        nm = "transient_solution_" + oss.str();

        return _system->system().get_vector(nm);
    }
    else
        return *_system->system().current_local_solution;
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::solution_sensitivity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;

    if (prev_iter) {
        // get references to the solution
        std::string
        nm = "transient_solution_sensitivity_" + oss.str();
        
        return _system->system().get_vector(nm);
    }
    else
        return _system->system().add_sensitivity_solution();
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::velocity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_velocity_" + oss.str();
    
    return _system->system().get_vector(nm);
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::velocity_sensitivity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_velocity_sensitivity_" + oss.str();
    
    return _system->system().get_vector(nm);
}




libMesh::NumericVector<Real>&
MAST::TransientSolverBase::acceleration(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_acceleration_" + oss.str();
    
    return _system->system().get_vector(nm);
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::acceleration_sensitivity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store);
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_acceleration_sensitivity_" + oss.str();
    
    return _system->system().get_vector(nm);
}



void
MAST::TransientSolverBase::solve(MAST::AssemblyBase& assembly) {
    
    // make sure that the system has been specified
    libmesh_assert_msg(_system, "System pointer is nullptr.");
    
    // ask the Newton solver to solve for the system solution
    _system->system().solve(*this, assembly);
    
}



void
MAST::TransientSolverBase::sensitivity_solve(MAST::AssemblyBase& assembly,
                                             const MAST::FunctionBase& f) {
    
    // make sure that the system has been specified
    libmesh_assert_msg(_system, "System pointer is nullptr.");
    
    // ask the Newton solver to solve for the system solution
    _system->system().sensitivity_solve(*this, assembly, f);
}



void
MAST::TransientSolverBase::
solve_highest_derivative_and_advance_time_step(MAST::AssemblyBase& assembly) {
    
    libmesh_assert(_first_step);
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(!_assembly);
    
    // tell the solver that the current solution being obtained is for the
    // highest time derivative
    _if_highest_derivative_solution = true;
    
    MAST::NonlinearSystem
    &sys = _system->system();
    
    // Build the residual and Jacobian
    assembly.set_elem_operation_object(*this);
    assembly.residual_and_jacobian(*sys.solution, sys.rhs, sys.matrix, sys);
    assembly.clear_elem_operation_object();
    
    // reset the solution flag
    _if_highest_derivative_solution = false;
    
    std::pair<unsigned int, Real>
    solver_params = sys.get_linear_solve_parameters();
    
    // Solve the linear system.
    libMesh::SparseMatrix<Real> *
    pc = sys.request_matrix("Preconditioner");
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    dvec(velocity().zero_clone().release());

    sys.linear_solver->solve (*sys.matrix, pc,
                              *dvec,
                              *sys.rhs,
                              solver_params.second,
                              solver_params.first);

    libMesh::NumericVector<Real> *vec = nullptr;
    
    switch (_ode_order) {
            
        case 1:
            vec = &velocity();
            break;
            
        case 2:
            vec = &acceleration();
            break;
            
        default:
            // higher than 2 derivative not implemented yet.
            libmesh_error();
    }
    
    vec->add(-1., *dvec);
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    sys.get_dof_map().enforce_constraints_exactly(sys, vec, /* homogeneous = */ true);
#endif
    
    sys.update();
    
    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store-1; i>0; i--) {
        this->solution(i).zero();
        this->solution(i).add(1., this->solution(i-1));
        this->solution(i).close();
        
        this->velocity(i).zero();
        this->velocity(i).add(1., this->velocity(i-1));
        this->velocity(i).close();
        
        if (_ode_order > 1) {
            
            this->acceleration(i).zero();
            this->acceleration(i).add(1., this->acceleration(i-1));
            this->acceleration(i).close();
        }
    }
    
    // finally, update the system time
    sys.time          += dt;
    _first_step        = false;
}




void
MAST::TransientSolverBase::
solve_highest_derivative_and_advance_time_step_with_sensitivity(MAST::AssemblyBase& assembly,
                                                                const MAST::FunctionBase& f) {
    
    libmesh_assert(_first_sensitivity_step);
    libmesh_assert(_system);
    libmesh_assert(_discipline);
    libmesh_assert(!_assembly);
    
    // tell the solver that the current solution being obtained is for the
    // highest time derivative
    _if_highest_derivative_solution = true;

    MAST::NonlinearSystem
    &sys = _system->system();
    
    // Build the Jacobian
    assembly.set_elem_operation_object(*this);
    assembly.residual_and_jacobian(*sys.solution, nullptr, sys.matrix, sys);
    // this should assembl the sensitivity of the residual
    assembly.sensitivity_assemble(f, *sys.rhs);
    assembly.clear_elem_operation_object();

    // reset the solution flag
    _if_highest_derivative_solution = false;
    
    std::pair<unsigned int, Real>
    solver_params = sys.get_linear_solve_parameters();
    
    // Solve the linear system.
    libMesh::SparseMatrix<Real> *
    pc = sys.request_matrix("Preconditioner");
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    dvec(velocity().zero_clone().release());
    
    sys.linear_solver->solve (*sys.matrix, pc,
                              *dvec,
                              *sys.rhs,
                              solver_params.second,
                              solver_params.first);
    
    libMesh::NumericVector<Real> *vec = nullptr;
    
    switch (_ode_order) {
            
        case 1:
            vec = &velocity_sensitivity();
            break;
            
        case 2:
            vec = &acceleration_sensitivity();
            break;
            
        default:
            // higher than 2 derivative not implemented yet.
            libmesh_error();
    }
    
    vec->add(-1., *dvec);
    
    // The linear solver may not have fit our constraints exactly
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    sys.get_dof_map().enforce_constraints_exactly(sys, vec, /* homogeneous = */ true);
#endif
    
    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store-1; i>0; i--) {
        
        this->solution_sensitivity(i).zero();
        this->solution_sensitivity(i).add(1., this->solution_sensitivity(i-1));
        this->solution_sensitivity(i).close();
        
        this->velocity_sensitivity(i).zero();
        this->velocity_sensitivity(i).add(1., this->velocity_sensitivity(i-1));
        this->velocity_sensitivity(i).close();
        
        if (_ode_order > 1) {
            
            this->acceleration_sensitivity(i).zero();
            this->acceleration_sensitivity(i).add(1., this->acceleration_sensitivity(i-1));
            this->acceleration_sensitivity(i).close();
        }
    }
    
    // finally, update the system time
    sys.time          += dt;
    _first_sensitivity_step        = false;
}



void
MAST::TransientSolverBase::
build_local_quantities(const libMesh::NumericVector<Real>& current_sol,
                       std::vector<libMesh::NumericVector<Real>*>& sol) {

    // make sure that the system has been specified
    libmesh_assert(_system);

    MAST::NonlinearSystem
    &sys = _system->system();
    
    // make sure there are no solutions in sol
    libmesh_assert(!sol.size());

    // resize the vector to store the quantities
    sol.resize(_ode_order+1);
    
    const std::vector<libMesh::dof_id_type>&
    send_list = sys.get_dof_map().get_send_list();
    

    for ( unsigned int i=0; i<=_ode_order; i++) {
        
        sol[i] = libMesh::NumericVector<Real>::build(sys.comm()).release();
        sol[i]->init(sys.n_dofs(),
                     sys.n_local_dofs(),
                     send_list,
                     false,
                     libMesh::GHOSTED);

        switch (i) {
                
            case 0: {
                
                if (!_if_highest_derivative_solution)
                    // copy the solution to this
                    current_sol.localize(*sol[i], send_list);
                else
                    solution().localize(*sol[i], send_list);
            }
                break;

            case 1: {
                
                // update the current local velocity vector
                libMesh::NumericVector<Real>& vel = this->velocity();

                if (!_if_highest_derivative_solution)
                    // calculate the velocity and localize it
                    update_velocity(vel, current_sol);
                
                vel.localize(*sol[i], send_list);
            }
                break;

            case 2: {
                
                // update the current local acceleration vector
                libMesh::NumericVector<Real>& acc = this->acceleration();
                
                if (!_if_highest_derivative_solution)
                    // calculate the acceleration and localize it
                    update_acceleration(acc, current_sol);
                
                acc.localize(*sol[i], send_list);
            }
                break;

            default:
                // should not get here
                libmesh_error();
                break;
        }
    }
}



void
MAST::TransientSolverBase::
build_sensitivity_local_quantities(unsigned int prev_iter,
                                   std::vector<libMesh::NumericVector<Real>*>& sol) {
    
    // make sure that the system has been specified
    libmesh_assert(_system);
    libmesh_assert_less_equal(prev_iter, _n_iters_to_store);
    
    MAST::NonlinearSystem
    &sys = _system->system();
    
    // make sure there are no solutions in sol
    libmesh_assert(!sol.size());
    
    // resize the vector to store the quantities
    sol.resize(_ode_order+1);
    
    const std::vector<libMesh::dof_id_type>&
    send_list = sys.get_dof_map().get_send_list();
    
    
    for ( unsigned int i=0; i<=_ode_order; i++) {
        
        sol[i] = libMesh::NumericVector<Real>::build(sys.comm()).release();
        sol[i]->init(sys.n_dofs(),
                     sys.n_local_dofs(),
                     send_list,
                     false,
                     libMesh::GHOSTED);
        
        switch (i) {
                
            case 0:
                solution_sensitivity(prev_iter).localize(*sol[i], send_list);
                break;
                
            case 1:
                velocity_sensitivity(prev_iter).localize(*sol[i], send_list);
                break;
                
            case 2:
                acceleration_sensitivity(prev_iter).localize(*sol[i], send_list);
                break;
                
            default:
                // should not get here
                libmesh_error();
                break;
        }
    }
}





void
MAST::TransientSolverBase::
build_perturbed_local_quantities(const libMesh::NumericVector<Real>& current_dsol,
                                 std::vector<libMesh::NumericVector<Real>*>& sol) {
    
    // make sure that the system has been specified
    libmesh_assert(_system);

    MAST::NonlinearSystem
    &sys = _system->system();
    
    // make sure there are no solutions in sol
    libmesh_assert(!sol.size());
    
    // resize the vector to store the quantities
    sol.resize(_ode_order+1);
    
    const std::vector<libMesh::dof_id_type>&
    send_list = sys.get_dof_map().get_send_list();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    tmp(this->solution().zero_clone().release());
    
    for ( unsigned int i=0; i<=_ode_order; i++) {
        
        sol[i] = libMesh::NumericVector<Real>::build(sys.comm()).release();
        sol[i]->init(sys.n_dofs(),
                     sys.n_local_dofs(),
                     send_list,
                     false,
                     libMesh::GHOSTED);
        
        switch (i) {
                
            case 0: {
                
                current_dsol.localize(*sol[i], send_list);
                if (_if_highest_derivative_solution)
                    // only the highest derivative is perturbed since
                    // all lower derivative quantities are known
                    sol[i]->zero();
            }
                break;
                
            case 1: {
                
                if (_ode_order == 1 && _if_highest_derivative_solution)
                    // the provided solution is the current velocity increment
                    current_dsol.localize(*sol[i], send_list);
                else if (_if_highest_derivative_solution)
                    sol[i]->zero();
                else {
                    update_delta_velocity(*tmp, current_dsol);
                    tmp->localize(*sol[i], send_list);
                }
            }
                break;
                
            case 2: {
                
                if (_ode_order == 2 && _if_highest_derivative_solution)
                    // the provided solution is the current velocity increment
                    current_dsol.localize(*sol[i], send_list);
                else if (_if_highest_derivative_solution)
                    sol[i]->zero();
                else {
                    update_delta_acceleration(*tmp, current_dsol);
                    tmp->localize(*sol[i], send_list);
                }
            }
                break;
                
            default:
                // should not get here
                libmesh_error();
                break;
        }
    }
}




void
MAST::TransientSolverBase::advance_time_step() {

    libmesh_assert(_system);
    libmesh_assert(_discipline);
    
    MAST::NonlinearSystem
    &sys = _system->system();
    
    // first ask the solver to update the velocity and acceleration vector
    sys.update();
    update_velocity(this->velocity(), *sys.current_local_solution);

    if (_ode_order > 1)
        update_acceleration(this->acceleration(), *sys.current_local_solution);

    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store-1; i>0; i--) {
        this->solution(i).zero();
        this->solution(i).add(1., this->solution(i-1));
        this->solution(i).close();
        
        this->velocity(i).zero();
        this->velocity(i).add(1., this->velocity(i-1));
        this->velocity(i).close();
        
        if (_ode_order > 1) {
            
            this->acceleration(i).zero();
            this->acceleration(i).add(1., this->acceleration(i-1));
            this->acceleration(i).close();
        }
    }

    // finally, update the system time
    sys.time          += dt;
    _first_step        = false;
}



void
MAST::TransientSolverBase::advance_time_step_with_sensitivity() {
    
    MAST::NonlinearSystem
    &sys = _system->system();

    // first ask the solver to update the velocity and acceleration vector
    update_sensitivity_velocity(this->velocity_sensitivity(), this->solution_sensitivity());
    
    if (_ode_order > 1)
        update_sensitivity_acceleration(this->acceleration_sensitivity(), this->solution_sensitivity());
    
    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store-1; i>0; i--) {
        this->solution_sensitivity(i).zero();
        this->solution_sensitivity(i).add(1., this->solution_sensitivity(i-1));
        this->solution_sensitivity(i).close();
        
        this->velocity_sensitivity(i).zero();
        this->velocity_sensitivity(i).add(1., this->velocity_sensitivity(i-1));
        this->velocity_sensitivity(i).close();
        
        if (_ode_order > 1) {
            
            this->acceleration_sensitivity(i).zero();
            this->acceleration_sensitivity(i).add(1., this->acceleration_sensitivity(i-1));
            this->acceleration_sensitivity(i).close();
        }
    }
    
    // finally, update the system time
    sys.time          += dt;
    _first_sensitivity_step   = false;
}



void
MAST::TransientSolverBase::set_elem_data(unsigned int dim,
                                         const libMesh::Elem& ref_elem,
                                         MAST::GeomElem &elem) const {
    
    libmesh_assert(_assembly_ops);
    _assembly_ops->set_elem_data(dim, ref_elem, elem);
}


void
MAST::TransientSolverBase::init(const MAST::GeomElem &elem) {
    
    libmesh_assert(_assembly_ops);
    _assembly_ops->init(elem);
}




void
MAST::TransientSolverBase::clear_elem() {
    
    libmesh_assert(_assembly_ops);
    _assembly_ops->clear_elem();
}


