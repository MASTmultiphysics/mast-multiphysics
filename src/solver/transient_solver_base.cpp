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

// MAST includes
#include "solver/transient_solver_base.h"
#include "base/transient_assembly.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/linear_solver.h"



MAST::TransientSolverBase::TransientSolverBase():
dt(0.),
_first_step(true),
_assembly(nullptr),
_system(nullptr),
_if_highest_derivative_solution(false) {

}




MAST::TransientSolverBase::~TransientSolverBase() {
    this->clear_assembly();
}



void
MAST::TransientSolverBase::set_assembly(MAST::TransientAssembly& assembly) {
    
    // make sure that the previous assembly association has been cleared
    this->clear_assembly();
    
    _assembly = &assembly;
    _system   = &assembly.system();
    
    // number of time steps to store
    unsigned int n_iters = _n_iters_to_store();
    
    // now, add the vectors
    std::string nm;
    for (unsigned int i=0; i<n_iters; i++) {
        std::ostringstream iter;
        iter << i;
        
        // add the solution only for previous iterations. The current
        // solution is obtained from system->solution
        if (i>0) {
            
            nm = "transient_solution_";
            nm += iter.str();
            _system->add_vector(nm);
        }
        
        // add the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        _system->add_vector(nm);
        
        if (this->ode_order() > 1) {
            // add the acceleration
            nm = "transient_acceleration_";
            nm += iter.str();
            _system->add_vector(nm);
        }
    }
    
    _first_step = true;
}



void
MAST::TransientSolverBase::clear_assembly() {
    
    // if no system has been set so far, nothing needs to be done
    if (!_system)
        return;
    
    // clear the transient solutions stored in system for solution
    // number of time steps to store
    unsigned int n_iters = _n_iters_to_store();
    
    // now localize the vectors
    // add the vectors
    std::string nm;
    for (unsigned int i=0; i<n_iters; i++) {
        std::ostringstream iter;
        iter << i;
        
        // add the solution
        nm = "transient_solution_";
        nm += iter.str();
        if (_system->have_vector(nm))
            _system->remove_vector(nm);
        
        // add the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        if (_system->have_vector(nm))
            _system->remove_vector(nm);
        
        if (this->ode_order() > 1) {
            // add the acceleration
            nm = "transient_acceleration_";
            nm += iter.str();
            if (_system->have_vector(nm))
                _system->remove_vector(nm);
            
        }
    }
    
    _assembly   = nullptr;
    _system     = nullptr;
    _first_step = true;
}




libMesh::NumericVector<Real>&
MAST::TransientSolverBase::solution(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    if (prev_iter) {
        // get references to the solution
        std::string
        nm = "transient_solution_" + oss.str();

        return _system->get_vector(nm);
    }
    else
        return *_system->solution;
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::velocity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_velocity_" + oss.str();
    
    return _system->get_vector(nm);
}



libMesh::NumericVector<Real>&
MAST::TransientSolverBase::acceleration(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    std::ostringstream oss;
    oss << prev_iter;
    
    // get references to the solution
    std::string
    nm = "transient_acceleration_" + oss.str();
    
    return _system->get_vector(nm);
}




void
MAST::TransientSolverBase::solve_highest_derivative_and_advance_time_step() {
    
    libmesh_assert(_first_step);
    
    // tell the solver that the current solution being obtained is for the
    // highest time derivative
    _if_highest_derivative_solution = true;
    
    // Build the residual and Jacobian
    _system->assembly(true, true);

    // reset the solution flag
    _if_highest_derivative_solution = false;
    
    // The sensitivity problem is linear
    libMesh::LinearSolver<Real> * linear_solver = _system->get_linear_solver();
    
    std::pair<unsigned int, Real>
    solver_params = _system->get_linear_solve_parameters();
    
    // Solve the linear system.
    libMesh::SparseMatrix<Real> *
    pc = _system->request_matrix("Preconditioner");
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    dvec(velocity().zero_clone().release());

    std::pair<unsigned int, Real>
    rval = linear_solver->solve (*_system->matrix, pc,
                                 *dvec,
                                 *_system->rhs,
                                 solver_params.second,
                                 solver_params.first);

    libMesh::NumericVector<Real> *vec = nullptr;
    
    switch (this->ode_order()) {
            
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
    _system->get_dof_map().enforce_constraints_exactly
    (*_system, vec, /* homogeneous = */ true);
#endif
    
    _system->release_linear_solver(linear_solver);
    
    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store()-1; i>0; i--) {
        this->solution(i).zero();
        this->solution(i).add(1., this->solution(i-1));
        this->solution(i).close();
        
        this->velocity(i).zero();
        this->velocity(i).add(1., this->velocity(i-1));
        this->velocity(i).close();
        
        if (this->ode_order() > 1) {
            
            this->acceleration(i).zero();
            this->acceleration(i).add(1., this->acceleration(i-1));
            this->acceleration(i).close();
        }
    }
    
    // finally, update the system time
    _system->time     += dt;
    _first_step        = false;
}



void
MAST::TransientSolverBase::
build_local_quantities(const libMesh::NumericVector<Real>& current_sol,
                       std::vector<libMesh::NumericVector<Real>*>& sol) {

    // make sure that the system has been specified
    libmesh_assert(_system);

    // make sure there are no solutions in sol
    libmesh_assert(!sol.size());

    // resize the vector to store the quantities
    sol.resize(this->ode_order()+1);
    
    const std::vector<libMesh::dof_id_type>&
    send_list = _system->get_dof_map().get_send_list();
    

    for ( unsigned int i=0; i<=this->ode_order(); i++) {
        
        sol[i] = libMesh::NumericVector<Real>::build(_system->comm()).release();
        sol[i]->init(_system->n_dofs(),
                     _system->n_local_dofs(),
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
build_perturbed_local_quantities(const libMesh::NumericVector<Real>& current_dsol,
                                 std::vector<libMesh::NumericVector<Real>*>& sol) {
    
    // make sure that the system has been specified
    libmesh_assert(_system);

    // make sure there are no solutions in sol
    libmesh_assert(!sol.size());
    
    // resize the vector to store the quantities
    sol.resize(this->ode_order()+1);
    
    const std::vector<libMesh::dof_id_type>&
    send_list = _system->get_dof_map().get_send_list();
    
    
    std::unique_ptr<libMesh::NumericVector<Real> >
    tmp(this->solution().zero_clone().release());
    
    for ( unsigned int i=0; i<=this->ode_order(); i++) {
        
        sol[i] = libMesh::NumericVector<Real>::build(_system->comm()).release();
        sol[i]->init(_system->n_dofs(),
                     _system->n_local_dofs(),
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
                
                if (ode_order() == 1 && _if_highest_derivative_solution)
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
                
                if (ode_order() == 2 && _if_highest_derivative_solution)
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

    // first ask the solver to update the velocity and acceleration vector
    update_velocity(this->velocity(), *_system->solution);

    if (this->ode_order() > 1)
        update_acceleration(this->acceleration(), *_system->solution);

    // next, move all the solutions and velocities into older
    // time step locations
    for (unsigned int i=_n_iters_to_store()-1; i>0; i--) {
        this->solution(i).zero();
        this->solution(i).add(1., this->solution(i-1));
        this->solution(i).close();
        
        this->velocity(i).zero();
        this->velocity(i).add(1., this->velocity(i-1));
        this->velocity(i).close();
        
        if (this->ode_order() > 1) {
            
            this->acceleration(i).zero();
            this->acceleration(i).add(1., this->acceleration(i-1));
            this->acceleration(i).close();
        }
    }

    // and localize solution so that if the user calls the solution and
    // velocity routines,
    
    // finally, update the system time
    _system->time     += dt;
    _first_step        = false;
}




