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
#include "solver/transient_solver_base.h"
#include "base/transient_assembly.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"



MAST::TransientSolverBase::TransientSolverBase():
dt(0.),
_first_step(false),
_assembly(NULL),
_system(NULL) {

}




MAST::TransientSolverBase::~TransientSolverBase() {
    this->clear_assembly();
}



void
MAST::TransientSolverBase::set_assembly(MAST::TransientAssembly& assembly) {
    
    // make sure that the previous assembly association has been cleared
    this->clear_assembly();
    
    _assembly = &assembly;
    _system   = dynamic_cast<libMesh::NonlinearImplicitSystem*>(&assembly.system());
    
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
    
    _assembly   = NULL;
    _system     = NULL;
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
MAST::TransientSolverBase::
build_local_quantities(const libMesh::NumericVector<Real>& current_sol,
                       std::vector<libMesh::NumericVector<Real>*>& sol) {

    // make sure that the system has been specified
    libmesh_assert(_system);

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
                
                // copy the solution to this
                current_sol.localize(*sol[i], send_list);
            }
                break;

            case 1: {
                
                // update the current local velocity vector
                libMesh::NumericVector<Real>& vel = this->velocity();

                // calculate the velocity and localize it
                _update_velocity(vel, current_sol);
                vel.localize(*sol[i], send_list);
            }
                break;

            case 2: {
                
                // update the current local acceleration vector
                libMesh::NumericVector<Real>& acc = this->acceleration();

                // calculate the acceleration and localize it
                _update_acceleration(acc, current_sol);
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
MAST::TransientSolverBase::set_initial_condition(Real val) {
    
    libmesh_assert(_system);
    libmesh_assert(_first_step);
    
    std::string nm;
    for (unsigned int i=0; i<_n_iters_to_store(); i++) {
        std::ostringstream oss;
        oss << i;
        nm = "transient_solution_" + oss.str();
        libMesh::NumericVector<Real> & sol =  _system->get_vector(nm);
        sol = val;
        sol.close();
    }

    _system->solution->operator=(val);
    _system->solution->close();
}


void
MAST::TransientSolverBase::advance_time_step() {


    // first ask the solver to update the acceleration vector
    _update_velocity(this->velocity(), *_system->solution);
    if (this->ode_order() > 1)
        _update_acceleration(this->acceleration(), *_system->solution);

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




