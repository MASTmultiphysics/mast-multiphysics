/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2015  Manav Bhatia
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
    
    // resize the vectors if needed
    _solution.resize(n_iters);
    _velocity.resize(n_iters);
    if (this->ode_order()>1)
        _acceleration.resize(n_iters);
        
    for (unsigned int i=0; i<n_iters; i++) {
        _solution[i] = _system->solution->zero_clone().release();
        _velocity[i] = _system->solution->zero_clone().release();
        if (this->ode_order()>1)
            _acceleration[i] = _system->solution->zero_clone().release();
    }
    
    // now localize the vectors
    // add the vectors
    std::string nm;
    for (unsigned int i=0; i<n_iters; i++) {
        std::ostringstream iter;
        iter << i;
        
        // add the solution
        nm = "transient_solution_";
        nm += iter.str();
        if (!_system->have_vector(nm))
            _system->add_vector(nm);
        
        // add the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        if (!_system->have_vector(nm))
            _system->add_vector(nm);
        
        if (this->ode_order() > 1) {
            // add the acceleration
            nm = "transient_acceleration_";
            nm += iter.str();
            if (!_system->have_vector(nm))
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
    
    if (_solution.size()) {
        
        for (unsigned int i=0; i<n_iters; i++) {
            delete _solution[i]; _solution[i] = NULL;
            delete _velocity[i]; _velocity[i] = NULL;
            if (this->ode_order()>1) {
                delete _acceleration[i]; _acceleration[i] = NULL;
            }
        }
        
        _solution.resize(0);
        _velocity.resize(0);
        if (this->ode_order()>1)
            _acceleration.resize(0);
    }

    
    _assembly   = NULL;
    _system     = NULL;
    _first_step = true;
}




const libMesh::NumericVector<Real>&
MAST::TransientSolverBase::solution(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    libmesh_assert_equal_to(_solution.size(), _n_iters_to_store());
    libmesh_assert(_solution[prev_iter]);
    
    return *(_solution[prev_iter]);
}



const libMesh::NumericVector<Real>&
MAST::TransientSolverBase::velocity(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    libmesh_assert_equal_to(_velocity.size(), _n_iters_to_store());
    libmesh_assert(_velocity[prev_iter]);
    
    return *(_velocity[prev_iter]);
}



const libMesh::NumericVector<Real>&
MAST::TransientSolverBase::acceleration(unsigned int prev_iter) const {
    
    // make sura that prev_iter is within acceptable bounds
    libmesh_assert_less(prev_iter, _n_iters_to_store());
    
    // make sure that the vectors have been initialized
    libmesh_assert_equal_to(_acceleration.size(), _n_iters_to_store());
    libmesh_assert(_acceleration[prev_iter]);
    
    return *(_acceleration[prev_iter]);
}



void
MAST::TransientSolverBase::
_localize_solutions(const libMesh::NumericVector<Real>& current_sol) {
    
    // make sure that the system has been specified
    libmesh_assert(_system);
    
    // number of time steps to store
    unsigned int n_iters = _n_iters_to_store();

    // send list
    const std::vector<libMesh::dof_id_type>& send_list =
    _system->get_dof_map().get_send_list();
    
    // now localize the vectors
    std::string nm;
    const libMesh::NumericVector<Real> *vec;
    for (unsigned int i=0; i<n_iters; i++) {
        
        std::ostringstream iter;
        iter << i;

        nm = "transient_solution_" + iter.str();

        ///////////////////////////////////////////
        // localize the solution
        ///////////////////////////////////////////
        if (i == 0) {
            // update both the localized and system stored solution
            _system->get_vector(nm) = current_sol;
            _system->get_vector(nm).close();
            vec = &current_sol;
        }
        else {
            // update the localized with the system stored solution
            vec = &(_system->get_vector(nm));
        }
        _solution[i]->init(_system->n_dofs(),
                           _system->n_local_dofs(),
                           send_list,
                           false,
                           libMesh::GHOSTED);
        vec->localize(*_solution[i], send_list);
        _solution[i]->close();

        ///////////////////////////////////////////
        // localize the velocity
        ///////////////////////////////////////////
        nm = "transient_velocity_" + iter.str();
        //
        if (i == 0) {
            // update the velocity using the current estimate
            _update_velocity(_system->get_vector(nm));
        }

        // now update localized vectors using the stored velocity
        vec = &(_system->get_vector(nm));
        _velocity[i]->init(_system->n_dofs(),
                           _system->n_local_dofs(),
                           send_list,
                           false,
                           libMesh::GHOSTED);
        vec->localize(*_velocity[i], send_list);
        _velocity[i]->close();
        
        if (this->ode_order()) {
            ///////////////////////////////////////////
            // localize the acceleration
            ///////////////////////////////////////////
            nm = "transient_acceleration_" + iter.str();
            //
            if (i == 0) {
                // update the acceleration using the current estimate
                _update_acceleration(_system->get_vector(nm));
            }
            
            // now update localized vectors using the stored acceleration
            vec = &(_system->get_vector(nm));
            _acceleration[i]->init(_system->n_dofs(),
                                   _system->n_local_dofs(),
                                   send_list,
                                   false,
                                   libMesh::GHOSTED);
            vec->localize(*_acceleration[i], send_list);
            _acceleration[i]->close();
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
    _update_velocity(_system->get_vector("transient_velocity_0"));
    if (this->ode_order() > 1)
        _update_acceleration(_system->get_vector("transient_acceleration_0"));

    // next, move all the solutions and velocities into older
    // time step locations
    std::string nm1, nm2;
    for (unsigned int i=_n_iters_to_store()-1; i>0; i--) {
        std::ostringstream oss1, oss2;
        oss1 << i;
        oss2 << i-1;
        
        // get references to the solution
        nm1 = "transient_solution_" + oss1.str();
        nm2 = "transient_solution_" + oss2.str();

        libMesh::NumericVector<Real>
        & old_sol     = _system->get_vector(nm1),
        & sol         = _system->get_vector(nm2);

        old_sol = sol;
        old_sol.close();
        
        // get references to the velocity
        nm1 = "transient_velocity_" + oss1.str();
        nm2 = "transient_velocity_" + oss2.str();
        
        libMesh::NumericVector<Real>
        & old_vel     = _system->get_vector(nm1),
        & vel         = _system->get_vector(nm2);

        old_vel = vel;
        old_vel.close();
        
        if (this->ode_order() > 1) {
            // get references to the acceleration
            nm1 = "transient_acceleration_" + oss1.str();
            nm2 = "transient_acceleration_" + oss2.str();
            
            libMesh::NumericVector<Real>
            & old_accel     = _system->get_vector(nm1),
            & accel         = _system->get_vector(nm2);
            
            old_accel = accel;
            old_accel.close();
        }
        

    }

    // and localize solution so that if the user calls the solution and
    // velocity routines,
    
    // finally, update the system time
    _system->time += dt;
    _first_step    = false;
}




