
// MAST includes
#include "solver/transient_solver_base.h"
#include "base/transient_assembly.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"



MAST::TransientSolverBase::TransientSolverBase():
dt(0.),
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
    _system   = dynamic_cast<libMesh::TransientNonlinearImplicitSystem*>(&assembly.system());
    
    // number of time steps to store
    unsigned int n_iters = _n_iters_to_store();
    
    // resize the vectors if needed
    _solution.resize(n_iters);
    _velocity.resize(n_iters);
        
    for (unsigned int i=0; i<n_iters; i++) {
        _solution[i] = _system->solution->zero_clone().release();
        _velocity[i] = _system->solution->zero_clone().release();
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
    }
}



void
MAST::TransientSolverBase::clear_assembly() {
    _assembly = NULL;
    _system   = NULL;
    
    // also clear the vectors
    if (_solution.size()) {
        for (unsigned int i=0; i<_n_iters_to_store(); i++) {
            delete _solution[i]; _solution[i] = NULL;
            delete _velocity[i]; _velocity[i] = NULL;
        }
        
        _solution.resize(0);
        _velocity.resize(0);
    }
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



void
MAST::TransientSolverBase::_localize_solutions() {
    
    // make sure that the system has been specified
    libmesh_assert(_system);
    
    // number of time steps to store
    unsigned int n_iters = _n_iters_to_store();

    // send list
    const std::vector<libMesh::dof_id_type>& send_list =
    _system->get_dof_map().get_send_list();
    
    // now localize the vectors
    std::string nm;
    libMesh::NumericVector<Real> *vec;
    for (unsigned int i=0; i<n_iters; i++) {
        std::ostringstream iter;
        iter << i;
        
        // localize the solution
        nm = "transient_solution_";
        nm += iter.str();
        vec = &(_system->get_vector(nm));
        _solution[i]->init(_system->n_dofs(),
                           _system->n_local_dofs(),
                           send_list,
                           false,
                           libMesh::GHOSTED);
        vec->localize(*_solution[i], send_list);
        _solution[i]->close();

        
        // localize the velocity
        nm = "transient_velocity_";
        nm += iter.str();
        vec = &(_system->get_vector(nm));
        _velocity[i]->init(_system->n_dofs(),
                           _system->n_local_dofs(),
                           send_list,
                           false,
                           libMesh::GHOSTED);
        vec->localize(*_velocity[i], send_list);
        _velocity[i]->close();
    }
}



void
MAST::TransientSolverBase::advance_time_step() {

    // first, move all the solutions and velocities into older
    // time step locations
    
    std::string nm1, nm2;
    for (unsigned int i=_n_iters_to_store()-1; i>0; i--) {
        std::ostringstream oss1, oss2;
        oss1 << i;
        oss2 << i-1;
        
        // get references to the solution
        nm1 = "transient_solution_"; nm1 += oss1.str();
        nm2 = "transient_solution_"; nm2 += oss2.str();

        libMesh::NumericVector<Real>
        & old_sol     = _system->get_vector(nm1),
        & sol = _system->get_vector(nm2);

        // get references to the velocity
        nm1 = "transient_velocity_"; nm1 += oss1.str();
        nm2 = "transient_velocity_"; nm2 += oss2.str();
        
        libMesh::NumericVector<Real>
        & old_vel     = _system->get_vector(nm1),
        & vel = _system->get_vector(nm2);

        old_sol = sol;
        old_vel = vel;
        
        old_sol.close();
        old_vel.close();
    }
    
    // finally, copy the latest solution and update velocity
    libMesh::NumericVector<Real>
    &old_sol = _system->get_vector("transient_solution_0"),
    &sol     = _system->get_vector("transient_solution_1"),
    &vel     = _system->get_vector("transient_velocity_0");
    
    // now use finite difference to approximate the velocity
    sol = *(_system->solution);
    vel = sol;
    vel.add(-1., old_sol);
    vel.scale(1./dt);
    
    vel.close();
    
    // finally, update the system time
    _system->time += dt;
}




