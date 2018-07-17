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

// C++ includes
#include <sstream>

// MAST includes
#include "examples/fluid/base/fluid_example_base.h"
#include "examples/base/input_wrapper.h"
#include "examples/base/plot_results.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/transient_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/element_property_card_base.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "fluid/flight_condition.h"
#include "fluid/integrated_force_output.h"

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/serial_mesh.h"




MAST::Examples::FluidExampleBase::
FluidExampleBase(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::ExampleBase(comm_in),
_flight_cond        (nullptr) {
    
}


MAST::Examples::FluidExampleBase::~FluidExampleBase() {
    
    if (_flight_cond) delete _flight_cond;
}


void
MAST::Examples::FluidExampleBase::init(MAST::Examples::GetPotWrapper& input,
                                       const std::string& prefix) {
    
    libmesh_assert(!_initialized);
    
    _prefix      = prefix;
    _input       = &input;
    
    // call the initialization routines for each component
    _init_fetype();
    _init_mesh();
    _init_system_and_discipline();
    _init_loads();
    _init_eq_sys();
    _init_material();
    
    _initialized = true;
}


void
MAST::Examples::FluidExampleBase::initialize_solution() {
    
    unsigned int
    n = _mesh->mesh_dimension();
    RealVectorX s = RealVectorX::Zero(n+2);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    if (n > 2)
        s(3) = _flight_cond->rho_u3();
    s(n+1) = _flight_cond->rho_e();
    
    _sys_init->initialize_solution(s);
}



void
MAST::Examples::FluidExampleBase::_init_system_and_discipline() {
    
    // make sure that the mesh has been initialized
    libmesh_assert(_mesh);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system and set the preferences for structural
    // eigenvalue problems
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    
    _sys_init       = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                      _sys->name(),
                                                                      _fetype,
                                                                      _mesh->mesh_dimension());
    _discipline     = new MAST::ConservativeFluidDiscipline(*_eq_sys);
}




void
MAST::Examples::FluidExampleBase::_init_eq_sys() {
    
    _eq_sys->init();
}




void
MAST::Examples::FluidExampleBase::
_init_boundary_conditions(const std::vector<unsigned int>& slip,
                          const std::vector<unsigned int>& no_slip,
                          const std::vector<unsigned int>& symm,
                          const std::vector<unsigned int>& far_field) {
    
    for (unsigned int i=0; i<slip.size(); i++) {
        
        MAST::BoundaryConditionBase* bc = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
        _discipline->add_side_load(slip[i], *bc);
        this->register_loading(*bc);
    }
    
    
    for (unsigned int i=0; i<symm.size(); i++) {
        
        MAST::BoundaryConditionBase* bc = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
        _discipline->add_side_load(symm[i], *bc);
        this->register_loading(*bc);
    }
    
    
    for (unsigned int i=0; i<far_field.size(); i++) {
        
        MAST::BoundaryConditionBase* bc = new MAST::BoundaryConditionBase(MAST::FAR_FIELD);
        _discipline->add_side_load(far_field[i], *bc);
        this->register_loading(*bc);
    }
    

    if (no_slip.size()) {

        std::vector<unsigned int> constrained_vars(2);
        constrained_vars[0] = _sys_init->vars()[1];
        constrained_vars[1] = _sys_init->vars()[2];
        if (_mesh->mesh_dimension() > 2)
            constrained_vars.push_back(_sys_init->vars()[3]);

        for (unsigned int i=0; i<no_slip.size(); i++) {
            
            MAST::DirichletBoundaryCondition* bc = new MAST::DirichletBoundaryCondition;
            bc->init(no_slip[i], constrained_vars);
            _discipline->add_dirichlet_bc(no_slip[i], *bc);
            this->register_loading(*bc);
        }
        _discipline->init_system_dirichlet_bc(*_sys);
    }
}




void
MAST::Examples::FluidExampleBase::_init_material() {
    
    _flight_cond    =  new MAST::FlightCondition;

    _flight_cond->flow_unit_vector(0)  =
    (*_input)(_prefix+"flow_unit_vector", "unit vector defining direction of flow", 1., 0);
    _flight_cond->flow_unit_vector(1)  =
    (*_input)(_prefix+"flow_unit_vector", "unit vector defining direction of flow", 0., 1);
    _flight_cond->flow_unit_vector(2)  =
    (*_input)(_prefix+"flow_unit_vector", "unit vector defining direction of flow", 0., 2);

    _flight_cond->ref_chord       =
    (*_input)(_prefix+"ref_c", "reference chord used in flutter analysis", 1.);
    _flight_cond->mach            =
    (*_input)(_prefix+"mach", "ambient flow Mach number for flow computation", .5);
    _flight_cond->gas_property.cp =
    (*_input)(_prefix+"cp", "ambient fluid heat capacity at constant pressure", 1003.);
    _flight_cond->gas_property.cv =
    (*_input)(_prefix+"cv", "ambient fluid heat capacity at constant volume", 716.);
    _flight_cond->gas_property.T  =
    (*_input)(_prefix+"temp", "ambient fluid temperature", 300.);
    _flight_cond->gas_property.rho=
    (*_input)(_prefix+"rho", "ambient fluid density", 1.05);
    _flight_cond->gas_property.if_viscous =
    (*_input)(_prefix+"if_viscous", "if the flow analysis should include viscosity", false);

    _flight_cond->init();

    // tell the discipline about the fluid values
    dynamic_cast<MAST::ConservativeFluidDiscipline*>(_discipline)->set_flight_condition(*_flight_cond);
}




void
MAST::Examples::FluidExampleBase::transient_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient.exo";
    
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                                  assembly;
    MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                   solver;
    RealVectorX
    nvec = RealVectorX::Zero(3);
    nvec(1) = 1.;
    MAST::IntegratedForceOutput                              force(nvec);
    std::set<libMesh::boundary_id_type> bids;
    bids.insert(3);
    force.set_participating_boundaries(bids);

    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    force.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_elem_operation_object(elem_ops);
    
    // initialize the solution to zero, or to something that the
    // user may have provided
    this->initialize_solution();
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO transient_output(*_mesh);
    std::ofstream force_output;
    force_output.open("force.txt");
    force_output
            << std::setw(10) << "t"
            << std::setw(30) << "force" << std::endl;
    
    // time solver parameters
    Real
    factor   = 0.,
    min_fac  = 1.5,
    vel_0    = 0.,
    vel_1    = 1.e12,
    p        = 0.5,
    tval     = 0.,
    max_dt   = (*_input)(_prefix+"max_dt", "maximum time-step size", 1.e-1);
    
    unsigned int
    t_step            = 0,
    iter_count_dt     = 0,
    n_iters_change_dt = (*_input)(_prefix+"n_iters_change_dt", "number of time-steps before dt is changed", 4),
    n_steps           = (*_input)(_prefix+"n_transient_steps", "number of transient time-steps", 100);
    solver.dt         = (*_input)(_prefix+"dt", "time-step size",    1.e-3);
    libMesh::out << "q_dyn = " << _flight_cond->q0() << std::endl;
    
    // ask the solver to update the initial condition for d2(X)/dt2
    // This is recommended only for the initial time step, since the time
    // integration scheme updates the velocity and acceleration at
    // each subsequent iterate
    solver.solve_highest_derivative_and_advance_time_step(assembly);
    
    // loop over time steps
    while (t_step < n_steps) {

        if (iter_count_dt == n_iters_change_dt) {

            libMesh::out
            << "Changing dt:  old dt = " << solver.dt
            << "    new dt = " ;

            factor        = std::pow(vel_0/vel_1, p);
            factor        = std::max(factor, min_fac);
            solver.dt     = std::min(solver.dt*factor, max_dt);

            libMesh::out << solver.dt << std::endl;

            iter_count_dt = 0;
            vel_0         = vel_1;
        }
        else
            iter_count_dt++;

        libMesh::out
        << "Time step: "    << t_step
        << " :  t = "       << tval
        << " :  dt = "      << solver.dt
        << " :  xdot-L2 = " << solver.velocity().l2_norm()
        << std::endl;
        
        // write the time-step
        if (output) {
            
            transient_output.write_timestep(transient_output_name,
                                            *_eq_sys,
                                            t_step+1,
                                            _sys->time);
            std::ostringstream oss;
            oss << output_name << "_sol_t_" << t_step;
            _sys->write_out_vector(*_sys->solution, "data", oss.str(), true);
        }
        
        // calculate the output quantity
        force.zero_for_analysis();
        assembly.calculate_output(solver.solution(), force);
        force_output
                << std::setw(10) << tval
                << std::setw(30) << force.output_total() << std::endl;

        //_sys->adjoint_solve(solver, force, assembly, true);
        //_sys->solution->swap(_sys->get_adjoint_solution(0));
        //adjoint_output.write_timestep("adjoint.exo", *_eq_sys, t_step+1, _sys->time);
        //_sys->solution->swap(_sys->get_adjoint_solution(0));
        
        // solve for the time-step
        solver.solve(assembly);
        solver.advance_time_step();
        
        // update time value
        tval  += solver.dt;
        t_step++;
    }
}




void
MAST::Examples::FluidExampleBase::transient_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    bool
    output                = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name           = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo";
    
    // the output from analysis should have been saved for sensitivity
    libmesh_assert(output);
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                                  assembly;
    MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                   solver;
    
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_elem_operation_object(elem_ops);
    
    // initialize the solution to zero, or to something that the
    // user may have provided
    //this->initialize_sensitivity_solution();
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // time solver parameters
    Real
    tval     = 0.;
    
    unsigned int
    t_step            = 0,
    n_steps           = (*_input)(_prefix+"n_transient_steps", "number of transient time-steps", 100);
    solver.dt         = (*_input)(_prefix+"dt", "time-step size",    1.e-3);
    
    
    // ask the solver to update the initial condition for d2(X)/dt2
    // This is recommended only for the initial time step, since the time
    // integration scheme updates the velocity and acceleration at
    // each subsequent iterate
    solver.solve_highest_derivative_and_advance_time_step_with_sensitivity(assembly, p);
    
    // loop over time steps
    while (t_step < n_steps) {
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << tval
        << " :  xdot-L2 = " << solver.velocity().l2_norm()
        << std::endl;
        
        // write the time-step
        if (output) {
            
            _sys->solution->swap(solver.solution_sensitivity());
            exodus_writer.write_timestep(transient_output_name,
                                         *_eq_sys,
                                         t_step+1,
                                         _sys->time);
            _sys->solution->swap(solver.solution_sensitivity());
        }
        
        std::ostringstream oss;
        oss << output_name << "_sol_t_" << t_step;
        _sys->read_in_vector(*_sys->solution, "data", oss.str(), true);
        
        // solve for the sensitivity time-step
        solver.sensitivity_solve(assembly, p);
        solver.advance_time_step_with_sensitivity();
        
        // update time value
        tval  += solver.dt;
        t_step++;
    }
    
}


