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
#include "examples/thermal/base/thermal_example_base.h"
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
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "heat_conduction/heat_conduction_transient_assembly.h"

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/serial_mesh.h"




MAST::Examples::ThermalExampleBase::
ThermalExampleBase(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::ExampleBase(comm_in) {
    
}


MAST::Examples::ThermalExampleBase::~ThermalExampleBase() {
    
    
}



void
MAST::Examples::ThermalExampleBase::initialize_solution() {
    
    _sys->solution->zero();
}



void
MAST::Examples::ThermalExampleBase::_init_system_and_discipline() {
    
    // make sure that the mesh has been initialized
    libmesh_assert(_mesh);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system and set the preferences for structural
    // eigenvalue problems
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);
    
    // initialize the system to the right set of variables
    _sys_init       = new MAST::HeatConductionSystemInitialization(*_sys,
                                                                   _sys->name(),
                                                                   _fetype);
    _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
}


void
MAST::Examples::ThermalExampleBase::_init_eq_sys() {
    
    _eq_sys->init();
}




void
MAST::Examples::ThermalExampleBase::_init_dirichlet_conditions() {
    // should be overloaded in the derived classes
}



void
MAST::Examples::ThermalExampleBase::
_init_boundary_dirichlet_constraint(const unsigned int bid,
                                    const std::string& tag) {
    
    int
    constraint_dofs = (*_input)(_prefix+tag, "dofs to constrain on side", 1), // by-default, constrain temperature
    dof             = 0;
    
    // nothing to constrain if the value is negative
    if (constraint_dofs < 0) return;
    
    // now, use this to add numbers
    std::vector<unsigned int>
    constr_dofs;
    
    while (constraint_dofs > 0) {
        dof             = constraint_dofs%10;  // gives the last integer
        constraint_dofs = constraint_dofs/10;  // removes the last integrer
        constr_dofs.push_back(dof-1);
    }
    
    // remove duplicates, if any
    std::sort(constr_dofs.begin(), constr_dofs.end());
    constr_dofs.erase(std::unique(constr_dofs.begin(), constr_dofs.end()), constr_dofs.end());
    
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;
    dirichlet->init (bid,  constr_dofs);
    _discipline->add_dirichlet_bc(bid,  *dirichlet);
    
    
    this->register_loading(*dirichlet);
}




void
MAST::Examples::ThermalExampleBase::_init_loads() {
    // should be overloaded in the derived classes
}


void
MAST::Examples::ThermalExampleBase::_init_material() {
    
    Real
    kval      = (*_input)(_prefix+"k", "thermal conductivity",  190.),
    rhoval    = (*_input)(_prefix+"rho", "material density",   2700.),
    cpval     = (*_input)(_prefix+"cp", "thermal capacitance",  864.);
    
    
    MAST::Parameter
    *k         = new MAST::Parameter("k",          kval),
    *rho       = new MAST::Parameter("rho",      rhoval),
    *cp        = new MAST::Parameter("cp",        cpval);
    
    MAST::ConstantFieldFunction
    *k_f     = new MAST::ConstantFieldFunction( "k_th",      *k),
    *rho_f   = new MAST::ConstantFieldFunction(  "rho",    *rho),
    *cp_f    = new MAST::ConstantFieldFunction(   "cp",     *cp);
    
    this->add_parameter(*k);
    this->add_parameter(*rho);
    this->add_parameter(*cp);
    this->register_field_function(*k_f);
    this->register_field_function(*rho_f);
    this->register_field_function(*cp_f);
    
    _m_card = new MAST::IsotropicMaterialPropertyCard;
    _m_card->add(*k_f);
    _m_card->add(*rho_f);
    _m_card->add(*cp_f);
}


void
MAST::Examples::ThermalExampleBase::_init_section_property() {
    // should be overloaded in the derived classes
}


void
MAST::Examples::ThermalExampleBase::
_init_flux_load(bool on_side, unsigned int id_num) {
    
    Real
    q_val    =  (*_input)(_prefix+"flux", "flux load on side of domain",   2.e1);
    
    MAST::Parameter
    *q       = new MAST::Parameter( "flux",  q_val);
    
    MAST::ConstantFieldFunction
    *q_f     = new MAST::ConstantFieldFunction("heat_flux", *q);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *q_load  = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
    q_load->add(*q_f);

    if (on_side)
        _discipline->add_side_load(id_num, *q_load);
    else
        _discipline->add_volume_load(id_num, *q_load);

    this->add_parameter(*q);
    this->register_field_function(*q_f);
    this->register_loading(*q_load);
}



void
MAST::Examples::ThermalExampleBase::
_init_convection_load(bool on_side, unsigned int id_num) {
    
    Real
    h_val    =  (*_input)(_prefix+"convection_coeff", "coefficient of thermal convection",   1.e-2),
    T_val    =  (*_input)(_prefix+"convection_T_inf", "ambient temperature for thermal convection",   30.0);

    MAST::Parameter
    *h       = new MAST::Parameter( "h_coeff_convection",  h_val),
    *T       = new MAST::Parameter( "T_inf_convection",    T_val);

    MAST::ConstantFieldFunction
    *h_f     = new MAST::ConstantFieldFunction("convection_coeff", *h),
    *T_f     = new MAST::ConstantFieldFunction("ambient_temperature", *T);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *q_load  = new MAST::BoundaryConditionBase(MAST::CONVECTION_HEAT_FLUX);
    q_load->add(*h_f);
    q_load->add(*T_f);
    if (on_side)
        _discipline->add_side_load(id_num, *q_load);
    else
        _discipline->add_volume_load(id_num, *q_load);

    this->add_parameter(*h);
    this->add_parameter(*T);
    this->register_field_function(*h_f);
    this->register_field_function(*T_f);
    this->register_loading(*q_load);
}



void
MAST::Examples::ThermalExampleBase::
_init_radiation_load(bool on_side, unsigned int id_num) {
    
    Real
    s_val    =  (*_input)(_prefix+"sb_constant", "Stefan-Boltzmann constant",             5.670367e-8),
    e_val    =  (*_input)(_prefix+"emissivity", "radiation emissivity",                          0.85),
    T_val    =  (*_input)(_prefix+"radiation_T_inf", "ambient temperature for thermal radiation",  0.),
    T_zero   =  (*_input)(_prefix+"absolute_zero_T", "absolute zero temperature",                273.);
    
    MAST::Parameter
    *s       = new MAST::Parameter( "stefan_bolzmann_constant",     s_val),
    *e       = new MAST::Parameter( "emissivity",                   e_val),
    *T       = new MAST::Parameter( "ambient_temperature",          T_val),
    *T0      = new MAST::Parameter( "reference_zero_temperature",  T_zero);

    MAST::ConstantFieldFunction
    *e_f     = new MAST::ConstantFieldFunction("emissivity", *e);

    // initialize the load
    MAST::BoundaryConditionBase
    *q_load  = new MAST::BoundaryConditionBase(MAST::SURFACE_RADIATION_HEAT_FLUX);
    q_load->add(*e_f);
    q_load->add(*s);
    q_load->add(*T);
    q_load->add(*T0);
    
    if (on_side)
        _discipline->add_side_load(id_num, *q_load);
    else
        _discipline->add_volume_load(id_num, *q_load);

    this->add_parameter(*s);
    this->add_parameter(*e);
    this->add_parameter(*T);
    this->add_parameter(*T0);
    this->register_field_function(*e_f);
    this->register_loading(*q_load);
}



void
MAST::Examples::ThermalExampleBase::_init_source_load(unsigned int   domain_num) {
    
    Real
    q_val    =  (*_input)(_prefix+"heat_source", "heat source on domain",   2.e1);
    
    MAST::Parameter
    *q       = new MAST::Parameter( "q",  q_val);
    
    MAST::ConstantFieldFunction
    *q_f     = new MAST::ConstantFieldFunction("heat_source", *q);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *q_load  = new MAST::BoundaryConditionBase(MAST::HEAT_SOURCE);
    q_load->add(*q_f);
    
    _discipline->add_volume_load(domain_num, *q_load);
    
    this->add_parameter(*q);
    this->register_field_function(*q_f);
    this->register_loading(*q_load);
}





void
MAST::Examples::ThermalExampleBase::steady_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);
    
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_static.exo";
    
    // create the nonlinear assembly object
    MAST::NonlinearImplicitAssembly                      assembly;
    MAST::HeatConductionNonlinearAssemblyElemOperations  elem_ops;
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    
    
    // initialize the solution before solving
    this->initialize_solution();
    
    
    
    // solve the system
    _sys->solve(elem_ops, assembly);
    
    // output, if asked
    if (output)
        // write the solution for visualization
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(output_name, *_eq_sys);
    
    _sys->nonlinear_solver->nearnullspace_object = nullptr;
}




void
MAST::Examples::ThermalExampleBase::steady_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);
    
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_static_sens_"+p.name()+".exo";
    
    // create the nonlinear assembly object
    MAST::NonlinearImplicitAssembly                  assembly;
    MAST::HeatConductionNonlinearAssemblyElemOperations  elem_ops;
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    
    libMesh::NumericVector<Real>
    &dXdp      = _sys->add_sensitivity_solution();
    
    /////////////////////////////////////////////////////////////////////
    //   sensitivity of solution
    /////////////////////////////////////////////////////////////////////
    
    // zero the solution before solving
    // we are assuming that the nonlinear solve was completed before this.
    // So, we will reuse that matrix, as opposed to reassembling it
    _sys->sensitivity_solve(elem_ops, assembly, p, false);
    
    /////////////////////////////////////////////////////////////////////
    // write the solution for visualization
    /////////////////////////////////////////////////////////////////////
    if (output) {
        
        // swap solutions for output, since libMesh writes System::solution
        // to the output.
        _sys->solution->swap(dXdp);
        
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(output_name, *_eq_sys);
        
        _sys->solution->swap(dXdp);
    }
}






void
MAST::Examples::ThermalExampleBase::transient_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient.exo";
    
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                               assembly;
    MAST::HeatConductionTransientAssemblyElemOperations   elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                solver;
    
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_discipline_and_system(*_discipline, *_sys_init);
    solver.set_elem_operation_object(elem_ops);
    
    // initialize the solution to zero, or to something that the
    // user may have provided
    this->initialize_solution();
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO transient_output(*_mesh);
    
    
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
    solver.solve_highest_derivative_and_advance_time_step(assembly);
    
    // loop over time steps
    while (t_step < n_steps) {
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << tval
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
        
        // solve for the time-step
        solver.solve(assembly);
        solver.advance_time_step();
        
        // update time value
        tval  += solver.dt;
        t_step++;
    }
}




void
MAST::Examples::ThermalExampleBase::transient_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    bool
    output                = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name           = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo";
    
    // the output from analysis should have been saved for sensitivity
    libmesh_assert(output);
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                               assembly;
    MAST::HeatConductionTransientAssemblyElemOperations   elem_ops;
    MAST::FirstOrderNewmarkTransientSolver                solver;
    
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


