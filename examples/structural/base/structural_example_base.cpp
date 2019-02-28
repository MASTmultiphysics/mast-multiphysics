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

// C++ includes
#include <sstream>

// MAST includes
#include "examples/structural/base/structural_example_base.h"
#include "examples/structural/base/thermal_stress_jacobian_scaling_function.h"
#include "examples/base/input_wrapper.h"
#include "examples/base/plot_results.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_system.h"
#include "base/nonlinear_implicit_assembly.h"
#include "base/eigenproblem_assembly.h"
#include "base/transient_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/element_property_card_base.h"
#include "solver/second_order_newmark_transient_solver.h"
#include "solver/slepc_eigen_solver.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_transient_assembly.h"
#include "elasticity/stress_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/fluid_structure_assembly_elem_operations.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/flutter_root_base.h"

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/serial_mesh.h"




MAST::Examples::StructuralExampleBase::
StructuralExampleBase(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::ExampleBase(comm_in),
_flutter_solver  (nullptr),
_flutter_root    (nullptr) {
    
}


MAST::Examples::StructuralExampleBase::~StructuralExampleBase() {
    
    
    if (_flutter_solver) delete _flutter_solver;
}



void
MAST::Examples::StructuralExampleBase::initialize_solution() {
    
    _sys->solution->zero();
}



void
MAST::Examples::StructuralExampleBase::_init_system_and_discipline() {
    
    // make sure that the mesh has been initialized
    libmesh_assert(_mesh);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system and set the preferences for structural
    // eigenvalue problems
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);

    // initialize the system to the right set of variables
    _sys_init       = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               _fetype);
    _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
}


void
MAST::Examples::StructuralExampleBase::_init_eq_sys() {
    
    _eq_sys->init();
    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _sys->set_exchange_A_and_B(true);
}




void
MAST::Examples::StructuralExampleBase::_init_dirichlet_conditions() {
    // should be overloaded in the derived classes
}



void
MAST::Examples::StructuralExampleBase::
_init_boundary_dirichlet_constraint(const unsigned int bid,
                                    const std::string& tag) {
    
    int
    constraint_dofs = (*_input)(_prefix+tag, "dofs to constrain on side", 123456), // by-default, constrain everything
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
MAST::Examples::StructuralExampleBase::_init_loads() {
    // should be overloaded in the derived classes
}


void
MAST::Examples::StructuralExampleBase::_init_material() {

    Real
    Eval      = (*_input)(_prefix+"E", "modulus of elasticity", 72.e9),
    rhoval    = (*_input)(_prefix+"rho", "material density", 2700.),
    nu_val    = (*_input)(_prefix+"nu", "Poisson's ratio",  0.33),
    kappa_val = (*_input)(_prefix+"kappa", "shear correction factor",  5./6.);
    
    
    MAST::Parameter
    *E         = new MAST::Parameter("E",          Eval),
    *rho       = new MAST::Parameter("rho",      rhoval),
    *nu        = new MAST::Parameter("nu",       nu_val),
    *kappa     = new MAST::Parameter("kappa", kappa_val);
    
    MAST::ConstantFieldFunction
    *E_f     = new MAST::ConstantFieldFunction(    "E",      *E),
    *rho_f   = new MAST::ConstantFieldFunction(  "rho",    *rho),
    *nu_f    = new MAST::ConstantFieldFunction(   "nu",     *nu),
    *kappa_f = new MAST::ConstantFieldFunction("kappa",  *kappa);
    
    this->add_parameter(*E);
    this->add_parameter(*rho);
    this->add_parameter(*nu);
    this->add_parameter(*kappa);
    this->register_field_function(*E_f);
    this->register_field_function(*rho_f);
    this->register_field_function(*nu_f);
    this->register_field_function(*kappa_f);
    
    _m_card = new MAST::IsotropicMaterialPropertyCard;
    _m_card->add(*E_f);
    _m_card->add(*rho_f);
    _m_card->add(*nu_f);
    _m_card->add(*kappa_f);
}


void
MAST::Examples::StructuralExampleBase::_init_section_property() {
    // should be overloaded in the derived classes
}


void
MAST::Examples::StructuralExampleBase::_init_pressure_load(bool on_side,
                                                           unsigned int id_num) {
    
    Real
    p_val    =  (*_input)(_prefix+"pressure", "pressure on side of domain",   2.e4);
    
    MAST::Parameter
    *press   = new MAST::Parameter( "p",  p_val);
    
    MAST::ConstantFieldFunction
    *press_f         = new MAST::ConstantFieldFunction("pressure", *press);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    
    p_load->add(*press_f);
    if (on_side)
        _discipline->add_side_load(id_num, *p_load);
    else
        _discipline->add_volume_load(id_num, *p_load);

    
    this->add_parameter(*press);
    this->register_field_function(*press_f);
    this->register_loading(*p_load);
    this->add_load_parameter(*press);
}


void
MAST::Examples::StructuralExampleBase::_init_temperature_load() {
    
    Real
    tval        =  (*_input)(_prefix+"temperature", "temperature on domain",       60.),
    alpha_val   =  (*_input)(_prefix+"alpha_expansion", "coefficient of thermal expansion", 2.5e-5);
    
    MAST::Parameter
    *temp            = new MAST::Parameter( "temperature",       tval),
    *alpha           = new MAST::Parameter(       "alpha",  alpha_val),
    &zero            = this->get_parameter("zero");
    
    MAST::ConstantFieldFunction
    *alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *alpha),
    *temp_f          = new MAST::ConstantFieldFunction("temperature",      *temp),
    *ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature",   zero);
    
    MAST::BoundaryConditionBase
    *T_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);

    MAST::FieldFunction<Real>
    *jac_scaling = new MAST::Examples::ThermalJacobianScaling;

    T_load->add(*temp_f);
    T_load->add(*ref_temp_f);
    T_load->add(*jac_scaling);
    _m_card->add(*alpha_f);

    _discipline->add_volume_load(0, *T_load);

    
    this->add_parameter(*alpha);
    this->add_parameter(*temp);
    this->register_field_function(*alpha_f);
    this->register_field_function(*temp_f);
    this->register_field_function(*ref_temp_f);
    this->register_field_function(*jac_scaling);
    this->register_loading(*T_load);
    this->add_load_parameter(*temp);
}


void
MAST::Examples::StructuralExampleBase::_init_piston_theory_load() {

    Real
    M                = (*_input)(_prefix+"mach", "Mach number", 3.),
    rho              = (*_input)(_prefix+"rho_fluid", "fluid density", 1.05),
    gamma            = (*_input)(_prefix+"gamma_fluid", "ratio of specific heats of fluid", 1.4);
    
    MAST::Parameter
    *velocity        = new MAST::Parameter("V"   ,     0.),
    *mach            = new MAST::Parameter("mach",      M),
    *rho_air         = new MAST::Parameter("rho" ,    rho),
    *gamma_air       = new MAST::Parameter("gamma", gamma);
    
    MAST::ConstantFieldFunction
    *velocity_f      = new MAST::ConstantFieldFunction("V",      *velocity),
    *mach_f          = new MAST::ConstantFieldFunction("mach",       *mach),
    *rho_air_f       = new MAST::ConstantFieldFunction("rho",     *rho_air),
    *gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *gamma_air);
    
    
    // now initialize the piston theory boundary conditions
    RealVectorX  vel = RealVectorX::Zero(3);
    vel(0)           = 1.;  // flow along the x-axis
    
    MAST::PistonTheoryBoundaryCondition
    *piston_bc       = new MAST::PistonTheoryBoundaryCondition(1,     // order
                                                               vel);  // vel vector
    piston_bc->add(*velocity_f);
    piston_bc->add(*mach_f);
    piston_bc->add(*rho_air_f);
    piston_bc->add(*gamma_air_f);
    
    _discipline->add_volume_load(0, *piston_bc);
    
    
    this->add_parameter(*velocity);
    this->add_parameter(*mach);
    this->add_parameter(*rho_air);
    this->add_parameter(*gamma_air);
    this->register_field_function(*velocity_f);
    this->register_field_function(*mach_f);
    this->register_field_function(*rho_air_f);
    this->register_field_function(*gamma_air_f);
    this->register_loading(*piston_bc);
}



void
MAST::Examples::StructuralExampleBase::static_solve() {
    
    libmesh_assert(_initialized);
    
    _eq_sys->print_info();
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false),
    nonlinear  = (*_input)(_prefix+"if_nonlinear", "flag to turn on/off nonlinear strain", false);
    
    unsigned int
    n_steps    = (*_input)(_prefix+"nonlinear_load_steps", "number of load steps in nonlinear problem", 20);
    
    if (!nonlinear) n_steps = 1;
    
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_static.exo";

    // create the nonlinear assembly object
    MAST::NonlinearImplicitAssembly                  assembly;
    MAST::StructuralNonlinearAssemblyElemOperations  elem_ops;
    MAST::StressAssembly                             stress_assembly;
    MAST::StressStrainOutputBase                     stress_elem_ops;
    stress_elem_ops.set_participating_elements_to_all();
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    MAST::Examples::ThermalJacobianScaling
    *jac_scaling = nullptr;
    
    if (this->has_field_function("thermal_jacobian_scaling")) {
        
        jac_scaling = &(dynamic_cast<MAST::Examples::ThermalJacobianScaling&>
        (this->get_field_function("thermal_jacobian_scaling")));
        jac_scaling->set_assembly(assembly);
    }

    // initialize the solution before solving
    this->initialize_solution();
    
    MAST::StructuralNearNullVectorSpace nsp;
    _sys->nonlinear_solver->nearnullspace_object = &nsp;
    
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // now iterate over the load steps
    for (unsigned int i=0; i<n_steps; i++) {
        
        if (jac_scaling)
            jac_scaling->set_enable(true);
        assembly.reset_residual_norm_history();
        
        // update the load values
        this->update_load_parameters((i+1.)/(1.*n_steps));
        
        // solve the system
        _sys->solve(elem_ops, assembly);
        
        // output, if asked
        if (output) {
            
            // update stress values
            stress_assembly.update_stress_strain_data(stress_elem_ops, *_sys->solution);
            
            // write the solution for visualization
            exodus_writer.write_timestep(output_name,
                                         *_eq_sys,
                                         i+1,
                                         (1.*i)/(1.*(n_steps-1)));
        }
    }

    if (jac_scaling)
        jac_scaling->set_enable(false);
    _sys->nonlinear_solver->nearnullspace_object = nullptr;
}




void
MAST::Examples::StructuralExampleBase::static_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);

    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);

    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_static_sens_"+p.name()+".exo";

    // create the nonlinear assembly object
    MAST::NonlinearImplicitAssembly                  assembly;
    MAST::StructuralNonlinearAssemblyElemOperations  elem_ops;
    MAST::StressAssembly                             stress_assembly;
    MAST::StressStrainOutputBase                     stress_elem_ops;
    stress_elem_ops.set_participating_elements_to_all();
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_aggregation_coefficients(2., 10., 2.e8);
    
    MAST::StructuralSystemInitialization
    &str_sys_init = dynamic_cast<MAST::StructuralSystemInitialization&>(*_sys_init);

    libMesh::System
    &stress_sys = str_sys_init.get_stress_sys();

    libMesh::NumericVector<Real>
    &dXdp      = _sys->add_sensitivity_solution(),
    &dsigma_dp = stress_sys.add_sensitivity_solution();
    
    /////////////////////////////////////////////////////////////////////
    //   sensitivity of solution
    /////////////////////////////////////////////////////////////////////
    MAST::StructuralNearNullVectorSpace nsp;
    _sys->nonlinear_solver->nearnullspace_object = &nsp;

    // zero the solution before solving
    // we are assuming that the nonlinear solve was completed before this.
    // So, we will reuse that matrix, as opposed to reassembling it
    _sys->sensitivity_solve(elem_ops, assembly, p, false);

    // evaluate output before calculation of sensitivity, since the
    // object requires the data
    assembly.calculate_output(*_sys->solution, stress_elem_ops);
    stress_elem_ops.output_total(); // this evaluates the total output
    assembly.calculate_output_direct_sensitivity(*_sys->solution,
                                                 dXdp,
                                                 p,
                                                 stress_elem_ops);
    libMesh::out << "dq/dp: direct: " << stress_elem_ops.output_sensitivity_total(p) << std::endl;
    
    /////////////////////////////////////////////////////////////////////
    // write the solution for visualization
    /////////////////////////////////////////////////////////////////////
    if (output) {
        
        /////////////////////////////////////////////////////////////////////
        //   sensitivity of stress for plotting
        /////////////////////////////////////////////////////////////////////
        stress_assembly.update_stress_strain_sensitivity_data(stress_elem_ops,
                                                              *_sys->solution,
                                                              dXdp,
                                                              p,
                                                              dsigma_dp);
        
        // swap solutions for output, since libMesh writes System::solution
        // to the output.
        _sys->solution->swap(dXdp);
        stress_sys.solution->swap(dsigma_dp);
        
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(output_name,
                                                            *_eq_sys);
        
        _sys->solution->swap(dXdp);
        stress_sys.solution->swap(dsigma_dp);
    }
    
    _sys->nonlinear_solver->nearnullspace_object = nullptr;
}




void
MAST::Examples::StructuralExampleBase::
static_adjoint_sensitivity_solve(//MAST::OutputAssemblyElemOperations& q,
                                 MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);

    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_adjoint.exo";

    // create the nonlinear assembly object
    MAST::NonlinearImplicitAssembly                  assembly;
    MAST::StructuralNonlinearAssemblyElemOperations  elem_ops;
    MAST::StressAssembly                             stress_assembly;
    MAST::StressStrainOutputBase                     stress_elem_ops;

    // we define the output as the p-norm stress over the whole structure
    stress_elem_ops.set_participating_elements_to_all();
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_aggregation_coefficients(2., 10., 2.e8);

    
    MAST::StructuralNearNullVectorSpace nsp;
    _sys->nonlinear_solver->nearnullspace_object = &nsp;

    
    // zero the solution before solving
    // we are assuming that the nonlinear solve was completed before this.
    // So, we will reuse that matrix, as opposed to reassembling it
    assembly.calculate_output(*_sys->solution, stress_elem_ops);
    libMesh::out << "output : " << stress_elem_ops.output_total() << std::endl;
    _sys->adjoint_solve(elem_ops, stress_elem_ops, assembly, false);
    
    Real dqdp =
    assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                  _sys->get_adjoint_solution(),
                                                  p,
                                                  elem_ops,
                                                  stress_elem_ops);
    libMesh::out << "dq/dp: adjoint: " << dqdp << std::endl;
    
    // write the solution for visualization
    if (output) {
        
        _sys->solution->swap(_sys->get_adjoint_solution(0));
        
        libMesh::ExodusII_IO(*_mesh).write_equation_systems(output_name,
                                                            *_eq_sys);
        
        _sys->solution->swap(_sys->get_adjoint_solution(0));
    }
    
    _sys->nonlinear_solver->nearnullspace_object = nullptr;
}




void
MAST::Examples::StructuralExampleBase::modal_solve(std::vector<Real>& eig) {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false),
    static_solve = (*_input)(_prefix+"modal_about_nonlinear_static", "modal eigensolution is performed about a nonlinear equilibrium", false);

    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_modal.exo";

    unsigned int
    n_req        = (*_input)(_prefix+"n_eig", "number of eigenvalues to compute", 10);
    _sys->set_n_requested_eigenvalues(n_req);
    
    // create the nonlinear assembly object
    MAST::EigenproblemAssembly                               assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations  elem_ops;
    _sys->initialize_condensed_dofs(*_discipline);
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    
    // perform static solution if requested
    if (static_solve) {
        
        libMesh::NumericVector<Real>&
        base_sol = _sys->add_vector("base_solution");
        this->static_solve();
        base_sol = *_sys->solution;
        assembly.set_base_solution(base_sol);
    }

    
    _sys->eigenproblem_solve(elem_ops, assembly);
    
    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    eig.resize(nconv);
    _basis.resize(nconv);

    libMesh::ExodusII_IO*
    writer = nullptr;
    
    if (output)
        writer = new libMesh::ExodusII_IO(*_mesh);
    
    for (unsigned int i=0; i<nconv; i++) {
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        
        std::ostringstream oss;
        oss << "mode_" << i;
        
        _basis[i] = &(_sys->add_vector(oss.str()));
        _sys->get_eigenpair(i, re, im, *_basis[i]);
        eig[i] = re;
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        if (output) {
            
            // We write the file in the ExodusII format.
            _sys->solution->swap(*_basis[i]);
            writer->write_timestep(output_name,
                                   *_eq_sys,
                                   i+1, i);
            _sys->solution->swap(*_basis[i]);
        }
    }
}





void
MAST::Examples::StructuralExampleBase::modal_sensitivity_solve(MAST::Parameter& p,
                                                               std::vector<Real>& deig_dp) {
    libmesh_assert(_initialized);

    bool
    static_solve = (*_input)(_prefix+"modal_about_nonlinear_static", "modal eigensolution is performed about a nonlinear equilibrium", false);

    // Get the number of converged eigen pairs.
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    deig_dp.resize(nconv);
    
    // create the modal assembly object
    MAST::EigenproblemAssembly                               assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations  elem_ops;
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    // perform static solution sensitivity if requested
    if (static_solve) {
        
        libMesh::NumericVector<Real>&
        base_sol = _sys->get_vector("base_solution");

        (*_sys->solution) = base_sol;
        this->static_sensitivity_solve(p);

        assembly.set_base_solution(base_sol);
        assembly.set_base_solution(_sys->get_sensitivity_solution(0), true);
    }
    
    _sys->eigenproblem_sensitivity_solve(elem_ops, assembly, p, deig_dp);
}




void
MAST::Examples::StructuralExampleBase::modal_solve_with_nonlinear_load_stepping() {
    
    libmesh_assert(_initialized);
    
    bool
    output       = (*_input)(_prefix+"if_output", "if write output to a file", false),
    nonlinear  = (*_input)(_prefix+"if_nonlinear", "flag to turn on/off nonlinear strain", false);

    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    static_output_name = output_name + "_static.exo";

    libmesh_assert(nonlinear); // make sure that the user asked for nonlinear
    
    // set the number of load steps
    unsigned int
    n_steps    = (*_input)(_prefix+"nonlinear_load_steps", "number of load steps in nonlinear problem", 20),
    n_perturb = 1,
    n_eig_req = (*_input)(_prefix+"n_eig", "number of eigenvalues to compute", 10);
    _sys->set_n_requested_eigenvalues(n_eig_req);
    
    std::ofstream freq;
    
    // write only on the zero-th processor
    if (this->comm().rank() == 0) {
        
        freq.open("freq.txt", std::ofstream::out);
        
        // write the header to the freq file
        freq
        << std::setw(10) << "Step"
        << std::setw(35) << "Load-Factor";
        for (unsigned int i=0; i<n_eig_req; i++) {
            std::ostringstream nm;
            nm << "mode_" << i;
            freq << std::setw(35) << nm.str();
        }
    }
    
    // create the assembly object
    MAST::NonlinearImplicitAssembly                          nonlinear_assembly;
    MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
    MAST::StressAssembly                                     stress_assembly;
    MAST::StressStrainOutputBase                             stress_elem_ops;
    MAST::EigenproblemAssembly                               modal_assembly;
    MAST::StructuralModalEigenproblemAssemblyElemOperations  modal_elem_ops;
    _sys->initialize_condensed_dofs(*_discipline);
    stress_elem_ops.set_participating_elements_to_all();

    
    nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    modal_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    modal_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    MAST::Examples::ThermalJacobianScaling
    *jac_scaling = nullptr;
    
    if (this->has_field_function("thermal_jacobian_scaling")) {
        
        jac_scaling = &(dynamic_cast<MAST::Examples::ThermalJacobianScaling&>
                        (this->get_field_function("thermal_jacobian_scaling")));
        jac_scaling->set_assembly(nonlinear_assembly);
    }

    
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    // writer for the modes
    std::vector<libMesh::ExodusII_IO*>
    mode_writer(n_eig_req);

    for (unsigned int i=0; i<n_eig_req; i++)
        mode_writer[i] = new libMesh::ExodusII_IO(*_mesh);
    
    //  store the base solution
    libMesh::NumericVector<Real>&
    base_sol = _sys->add_vector("base_solution");
    
    // initialize the solution before solving
    this->initialize_solution();
    
    // now iterate over the load steps
    for (int i_step=0; i_step<n_steps; i_step++) {

        if (jac_scaling)
            jac_scaling->set_enable(true);
        nonlinear_assembly.reset_residual_norm_history();

        Real
        eta = (i_step+1.)/(1.*n_steps);
        this->update_load_parameters(eta);
        
        libMesh::out
        << "Load step: " << i_step << "  Load-factor = " << eta << std::endl;

        if (this->comm().rank() == 0)
            freq << std::endl
            << std::setw(10) << i_step
            << std::setw(35) << std::fixed << std::setprecision(15) << eta;
        

        ///////////////////////////////////////////////////////////////////
        // nonlinear solution
        ///////////////////////////////////////////////////////////////////
        _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
        

        ///////////////////////////////////////////////////////////////////
        // stress update
        ///////////////////////////////////////////////////////////////////
        stress_assembly.update_stress_strain_data(stress_elem_ops, *_sys->solution);

        if (output) {
            exodus_writer.write_timestep(static_output_name,
                                         *_eq_sys,
                                         i_step+1,
                                         eta);
        }
        
        
        ///////////////////////////////////////////////////////////////////
        // solve for the natural frequencies
        //////////////////////////////////////////////////////////////////
        base_sol = *_sys->solution;
        
        modal_assembly.clear_base_solution();
        modal_assembly.set_base_solution(base_sol);
        _sys->eigenproblem_solve(modal_elem_ops, modal_assembly);
        
        // Get the number of converged eigen pairs.
        unsigned int
        nconv = std::min(_sys->get_n_converged_eigenvalues(),
                         _sys->get_n_requested_eigenvalues());

        for (unsigned int i=0; i<nconv; i++) {
            
            // now write the eigenvalue
            Real
            re = 0.,
            im = 0.;
            _sys->get_eigenpair(i, re, im, base_sol);
            
            libMesh::out
            << std::setw(35) << std::fixed << std::setprecision(15)
            << re << std::endl;

            if (this->comm().rank() == 0)
                freq << std::setw(35) << std::fixed << std::setprecision(15) << re;
            
            // if any of the eigenvalues is negative, then use the mode of the
            // smallest eigenvalue to perturb the deformation shape
            if (re < 0. && n_perturb == 0) {
                libMesh::out
                << "** Negative Frequency: Perturbing panel in first mode. **"
                << std::endl;
                
                *_sys->solution = base_sol;
                _sys->solution->scale(1./_sys->solution->linfty_norm()*1.e-1);
                n_perturb++;
            }
            
            if (output) {
                
                std::ostringstream file_name;
                
                // We write the file in the ExodusII format.
                file_name
                << output_name
                << "_mode_"
                << std::setw(3) << std::setfill('0') << std::right << i
                << ".exo";
                
                // We write the file in the ExodusII format.
                base_sol.swap(*_sys->solution);
                mode_writer[i]->write_timestep(file_name.str(),
                                               *_eq_sys,
                                               i_step+1,
                                               eta);
                base_sol.swap(*_sys->solution);
            }
        }
        
        // copy the solution back to base_sol
        base_sol = *_sys->solution;
    }
    
    // delete the exodus writes for the modes
    for (unsigned int i=0; i<n_eig_req; i++)
        delete mode_writer[i];

    if (jac_scaling)
        jac_scaling->set_enable(false);
}



void
MAST::Examples::StructuralExampleBase::transient_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output     = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient.exo";

    
    // create the nonlinear assembly object
    MAST::TransientAssembly                           assembly;
    MAST::StructuralTransientAssemblyElemOperations   elem_ops;
    MAST::StressAssembly                              stress_assembly;
    MAST::StressStrainOutputBase                      stress_elem_ops;
    MAST::SecondOrderNewmarkTransientSolver           solver;

    stress_elem_ops.set_participating_elements_to_all();
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
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
        
        // update stress values
        stress_assembly.update_stress_strain_data(stress_elem_ops, *_sys->solution);

        // update time value
        tval  += solver.dt;
        t_step++;
    }
}




void
MAST::Examples::StructuralExampleBase::transient_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    bool
    output                = (*_input)(_prefix+"if_output", "if write output to a file", false);
    std::string
    output_name           = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    transient_output_name = output_name + "_transient_sensitivity_" + p.name() + ".exo";

    // the output from analysis should have been saved for sensitivity
    libmesh_assert(output);
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                           assembly;
    MAST::StructuralTransientAssemblyElemOperations   elem_ops;
    MAST::StressAssembly                              stress_assembly;
    MAST::StressStrainOutputBase                      stress_elem_ops;
    MAST::SecondOrderNewmarkTransientSolver solver;

    stress_elem_ops.set_participating_elements_to_all();
    assembly.set_discipline_and_system(*_discipline, *_sys_init);
    elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    stress_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
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
        
        // update stress values
        MAST::StructuralSystemInitialization
        &str_sys_init = dynamic_cast<MAST::StructuralSystemInitialization&>(*_sys_init);
        
        stress_assembly.update_stress_strain_sensitivity_data(stress_elem_ops,
                                                              solver.solution(),
                                                              solver.solution_sensitivity(),
                                                              p,
                                                              *str_sys_init.get_stress_sys().solution);
        
        // update time value
        tval  += solver.dt;
        t_step++;
    }
    
}





void
MAST::Examples::StructuralExampleBase::piston_theory_flutter_solve() {

    libmesh_assert(_initialized);

    bool
    output              = (*_input)(_prefix+"if_output", "if write output to a file", false);

    std::string
    output_name         = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    flutter_mode_name   = output_name + "_flutter_mode.exo",
    flutter_output_name = output_name + "_flutter.txt";


    Real
    tol        = (*_input)(_prefix+"flutter_solver_tol", "tolerance of flutter root convergence", 1.e-5);
    
    unsigned int
    max_iters  = (*_input)(_prefix+"flutter_solver_max_search_iters", "maximum iterations for convergence of flutter root after cross-over identification", 10);
    
    // clear out the data structures of the flutter solver before
    // this solution
    _flutter_root = nullptr;
    if (this->comm().rank() == 0)
        _flutter_solver->set_output_file(flutter_output_name);
    
    
    // set the velocity of piston theory to zero for modal analysis
    MAST::Parameter
    &velocity = this->get_parameter("velocity");
    velocity = 0.;
    
    //////////////////////////////////////////////////////////////////
    // modal solution will provide the modal basis for flutter solution
    //////////////////////////////////////////////////////////////////
    std::vector<Real> eig; // this is a dummy argument, not used later
    this->modal_solve(eig);

    //////////////////////////////////////////////////////////////////
    // now initialize the flutter solver
    //////////////////////////////////////////////////////////////////
    if (_flutter_solver)
        _flutter_solver->clear();
    else
        _flutter_solver = new MAST::TimeDomainFlutterSolver;
    
    Real
    V_low    =  (*_input)(_prefix+"V_lower", "time-domain flutter solver search: lower limit on speed",  1.e3),
    V_up     =  (*_input)(_prefix+"V_upper", "time-domain flutter solver search: upper limit on speed", 1.2e3);
    unsigned int
    n_V_divs =  (*_input)(_prefix+"n_V_divs", "time-domain flutter solver search: number of divisions between speed limits",   10);
    
    MAST::StructuralFluidInteractionAssembly fsi_assembly;
    MAST::FluidStructureAssemblyElemOperations fsi_elem_ops;
    fsi_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    fsi_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    
    _flutter_solver->attach_assembly(fsi_assembly);
    _flutter_solver->initialize(velocity,
                                V_low,
                                V_up,
                                n_V_divs,
                                _basis);
    
    std::pair<bool, MAST::FlutterRootBase*>
    sol = _flutter_solver->analyze_and_find_critical_root_without_tracking(tol,
                                                                           max_iters);
    _flutter_solver->print_sorted_roots();
    _flutter_solver->clear_assembly_object();

    // make sure solution was found
    libmesh_assert(sol.first);
    _flutter_root = sol.second;
    
    if (sol.first && output) {
        
        MAST::plot_structural_flutter_solution(flutter_mode_name,
                                               *_sys,
                                               sol.second->eig_vec_right,
                                               _basis);
    }
}





void
MAST::Examples::StructuralExampleBase::
piston_theory_flutter_sensitivity_solve(MAST::Parameter& p) {
    
    libmesh_assert(_initialized);
    
    //Make sure that  a solution is available for sensitivity
    libmesh_assert(_flutter_root);
    
    MAST::Parameter
    &velocity = this->get_parameter("velocity");
    
    // flutter solver will need velocity to be defined as a parameter for
    // sensitivity analysis
    // initialize the flutter solver for sensitivity.
    MAST::StructuralFluidInteractionAssembly fsi_assembly;
    MAST::FluidStructureAssemblyElemOperations fsi_elem_ops;
    
    fsi_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    fsi_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);

    _flutter_solver->attach_assembly(fsi_assembly);

    _flutter_solver->calculate_sensitivity(*_flutter_root, p);

    _flutter_solver->clear_assembly_object();
}


