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


// MAST includes
#include "examples/structural/thermoelasticity_topology_optim_2D/thermoelasticity_topology_optim_2D.h"
#include "examples/base/input_wrapper.h"
#include "level_set/level_set_discipline.h"
#include "level_set/level_set_system_initialization.h"
#include "level_set/level_set_eigenproblem_assembly.h"
#include "level_set/level_set_transient_assembly.h"
#include "level_set/level_set_nonlinear_implicit_assembly.h"
#include "level_set/level_set_reinitialization_transient_assembly.h"
#include "level_set/level_set_volume_output.h"
#include "level_set/level_set_boundary_velocity.h"
#include "level_set/indicator_function_constrain_dofs.h"
#include "level_set/level_set_constrain_dofs.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/stress_output_base.h"
#include "elasticity/stress_assembly.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/stress_temperature_adjoint.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "property_cards/material_property_card_base.h"
#include "property_cards/element_property_card_2D.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/petsc_nonlinear_solver.h"



namespace MAST {
    namespace Examples {
        class TemperatureFunction:
        public MAST::FieldFunction<Real> {
        public:
            
            TemperatureFunction(MAST::SystemInitialization& sys, const std::string& nm):
            MAST::FieldFunction<Real>(nm),
            _function(new MAST::MeshFieldFunction(sys, nm)) { }
            virtual ~TemperatureFunction() { delete _function;}
            void init(const libMesh::NumericVector<Real>& sol) { _function->init(sol);}
            void clear() { _function->clear();}
            virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {
                RealVectorX vec; (*_function)(p, t, vec); v = vec(0);
            }
            virtual void derivative (const MAST::FunctionBase& f,
                                     const libMesh::Point& p,
                                     const Real t,
                                     Real& v) const {
                
                v = 0.;
            }
        protected:
            
            MAST::MeshFieldFunction *_function;
        };
    }
}



MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
ThermoelasticityTopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::TopologyOptimizationLevelSet2D(comm_in),
_conduction_sys        (nullptr),
_conduction_sys_init   (nullptr),
_conduction_discipline (nullptr),
_temp_function  (nullptr) {
    
}



MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
~ThermoelasticityTopologyOptimizationLevelSet2D() {
    
    if (!_initialized)
        return;
    
    delete _conduction_sys_init;
    delete _conduction_discipline;
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
evaluate(const std::vector<Real>& dvars,
         Real& obj,
         bool eval_obj_grad,
         std::vector<Real>& obj_grad,
         std::vector<Real>& fvals,
         std::vector<bool>& eval_grads,
         std::vector<Real>& grads) {
    
    
    libMesh::out << "New Evaluation" << std::endl;
    
    // copy DVs to level set function
    for (unsigned int i=0; i<_n_vars; i++)
        _level_set_sys->solution->set(_dv_params[i].first, dvars[i]);
    _level_set_sys->solution->close();
    _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
    _sys->solution->zero();
    
    /**********************************************************************
     * DO NOT zero out the gradient vector, since GCMMA needs it for the  *
     * subproblem solution                                                *
     **********************************************************************/
    MAST::LevelSetNonlinearImplicitAssembly                  nonlinear_assembly;
    MAST::LevelSetNonlinearImplicitAssembly                  conduction_assembly;
    MAST::LevelSetNonlinearImplicitAssembly                  level_set_assembly;
    MAST::LevelSetEigenproblemAssembly                       eigen_assembly;
    MAST::StressAssembly                                     stress_assembly;
    MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
    MAST::HeatConductionNonlinearAssemblyElemOperations      conduction_elem_ops;
    MAST::StructuralModalEigenproblemAssemblyElemOperations  modal_elem_ops;
    
    // reinitialize the dof constraints before solution of the linear system
    // FIXME: we should be able to clear the constraint object from the
    // system before it goes out of scope, but libMesh::System does not
    // have a clear method. So, we are going to leave it as is, hoping
    // that libMesh::System will not attempt to use it (most likely, we
    // shoudl be ok).
    
    /////////////////////////////////////////////////////////////////////
    // first constrain the indicator function and solve
    /////////////////////////////////////////////////////////////////////
    SNESConvergedReason r;
    {
        libMesh::out << "Indicator Function" << std::endl;
        nonlinear_assembly.set_discipline_and_system(*_indicator_discipline, *_indicator_sys_init);
        conduction_elem_ops.set_discipline_and_system(*_indicator_discipline, *_indicator_sys_init);
        nonlinear_assembly.set_level_set_function(*_level_set_function);
        
        MAST::LevelSetConstrainDofs constrain(*_indicator_sys_init, *_level_set_function);
        constrain.constrain_all_negative_indices(true);
        _indicator_sys->attach_constraint_object(constrain);
        _indicator_sys->reinit_constraints();
        _indicator_sys->solve(conduction_elem_ops, nonlinear_assembly);
        r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
        (*_indicator_sys->nonlinear_solver).get_converged_reason();
        nonlinear_assembly.clear_discipline_and_system();
        nonlinear_assembly.clear_level_set_function();
        conduction_elem_ops.clear_discipline_and_system();
    }
    // if the solver diverged due to linear solve, then there is a problem with
    // this geometry and we need to return with a high value set for the
    // constraints
    if (r == SNES_DIVERGED_LINEAR_SOLVE) {
        
        obj = 1.e10;
        for (unsigned int i=0; i<_n_ineq; i++)
            fvals[i] = 1.e10;
        return;
    }
    
    
    /////////////////////////////////////////////////////////////////////
    // now, use the indicator function to constrain dofs in the structural
    // system
    /////////////////////////////////////////////////////////////////////
    MAST::MeshFieldFunction indicator(*_indicator_sys_init, "indicator");
    indicator.init(*_indicator_sys->solution);
    //MAST::IndicatorFunctionConstrainDofs constrain(*_sys_init, *_level_set_function, indicator);
    MAST::LevelSetConstrainDofs str_constrain (*_sys_init, *_level_set_function);
    MAST::LevelSetConstrainDofs cond_constrain(*_conduction_sys_init, *_level_set_function);
    _sys->attach_constraint_object(str_constrain);
    _sys->reinit_constraints();
    _sys->initialize_condensed_dofs(*_discipline);
    _conduction_sys->attach_constraint_object(cond_constrain);
    _conduction_sys->reinit_constraints();
    _conduction_sys->initialize_condensed_dofs(*_conduction_discipline);

    
    /////////////////////////////////////////////////////////////////////
    // first solve the conduction problem
    /////////////////////////////////////////////////////////////////////
    libMesh::out << "Conduction Solve" << std::endl;
    conduction_assembly.set_discipline_and_system(*_conduction_discipline, *_conduction_sys_init);
    conduction_assembly.set_level_set_function(*_level_set_function);
    conduction_assembly.set_level_set_velocity_function(*_level_set_vel);
    conduction_elem_ops.set_discipline_and_system(*_conduction_discipline, *_conduction_sys_init);
    _conduction_sys->solve(conduction_elem_ops, conduction_assembly);
    _temp_function->init(*_conduction_sys->solution);
    
    /////////////////////////////////////////////////////////////////////
    // first constrain the indicator function and solve
    /////////////////////////////////////////////////////////////////////
    nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    nonlinear_assembly.set_level_set_function(*_level_set_function);
    nonlinear_assembly.set_level_set_velocity_function(*_level_set_vel);
    //nonlinear_assembly.set_indicator_function(indicator);
    eigen_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    eigen_assembly.set_level_set_function(*_level_set_function);
    eigen_assembly.set_level_set_velocity_function(*_level_set_vel);
    stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
    level_set_assembly.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
    level_set_assembly.set_level_set_function(*_level_set_function);
    level_set_assembly.set_level_set_velocity_function(*_level_set_vel);
    nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    modal_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    //nonlinear_assembly.plot_sub_elems(true, false, true);
    
    
    
    
    MAST::LevelSetVolume                            volume(level_set_assembly.get_intersection());
    MAST::StressStrainOutputBase                    stress;
    volume.set_discipline_and_system(*_level_set_discipline, *_level_set_sys_init);
    stress.set_discipline_and_system(*_discipline, *_sys_init);
    volume.set_participating_elements_to_all();
    stress.set_participating_elements_to_all();
    stress.set_aggregation_coefficients(_p_val, _vm_rho, _stress_lim);
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the objective
    //////////////////////////////////////////////////////////////////////
    level_set_assembly.calculate_output(*_level_set_sys->solution, volume);
    obj       = volume.output_total();
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the stress constraint
    //////////////////////////////////////////////////////////////////////
    libMesh::out << "Static Solve" << std::endl;
    _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
    r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
    (*_sys->nonlinear_solver).get_converged_reason();
    
    // if the solver diverged due to linear solve, then there is a problem with
    // this geometry and we need to return with a high value set for the
    // constraints
    if (r == SNES_DIVERGED_LINEAR_SOLVE ||
        _sys->final_nonlinear_residual() > 1.e-1) {
        
        obj = 1.e10;
        for (unsigned int i=0; i<_n_ineq; i++)
            fvals[i] = 1.e10;
        return;
    }
    
    nonlinear_assembly.calculate_output(*_sys->solution, stress);
    fvals[0]  =  stress.output_total()/_stress_lim - 1.;  // g = sigma/sigma0-1 <= 0
    
    //stress_assembly.update_stress_strain_data(stress, *_sys->solution);
    //libMesh::ExodusII_IO(*_mesh).write_equation_systems("indicator.exo", *_eq_sys);
    //libMesh::ExodusII_IO(*_level_set_mesh).write_equation_systems("phi.exo", *_level_set_eq_sys);
    
    if (_n_eig_vals) {
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the eigenvalue constraint
        //////////////////////////////////////////////////////////////////////
        libMesh::out << "Eigen Solve" << std::endl;
        _sys->eigenproblem_solve(modal_elem_ops, eigen_assembly);
        Real eig_imag = 0.;
        //
        // hopefully, the solver found the requested number of eigenvalues.
        // if not, then we will set zero values for the ones it did not.
        //
        unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());
        std::vector<Real> eig(_n_eig_vals, 0.);
        
        // get the converged eigenvalues
        for (unsigned int i=0; i<n_conv; i++)      _sys->get_eigenvalue(0, eig[i], eig_imag);
        //
        //  eig > eig0
        //  -eig < -eig0
        //  -eig/eig0 < -1
        // -eig/eig0 + 1 < 0
        //
        for (unsigned int i=0; i<_n_eig_vals; i++)
            fvals[i+1] = -eig[i]/_ref_eig_val + 1.;
    }
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the objective sensitivities, if requested
    //////////////////////////////////////////////////////////////////////
    if (eval_obj_grad)
        _evaluate_volume_sensitivity(volume, level_set_assembly, obj_grad);
    
    //////////////////////////////////////////////////////////////////////
    // check to see if the sensitivity of constraint is requested
    //////////////////////////////////////////////////////////////////////
    bool if_grad_sens = false;
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_grad_sens = (if_grad_sens || eval_grads[i]);
    
    //////////////////////////////////////////////////////////////////////
    // evaluate the sensitivities for constraints
    //////////////////////////////////////////////////////////////////////
    if (if_grad_sens)
        _evaluate_constraint_sensitivity(stress,
                                         conduction_elem_ops,
                                         nonlinear_elem_ops,
                                         conduction_assembly,
                                         nonlinear_assembly,
                                         modal_elem_ops,
                                         eigen_assembly,
                                         eval_grads,
                                         grads);
    
    // also the stress data for plotting
    stress_assembly.update_stress_strain_data(stress, *_sys->solution);
    _temp_function->clear();
}




void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
_evaluate_constraint_sensitivity
(MAST::StressStrainOutputBase& stress,
 MAST::AssemblyElemOperations& conduction_elem_ops,
 MAST::AssemblyElemOperations& nonlinear_elem_ops,
 MAST::LevelSetNonlinearImplicitAssembly& conduction_assembly,
 MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
 MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
 MAST::LevelSetEigenproblemAssembly& eigen_assembly,
 const std::vector<bool>& eval_grads,
 std::vector<Real>& grads) {
    
    unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());

    //////////////////////////////////////////////////////////////////
    // first the structural adjoint solution
    //////////////////////////////////////////////////////////////////
    _sys->adjoint_solve(nonlinear_elem_ops, stress, nonlinear_assembly, false);

    
    //////////////////////////////////////////////////////////////////
    // next, use the structural adjoint to solve for the thermal adjoint
    //////////////////////////////////////////////////////////////////
    MAST::StressTemperatureAdjoint temperature_adj_elem_ops(stress);
    temperature_adj_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
    temperature_adj_elem_ops.set_thermal_assembly(conduction_assembly);
    temperature_adj_elem_ops.set_structural_solutions(*_sys->solution,
                                                      _sys->get_adjoint_solution());
    _conduction_sys->adjoint_solve(conduction_elem_ops,
                                   temperature_adj_elem_ops,
                                   conduction_assembly,
                                   false);

    
    std::unique_ptr<libMesh::NumericVector<Real>>
    dphi(_level_set_sys->solution->zero_clone().release());
    
    //////////////////////////////////////////////////////////////////
    // indices used by GCMMA follow this rule:
    // grad_k = dfi/dxj  ,  where k = j*NFunc + i
    //////////////////////////////////////////////////////////////////
    for (unsigned int i=0; i<_n_vars; i++) {
        
        dphi->zero();
        dphi->set(_dv_params[i].first, 1.);
        dphi->close();
        
        // initialize the level set perturbation function to create a velocity
        // field
        _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
        
        //////////////////////////////////////////////////////////////////////
        // stress sensitivity
        //////////////////////////////////////////////////////////////////////
        Real
        str_sens     = 0.,
        thermal_sens = 0.;
        
        str_sens =
        nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                _sys->get_adjoint_solution(),
                                                                *_dv_params[i].second,
                                                                nonlinear_elem_ops,
                                                                stress,
                                                                true);
        
        // do not include par_sigma/par_DV term since it was already included
        // in the structural contribution
        thermal_sens =
        conduction_assembly.calculate_output_adjoint_sensitivity(*_conduction_sys->solution,
                                                                 _conduction_sys->get_adjoint_solution(),
                                                                 *_dv_params[i].second,
                                                                 conduction_elem_ops,
                                                                 temperature_adj_elem_ops,
                                                                 false);
        grads[_n_ineq*i+0] = 1./_stress_lim*(str_sens + thermal_sens);
        stress.clear_sensitivity_data();
        
        //////////////////////////////////////////////////////////////////////
        // eigenvalue sensitivity, only if the values were requested
        //////////////////////////////////////////////////////////////////////
        if (_n_eig_vals) {
            
            std::vector<Real> sens;
            _sys->eigenproblem_sensitivity_solve(eigen_elem_ops,
                                                 eigen_assembly,
                                                 *_dv_params[i].second,
                                                 sens);
            for (unsigned int j=0; j<n_conv; j++)
                grads[_n_ineq*i+j+1] = -sens[j]/_ref_eig_val;
        }
    }
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
_init_system_and_discipline() {
    
    MAST::Examples::TopologyOptimizationLevelSet2D::_init_system_and_discipline();
    _conduction_sys         = &(_eq_sys->add_system<MAST::NonlinearSystem>("conduction"));
    
    _conduction_sys_init    = new MAST::HeatConductionSystemInitialization(*_conduction_sys,
                                                                           _conduction_sys->name(),
                                                                           _fetype);
    _conduction_discipline  = new MAST::PhysicsDisciplineBase(*_eq_sys);
}


void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::_init_dirichlet_conditions() {

    MAST::Examples::TopologyOptimizationLevelSet2D::_init_dirichlet_conditions();
    _init_conduction_boundary_dirichlet_constraint(1, "right_temperature_dirichlet");
    _init_conduction_boundary_dirichlet_constraint(3, "left_temperature_dirichlet");
    _conduction_discipline->init_system_dirichlet_bc(*_conduction_sys);
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
_init_conduction_boundary_dirichlet_constraint(const unsigned int bid,
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
    _conduction_discipline->add_dirichlet_bc(bid,  *dirichlet);
    
    
    this->register_loading(*dirichlet);
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
_init_loads() {
    
    MAST::Examples::TopologyOptimizationLevelSet2D::_init_loads();
    
    {
        Real
        length   = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
        frac     = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
        q_val    =  (*_input)(_prefix+"flux", "flux load on side of domain",  -2.e1);
        
        MAST::Examples::FluxLoad
        *q_f     = new MAST::Examples::FluxLoad("heat_flux", q_val, length, frac);;
        
        // initialize the load
        MAST::BoundaryConditionBase
        *q_load  = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
        q_load->add(*q_f);
        
        _conduction_discipline->add_side_load(2, *q_load);
        
        this->register_field_function(*q_f);
        this->register_loading(*q_load);
    }
    
    /*{
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
        
        _conduction_discipline->add_volume_load(0, *q_load);
        
        this->add_parameter(*s);
        this->add_parameter(*e);
        this->add_parameter(*T);
        this->add_parameter(*T0);
        this->register_field_function(*e_f);
        this->register_loading(*q_load);
    }*/
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::_init_temperature_load() {
    
    Real
    alpha_val   =  (*_input)(_prefix+"alpha_expansion", "coefficient of thermal expansion", 2.5e-5);
    
    MAST::Parameter
    *alpha           = new MAST::Parameter(       "alpha",  alpha_val),
    &zero            = this->get_parameter("zero");
    

    MAST::ConstantFieldFunction
    *alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *alpha),
    *ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature",   zero);
    
    _temp_function   = new MAST::Examples::TemperatureFunction(*_conduction_sys_init, "temperature");
    
    MAST::BoundaryConditionBase
    *T_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    
    T_load->add(*_temp_function);
    T_load->add(*ref_temp_f);
    _m_card->add(*alpha_f);
    
    _discipline->add_volume_load(0, *T_load);
    
    
    this->add_parameter(*alpha);
    this->register_field_function(*alpha_f);
    this->register_field_function(*_temp_function);
    this->register_field_function(*ref_temp_f);
    this->register_loading(*T_load);
}



void
MAST::Examples::ThermoelasticityTopologyOptimizationLevelSet2D::
_init_section_property() {
    
    MAST::Examples::TopologyOptimizationLevelSet2D::_init_section_property();
    _conduction_discipline->set_property_for_subdomain(0, _discipline->get_property_card(0));
}










