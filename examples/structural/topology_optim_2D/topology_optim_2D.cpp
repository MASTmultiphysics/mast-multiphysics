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
#include <structural/base/thermal_stress_jacobian_scaling_function.h>
#include "examples/structural/topology_optim_2D/topology_optim_2D.h"
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
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/structural_system_initialization.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "property_cards/material_property_card_base.h"
#include "property_cards/element_property_card_2D.h"
#include "optimization/optimization_interface.h"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/petsc_nonlinear_solver.h"


namespace MAST {
    
    class Phi:
    public MAST::FieldFunction<RealVectorX> {
        
    public:
        Phi(Real l1,
            Real l2,
            Real nx_mesh,
            Real ny_mesh,
            Real nx_holes,
            Real ny_holes):
        MAST::FieldFunction<RealVectorX>("Phi"),
        _l1  (l1),
        _l2  (l2),
        _nx_mesh  (nx_mesh),
        _ny_mesh  (ny_mesh),
        _nx_holes (nx_holes),
        _ny_holes (ny_holes),
        _pi  (acos(-1.)) {
         
            // initialize the locations at which the holes will be nucleated
            // first, along the x-axis
            if (_nx_holes == 1)
                _x_axis_hole_locations.insert(l1 * 0.5);
            else if (_nx_holes >= 2) {
                
                // add holes at the beginning and end
                _x_axis_hole_locations.insert(0.);
                _x_axis_hole_locations.insert(_l1);
                
                // now, add holes at uniformly spaced locations
                // in the domain
                Real
                dx = _l1/(1.*(_nx_holes-1));
                for (unsigned int i=2; i<_nx_holes; i++)
                    _x_axis_hole_locations.insert(dx*(i-1));
            }

            
            // now, along the y-axis
            if (_ny_holes == 1)
                _y_axis_hole_locations.insert(l2 * 0.5);
            else if (_ny_holes >= 2) {
                
                // add holes at the beginning and end
                _y_axis_hole_locations.insert(0.);
                _y_axis_hole_locations.insert(_l2);
                
                // now, add holes at uniformly spaced locations
                // in the domain
                Real
                dx = _l2/(1.*(_ny_holes-1));
                for (unsigned int i=2; i<_ny_holes; i++)
                    _y_axis_hole_locations.insert(dx*(i-1));
            }
        }
        virtual ~Phi() {}
        virtual void operator()(const libMesh::Point& p,
                                const Real t,
                                RealVectorX& v) const {
            
            libmesh_assert_less_equal(t, 1);
            libmesh_assert_equal_to(v.size(), 1);
            
            
            /*// circle
            //v(0) = -(pow(p(0)-_l1*.5, 2) + pow(p(1)-_l2*.5, 2) - pow(_r, 2));
            
            // waves
            Real
            c  = 0.5,
            pi = acos(-1.),
            x  = p(0)-.5*_l1,
            y  = p(1)-.5*_l2,
            r  = pow(pow(x,2)+pow(y,2),.5);
            //v(0) = 1.*cos(nx*r*_pi/_l1);
            v(0) = cos(2.*_nx*pi*x/_l1)+cos(2.*_ny*pi*y/_l2)+c;
            
            // linear
            //v(0) = (p(0)-_l1*0.5)*(-10.);
            //v(0) = (p(0)+p(1)-_l1)*(-10.);*/
            
            // the libMesh solution projection routine for Lagrange elements
            // will query the function value at the nodes. So, we figure
            // out which nodes should have zero values set to them.
            // if there is one hole in any direction, it will be in the
            // center of the domain. If there are more than 1, then two of
            // the holes will be on the boundary and others will fill the
            // interior evenly.
            
            const Real
            tol     = 1.e-6*std::min(_l1, _l2),
            dx_mesh = _l1/(1.*_nx_mesh),
            dy_mesh = _l2/(1.*_ny_mesh);
            
            std::set<Real>::const_iterator
            x_it_low = _x_axis_hole_locations.lower_bound(p(0)-dx_mesh),
            y_it_low = _y_axis_hole_locations.lower_bound(p(1)-dy_mesh);
            
            unsigned int
            n = 0;
            // see if the x-location needs a hole
            for ( ; x_it_low != _x_axis_hole_locations.end(); x_it_low++) {
                if (std::fabs(*x_it_low - p(0)) <= dx_mesh*0.5) {
                    n++;
                    break;
                }
            }
            
            // now check the y-location
            for ( ; y_it_low != _y_axis_hole_locations.end(); y_it_low++) {
                if (std::fabs(*y_it_low - p(1)) <= dy_mesh*0.5) {
                    n++;
                    break;
                }
            }

            if (n == 2)
                v(0) = -0.01;
            else
                v(0) = 0.01;
        }
    protected:
        Real
        _l1,
        _l2,
        _nx_mesh,
        _ny_mesh,
        _nx_holes,
        _ny_holes,
        _pi;
        std::set<Real> _x_axis_hole_locations;
        std::set<Real> _y_axis_hole_locations;
    };
    
    class Vel: public MAST::FieldFunction<Real> {
    public:
        Vel(): MAST::FieldFunction<Real>("vel") {}
        
        virtual void operator() (const libMesh::Point& p,
                                 const Real t,
                                 Real& v) const {
            
            // waves
            Real
            nt = 8.,
            th = atan2(p(0)-.15, p(1)-.15);
            v  = sin(nt*th/2.);
            
            // constant
            // v    = 1.;
        }
        
    protected:
        
    };
    
}


MAST::Examples::TopologyOptimizationLevelSet2D::
TopologyOptimizationLevelSet2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::StructuralExample2D  (comm_in),
MAST::FunctionEvaluation             (comm_in),
_obj_scaling                         (0.),
_stress_lim                          (0.),
_p_val                               (0.),
_vm_rho                              (0.),
_ref_eig_val                         (0.),
_n_eig_vals                          (0),
_level_set_mesh                      (nullptr),
_level_set_eq_sys                    (nullptr),
_level_set_sys                       (nullptr),
_level_set_sys_on_str_mesh           (nullptr),
_level_set_sys_init_on_str_mesh      (nullptr),
_indicator_sys                       (nullptr),
_level_set_sys_init                  (nullptr),
_indicator_sys_init                  (nullptr),
_indicator_discipline                (nullptr),
_level_set_discipline                (nullptr),
_level_set_function                  (nullptr),
_level_set_vel                       (nullptr),
_output                              (nullptr) {
    
}


MAST::Examples::TopologyOptimizationLevelSet2D::~TopologyOptimizationLevelSet2D() {

    if (!_initialized)
        return;
    
    delete _level_set_function;
    delete _level_set_vel;
    delete _level_set_sys_init;
    delete _indicator_sys_init;
    delete _indicator_discipline;
    delete _level_set_discipline;
    delete _level_set_eq_sys;
    delete _level_set_mesh;
    delete _output;
    delete _level_set_sys_init_on_str_mesh;
    
    for (unsigned int i=0; i<_dv_params.size(); i++)
        delete _dv_params[i].second;
}


void
MAST::Examples::TopologyOptimizationLevelSet2D::initialize_solution() {

    // initialize solution of the structural problem
    MAST::Examples::StructuralExample2D::initialize_solution();
    
    
    // initialize solution of the level set problem
    unsigned int
    nx_h    = (*_input)(_prefix+ "initial_level_set_n_holes_in_x",
                        "number of holes along x-direction for initial level-set field", 2.),
    ny_h    = (*_input)(_prefix+ "initial_level_set_n_holes_in_y",
                        "number of holes along y-direction for initial level-set field", 2.),
    nx_m    = (*_input)(_prefix+"level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_m    = (*_input)(_prefix+"level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);

    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3);

    Phi phi(length, height, nx_m, ny_m, nx_h, ny_h);
    _level_set_sys_init->initialize_solution(phi);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::init(MAST::Examples::GetPotWrapper& input,
                                                     const std::string& prefix) {
    
    libmesh_assert(_optimization_interface);
    
    // let all other data structures be initialized
    MAST::Examples::StructuralExample2D::init(input, prefix);

    // ask structure to use Mindlin bending operator
    dynamic_cast<MAST::ElementPropertyCard2D&>(*_p_card).set_bending_model(MAST::MINDLIN);
    
    /////////////////////////////////////////////////
    // now initialize the design data.
    /////////////////////////////////////////////////

    // first, initialize the level set functions over the domain
    this->initialize_solution();

    // next, define a new parameter to define design variable for nodal level-set
    // function value
    this->_init_phi_dvs();

    unsigned int
    max_inner_iters        = (*_input)(_prefix+"max_inner_iters", "maximum inner iterations in GCMMA", 15);
    
    Real
    constr_penalty         = (*_input)(_prefix+"constr_penalty", "constraint penalty in GCMMA", 50.),
    length                 = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height                 = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3);

    _optimization_interface->set_real_parameter   ( "constr_penalty",  constr_penalty);
    _optimization_interface->set_integer_parameter("max_inner_iters", max_inner_iters);
    
    _obj_scaling           = 100./length/height;
    _stress_lim            = (*_input)(_prefix+"vm_stress_limit", "limit von-mises stress value", 2.e8);
    _p_val                 = (*_input)(_prefix+"constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0);
    _vm_rho                = (*_input)(_prefix+"constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0);
    _level_set_vel         = new MAST::LevelSetBoundaryVelocity(2);
    _level_set_function    = new PhiMeshFunction;
    _output                = new libMesh::ExodusII_IO(*_mesh);

    _n_eig_vals            = (*_input)(_prefix+"n_eig", "number of eigenvalues to constrain", 5);
    if (_n_eig_vals) {
        // set only if the user requested eigenvalue constraints
        _ref_eig_val           = (*_input)(_prefix+"eigenvalue_low_bound", "lower bound enforced on eigenvalue constraints", 1.e3);
        _sys->set_n_requested_eigenvalues(_n_eig_vals);
    }

    // two inequality constraints: stress and eigenvalue.
    _n_ineq = 1+_n_eig_vals;

    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_optim_history.txt";
    this->set_output_file(output_name);
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::init_dvar(std::vector<Real>& x,
                                                          std::vector<Real>& xmin,
                                                          std::vector<Real>& xmax) {
    
    // one DV for each element
    x.resize(_n_vars);
    xmin.resize(_n_vars);
    xmax.resize(_n_vars);
    
    std::fill(xmin.begin(), xmin.end(),   -1.);
    std::fill(xmax.begin(), xmax.end(),    1.);
    
    // now, check if the user asked to initialize dvs from a previous file
    std::string
    nm    =  (*_input)(_prefix+"restart_optimization_file", "filename with optimization history for restart", "");
    
    if (nm.length()) {
        
        unsigned int
        iter = (*_input)(_prefix+"restart_optimization_iter", "restart iteration number from file", 0);
        this->initialize_dv_from_output_file(nm, iter, x);
    }
    else {
        
        for (unsigned int i=0; i<_n_vars; i++)
            x[i] = (*_dv_params[i].second)();
    }
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::evaluate(const std::vector<Real>& dvars,
                                                         Real& obj,
                                                         bool eval_obj_grad,
                                                         std::vector<Real>& obj_grad,
                                                         std::vector<Real>& fvals,
                                                         std::vector<bool>& eval_grads,
                                                         std::vector<Real>& grads) {
    
    libMesh::out << "New Evaluation" << std::endl;
    
    // copy DVs to level set function
    for (unsigned int i=0; i<_n_vars; i++)
        if (_dv_params[i].first >= _level_set_sys->solution->first_local_index() &&
            _dv_params[i].first <  _level_set_sys->solution->last_local_index())
        _level_set_sys->solution->set(_dv_params[i].first, dvars[i]);
    _level_set_sys->solution->close();
    _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
    _sys->solution->zero();
    
    /**********************************************************************
     * DO NOT zero out the gradient vector, since GCMMA needs it for the  *
     * subproblem solution                                                *
     **********************************************************************/
    MAST::LevelSetNonlinearImplicitAssembly                  nonlinear_assembly;
    MAST::LevelSetNonlinearImplicitAssembly                  level_set_assembly;
    MAST::LevelSetEigenproblemAssembly                       eigen_assembly;
    MAST::LevelSetStressAssembly                             stress_assembly;
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
        nonlinear_assembly.clear_level_set_function();
        nonlinear_assembly.clear_discipline_and_system();
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
    //MAST::LevelSetConstrainDofs constrain(*_sys_init, *_level_set_function);
    //_sys->attach_constraint_object(constrain);
    //_sys->reinit_constraints();
    //_sys->initialize_condensed_dofs(*_discipline);
    
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
    stress_assembly.init(*_level_set_function, nonlinear_assembly.get_dof_handler());
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
    obj       = volume.output_total() * _obj_scaling;

    //////////////////////////////////////////////////////////////////////
    // evaluate the stress constraint
    //////////////////////////////////////////////////////////////////////
    // tell the thermal jacobian scaling object about the assembly object
    MAST::Examples::ThermalJacobianScaling&
    scaling = dynamic_cast<MAST::Examples::ThermalJacobianScaling&>
            (this->get_field_function("thermal_jacobian_scaling"));
    scaling.clear_assembly();
    scaling.set_assembly(nonlinear_assembly);

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
                                         nonlinear_elem_ops,
                                         nonlinear_assembly,
                                         modal_elem_ops,
                                         eigen_assembly,
                                         eval_grads,
                                         grads);

    // also the stress data for plotting
    stress_assembly.update_stress_strain_data(stress, *_sys->solution);
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::
_evaluate_volume_sensitivity(MAST::LevelSetVolume& volume,
                             MAST::LevelSetNonlinearImplicitAssembly& assembly,
                             std::vector<Real>& obj_grad) {
    
    // iterate over each DV, create a sensitivity vector and calculate the
    // volume sensitivity explicitly
    std::unique_ptr<libMesh::NumericVector<Real>>
    dphi(_level_set_sys->solution->zero_clone().release());
    
    for (unsigned int i=0; i<_n_vars; i++) {
        
        dphi->zero();
        // set the value only if the dof corresponds to a local node
        if (_dv_params[i].first >=  dphi->first_local_index() &&
            _dv_params[i].first <   dphi->last_local_index())
            dphi->set(_dv_params[i].first, 1.);
        dphi->close();

        _level_set_vel->init(*_level_set_sys_init, *_level_set_sys->solution, *dphi);
        
        assembly.calculate_output_direct_sensitivity(*_level_set_sys->solution,
                                                     *dphi,
                                                     *_dv_params[i].second,
                                                     volume);
        obj_grad[i] = _obj_scaling * volume.output_sensitivity_total(*_dv_params[i].second);
    }
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::
_evaluate_constraint_sensitivity
(MAST::StressStrainOutputBase& stress,
 MAST::AssemblyElemOperations& nonlinear_elem_ops,
 MAST::LevelSetNonlinearImplicitAssembly& nonlinear_assembly,
 MAST::StructuralModalEigenproblemAssemblyElemOperations& eigen_elem_ops,
 MAST::LevelSetEigenproblemAssembly& eigen_assembly,
 const std::vector<bool>& eval_grads,
 std::vector<Real>& grads) {

    unsigned int n_conv = std::min(_n_eig_vals, _sys->get_n_converged_eigenvalues());

    _sys->adjoint_solve(nonlinear_elem_ops, stress, nonlinear_assembly, false);
    
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
        grads[_n_ineq*i+0] = 1./_stress_lim*
        nonlinear_assembly.calculate_output_adjoint_sensitivity(*_sys->solution,
                                                                _sys->get_adjoint_solution(),
                                                                *_dv_params[i].second,
                                                                nonlinear_elem_ops,
                                                                stress);
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
MAST::Examples::TopologyOptimizationLevelSet2D::_init_mesh() {
    
    // first call the parent method
    MAST::Examples::StructuralExample2D::_init_mesh();
    
    _level_set_mesh = new libMesh::SerialMesh(MAST::Examples::StructuralExample2D::comm());
    
    // identify the element type from the input file or from the order
    // of the element
    
    unsigned int
    nx_divs = (*_input)(_prefix+"level_set_nx_divs", "number of elements of level-set mesh along x-axis", 10),
    ny_divs = (*_input)(_prefix+"level_set_ny_divs", "number of elements of level-set mesh along y-axis", 10);
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3);
    
    std::string
    t = (*_input)(_prefix+"level_set_elem_type", "type of geometric element in the level set mesh", "quad4");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::QUAD4)
        e_type = libMesh::QUAD9;
    else if (_fetype.order > 1 && e_type == libMesh::TRI3)
        e_type = libMesh::TRI6;
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_level_set_mesh,
                                                 nx_divs, ny_divs,
                                                 0, length,
                                                 0, height,
                                                 e_type);
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_eq_sys() {
    
    // first call the parent method
    MAST::Examples::StructuralExample2D::_init_eq_sys();
    _level_set_eq_sys->init();
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_system_and_discipline() {
    
    // first initialize the structural system and discipline
    MAST::Examples::StructuralExample2D::_init_system_and_discipline();

    // FEType to initialize the system
    // get the order and type of element
    std::string
    order_str   = (*_input)(_prefix+ "level_set_fe_order", "order of finite element shape basis functions for level set method",    "first");
    
    libMesh::Order
    o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
    _level_set_fetype = libMesh::FEType(o, libMesh::LAGRANGE);
    

    // now initialize the level set related data structures
    _level_set_eq_sys      = new libMesh::EquationSystems(*_level_set_mesh);
    _level_set_sys         = &(_level_set_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
    _level_set_sys_init    = new MAST::LevelSetSystemInitialization(*_level_set_sys,
                                                                    _level_set_sys->name(),
                                                                    _level_set_fetype);
    _level_set_discipline  = new MAST::LevelSetDiscipline(*_eq_sys);
    
    _level_set_sys_on_str_mesh      = &(_eq_sys->add_system<MAST::NonlinearSystem>("level_set"));
    _indicator_sys                  = &(_eq_sys->add_system<MAST::NonlinearSystem>("indicator"));
    _level_set_sys_init_on_str_mesh = new MAST::LevelSetSystemInitialization(*_level_set_sys_on_str_mesh,
                                                                             _level_set_sys->name(),
                                                                             _level_set_fetype);
    _indicator_sys_init             = new MAST::HeatConductionSystemInitialization(*_indicator_sys,
                                                                                   _indicator_sys->name(),
                                                                                   _fetype);
    _indicator_discipline           = new MAST::PhysicsDisciplineBase(*_eq_sys);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_dirichlet_conditions() {
    
    // constrain only the left and right boundaries
    this->_init_boundary_dirichlet_constraint(1, "right_constraint");
    this->_init_boundary_dirichlet_constraint(3, "left_constraint");
    
    ///////////////////////////////////////////////////////////////////////
    // initialize Dirichlet conditions for structural system
    ///////////////////////////////////////////////////////////////////////
    _discipline->init_system_dirichlet_bc(*_sys);

    _init_indicator_system_dirichlet_conditions();
}




void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_indicator_system_dirichlet_conditions() {
    
    ///////////////////////////////////////////////////////////////////////
    // initialize Dirichlet conditions for indicator system
    ///////////////////////////////////////////////////////////////////////
    MAST::DirichletBoundaryCondition
    *dirichlet  = new MAST::DirichletBoundaryCondition;   // right boundary
    dirichlet->init(1, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(1,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    dirichlet   = new MAST::DirichletBoundaryCondition;   // left boundary
    dirichlet->init(3, _indicator_sys_init->vars());
    _indicator_discipline->add_dirichlet_bc(3,  *dirichlet);
    this->register_loading(*dirichlet);                   // register, so that it will be deleted later
    _indicator_discipline->init_system_dirichlet_bc(*_indicator_sys);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_loads() {
    
    Real
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.2),
    p_val   =  (*_input)(_prefix+"pressure", "pressure on side of domain",   2.e4);
    
    MAST::Examples::FluxLoad
    *press_f         = new MAST::Examples::FluxLoad( "pressure", p_val, length, frac),
    *flux_f          = new MAST::Examples::FluxLoad("heat_flux", -2.e6, length, frac);
    
    // initialize the load
    MAST::BoundaryConditionBase
    *p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE),
    *f_load          = new MAST::BoundaryConditionBase(MAST::HEAT_FLUX);
    
    p_load->add(*press_f);
    _discipline->add_side_load(2, *p_load);
    
    f_load->add(*flux_f);
    _indicator_discipline->add_side_load(2, *f_load);
    
    this->register_field_function(*press_f);
    this->register_field_function(*flux_f);
    this->register_loading(*p_load);
    this->register_loading(*f_load);

    _init_temperature_load();
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_material() {
    
    MAST::Examples::StructuralExample2D::_init_material();

    // add the conduction property to this card
    Real
    kval      = (*_input)(_prefix+"k", "thermal conductivity",  1.e-2),
    cpval     = (*_input)(_prefix+"cp", "thermal capacitance",  864.);
    
    
    MAST::Parameter
    *k         = new MAST::Parameter("k",          kval),
    *cp        = new MAST::Parameter("cp",        cpval);
    
    MAST::ConstantFieldFunction
    *k_f     = new MAST::ConstantFieldFunction( "k_th",      *k),
    *cp_f    = new MAST::ConstantFieldFunction(   "cp",     *cp);
    
    this->add_parameter(*k);
    this->add_parameter(*cp);
    this->register_field_function(*k_f);
    this->register_field_function(*cp_f);
    
    _m_card->add(*k_f);
    _m_card->add(*cp_f);
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_section_property() {
    
    MAST::Examples::StructuralExample2D::_init_section_property();
    _indicator_discipline->set_property_for_subdomain(0, _discipline->get_property_card(0));
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::_init_phi_dvs() {

    libmesh_assert(_initialized);
    // this assumes that level set is defined using lagrange shape functions
    libmesh_assert_equal_to(_level_set_fetype.family, libMesh::LAGRANGE);
    
    Real
    tol     = 1.e-6,
    length  = (*_input)(_prefix+"length", "length of domain along x-axis", 0.3),
    height  = (*_input)(_prefix+"height", "length of domain along y-axis", 0.3),
    frac    = (*_input)(_prefix+"load_length_fraction", "fraction of boundary length on which pressure will act", 0.2);

    unsigned int
    dof_id  = 0;
    
    Real
    val     = 0.;

    // all ranks will have DVs defined for all variables. So, we should be
    // operating on a replicated mesh
    libmesh_assert(_level_set_mesh->is_replicated());
    
    std::vector<Real> local_phi(_level_set_sys->solution->size());
    _level_set_sys->solution->localize(local_phi);
    
    // iterate over all the node values
    libMesh::MeshBase::const_node_iterator
    it  = _level_set_mesh->nodes_begin(),
    end = _level_set_mesh->nodes_end();
    
    // maximum number of dvs is the number of nodes on the level set function
    // mesh. We will evaluate the actual number of dvs
    _dv_params.reserve(_level_set_mesh->n_nodes());
    _n_vars = 0;
    
    for ( ; it!=end; it++) {
        
        const libMesh::Node& n = **it;
        
        dof_id                     = n.dof_number(0, 0, 0);

        // only if node is not on the upper edge
        if ((std::fabs(n(1)-height) > tol) ||
            (n(0) > length*.5*(1.+frac))   ||
            (n(0) < length*.5*(1.-frac))) {
     
            std::ostringstream oss;
            oss << "dv_" << _n_vars;
            val  = local_phi[dof_id];
            
//            // on the traction free boundary, set everything to be zero, so that there
//            // is always a boundary there that the optimizer can move
//            if (n(1) < tol                     ||
//                std::fabs(n(1) - height) < tol) {
//
//                if (dof_id >= _level_set_sys->solution->first_local_index() &&
//                    dof_id <  _level_set_sys->solution->last_local_index())
//                    _level_set_sys->solution->set(dof_id, 0.);
//                val = 0.;
//            }

            _dv_params.push_back(std::pair<unsigned int, MAST::Parameter*>());
            _dv_params[_n_vars].first  = dof_id;
            _dv_params[_n_vars].second = new MAST::Parameter(oss.str(), val);
            _dv_params[_n_vars].second->set_as_topology_parameter(true);
            
            _n_vars++;
        }
        else {
            // set value at the material points to a small positive number
            if (dof_id >= _level_set_sys->solution->first_local_index() &&
                dof_id <  _level_set_sys->solution->last_local_index())
                _level_set_sys->solution->set(dof_id, 0.01);
        }
    }
    
    _level_set_sys->solution->close();
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::level_set_solve() {
    
    libmesh_assert(_initialized);
    
    bool
    output      = (*_input)(_prefix+"if_output", "if write output to a file", false),
    propagate   = (*_input)(_prefix+"if_propagate", "if propagate level set, or reinitialize it", true);
    
    
    std::string
    output_name = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output");
    output_name += "_level_set.exo";
    
    // create the nonlinear assembly object
    std::unique_ptr<MAST::TransientAssembly> level_set_assembly;
    if (propagate)
        level_set_assembly.reset(new MAST::TransientAssembly);
    else {
        MAST::LevelSetReinitializationTransientAssembly
        *assembly = new MAST::LevelSetReinitializationTransientAssembly;
        
        libMesh::NumericVector<Real>
        &base_sol = _level_set_sys->add_vector("base_sol");
        base_sol  = *_level_set_sys->solution;
        
        assembly->set_reference_solution(base_sol);
        _level_set_discipline->set_level_set_propagation_mode(false);
        level_set_assembly.reset(assembly);
    }
    MAST::LevelSetTransientAssemblyElemOperations            level_set_elem_ops;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  level_set_solver;
    
    // now solve the system
    level_set_assembly->set_discipline_and_system(*_level_set_discipline,
                                                  *_level_set_sys_init);
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
    // time solver parameters
    unsigned int
    t_step                         = 0,
    n_steps                        = (*_input)(_prefix+"level_set_n_transient_steps", "number of transient time-steps", 100);
    level_set_solver.dt            = (*_input)(_prefix+"level_set_dt", "time-step size",    1.e-3);
    level_set_solver.beta          = 0.5;
    
    // set the previous state to be same as the current state to account for
    // zero velocity as the initial condition
    level_set_solver.solution(1).zero();
    level_set_solver.solution(1).add(1., level_set_solver.solution());
    level_set_solver.solution(1).close();
    
    
    // loop over time steps
    while (t_step <= n_steps) {
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << _level_set_sys->time
        << std::endl;
        
        // write the time-step
        if (output) {
            
            exodus_writer.write_timestep(output_name,
                                         *_eq_sys,
                                         t_step+1,
                                         _level_set_sys->time);
        }
        
        level_set_solver.solve(*level_set_assembly);
        
        level_set_solver.advance_time_step();
        t_step++;
    }
    
    level_set_assembly->clear_discipline_and_system();
}



void
MAST::Examples::TopologyOptimizationLevelSet2D::output(unsigned int iter,
                                                       const std::vector<Real>& x,
                                                       Real obj,
                                                       const std::vector<Real>& fval,
                                                       bool if_write_to_optim_file) {

    libmesh_assert_equal_to(x.size(), _n_vars);
    
    std::string
    output_name  = (*_input)(_prefix+"output_file_root", "prefix of output file names", "output"),
    modes_name   = output_name + "modes.exo";
    output_name += "_optim.exo";

    // copy DVs to level set function
    for (unsigned int i=0; i<_n_vars; i++)
        if (_dv_params[i].first >= _level_set_sys->solution->first_local_index() &&
            _dv_params[i].first <  _level_set_sys->solution->last_local_index())
            _level_set_sys->solution->set(_dv_params[i].first, x[i]);
    _level_set_sys->solution->close();
    _level_set_function->init(*_level_set_sys_init, *_level_set_sys->solution);
    _level_set_sys_init_on_str_mesh->initialize_solution(_level_set_function->get_mesh_function());

    std::vector<bool> eval_grads(this->n_ineq(), false);
    std::vector<Real> f(this->n_ineq(), 0.), grads;
    this->evaluate(x, obj, false, grads, f, eval_grads, grads);
    
    _output->write_timestep(output_name, *_eq_sys, iter+1, (1.*iter));

    if (_n_eig_vals) {
        
        //////////////////////////////////////////////////////////////////////////
        // eigenvalue analysis: write modes to file
        //////////////////////////////////////////////////////////////////////////
        libMesh::ExodusII_IO writer(*_mesh);
        Real eig_r, eig_i;
        for (unsigned int i=0; i<_sys->get_n_converged_eigenvalues(); i++) {
            _sys->get_eigenpair(i, eig_r, eig_i, *_sys->solution);
            writer.write_timestep(modes_name, *_eq_sys, i+1, i);
        }
        _sys->solution->zero();
    }
    
    MAST::FunctionEvaluation::output(iter, x, obj/_obj_scaling, fval, if_write_to_optim_file);
}

