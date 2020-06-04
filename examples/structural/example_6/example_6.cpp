/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2020  Manav Bhatia and MAST authors
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
#include <iomanip>

// MAST includes
#include "examples/base/input_wrapper.h"
#include "examples/fluid/meshing/cylinder.h"
#include "level_set/filter_base.h"
#include "level_set/level_set_parameter.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/ks_stress_output.h"
#include "elasticity/smooth_ramp_stress_output.h"
#include "elasticity/level_set_stress_assembly.h"
#include "elasticity/compliance_output.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/mesh_field_function.h"
#include "base/nonlinear_system.h"
#include "base/transient_assembly.h"
#include "base/boundary_condition_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "solver/slepc_eigen_solver.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/isotropic_element_property_card_3D.h"
#include "optimization/gcmma_optimization_interface.h"
#include "optimization/npsol_optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "examples/structural/base/bracket_2d_model.h"
#include "examples/structural/base/bracket_3d_model.h"
#include "examples/structural/base/inplane_2d_model.h"
#include "examples/structural/base/truss_2d_model.h"
#include "examples/structural/base/eyebar_2d_model.h"


// libMesh includes
#include "libmesh/fe_type.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/parallel.h"


void
_optim_obj(int*    mode,
           int*    n,
           double* x,
           double* f,
           double* g,
           int*    nstate);
void
_optim_con(int*    mode,
           int*    ncnln,
           int*    n,
           int*    ldJ,
           int*    needc,
           double* x,
           double* c,
           double* cJac,
           int*    nstate);

//
// BEGIN_TRANSLATE SIMP topology optimization
//
//   \tableofcontents
//
//  This example computes the optimal topology of a structure subject to
//  specified boundary conditions (Dirichlet and Neumann). An element-wise density
//  is used to parameterize the topology.
//
//  Elasticity function with the penalty term
class ElasticityFunction:
public MAST::FieldFunction<Real> {
public:
    ElasticityFunction(Real E0, Real rho_min, Real penalty,
                       MAST::MeshFieldFunction& rho):
    MAST::FieldFunction<Real>("E"),
    _E0(E0),
    _rho_min(rho_min),
    _penalty(penalty),
    _rho(rho) { }
    virtual ~ElasticityFunction(){}
    void set_penalty_val(Real penalty) {_penalty = penalty;}
    
    virtual bool depends_on(const MAST::FunctionBase& f) const { return true;}
    virtual void operator() (const libMesh::Point& p, const Real t, Real& v) const {

        RealVectorX v1;
        _rho(p, t, v1);
        
        v = _E0 * (_rho_min + (1.-_rho_min) * pow(v1(0), _penalty));
    }

    virtual void derivative(const MAST::FunctionBase& f,
                            const libMesh::Point& p, const Real t, Real& v) const {
        
        RealVectorX v1, dv1;
        _rho(p, t, v1);
        _rho.derivative(f, p, t, dv1);
        
        v = _E0 * (1.-_rho_min) * _penalty * pow(v1(0), _penalty-1.) * dv1(0);
    }

protected:
    Real                    _E0; // value of the material Young's modulus
    Real                    _rho_min; // lower limit on density
    Real                    _penalty; // value of penalty term
    MAST::MeshFieldFunction &_rho;
};



class ElementParameterDependence:
public MAST::AssemblyBase::ElemParameterDependence {
public:
    ElementParameterDependence(const MAST::FilterBase& filter):
    MAST::AssemblyBase::ElemParameterDependence(true), _filter(filter) {}
    virtual ~ElementParameterDependence() {}
    virtual bool if_elem_depends_on_parameter(const libMesh::Elem& e,
                                              const MAST::FunctionBase& p) const {
        const MAST::LevelSetParameter
        &p_ls = dynamic_cast<const MAST::LevelSetParameter&>(p);
        
        return _filter.if_elem_in_domain_of_influence(e, *p_ls.level_set_node());
    }

private:
    const MAST::FilterBase& _filter;
};


template <typename T>
class TopologyOptimizationSIMP:
public MAST::FunctionEvaluation {
    
public:
    
    bool                                      _initialized;
    MAST::Examples::GetPotWrapper&            _input;
    
    std::string                               _problem;
    Real                                      _volume;
    Real                                      _obj_scaling;
    Real                                      _stress_penalty;
    Real                                      _perimeter_penalty;
    Real                                      _stress_lim;
    Real                                      _p_val, _vm_rho;
    Real                                      _vf;      // volume fraction
    Real                                      _rho_min; // lower limit on density

    ElasticityFunction*                       _Ef;
    libMesh::UnstructuredMesh*                _mesh;
    
    libMesh::EquationSystems*                 _eq_sys;
    
    MAST::NonlinearSystem*                    _sys;
    libMesh::System*                          _density_sys;
    
    MAST::StructuralSystemInitialization*     _sys_init;
    
    MAST::PhysicsDisciplineBase*              _discipline;

    MAST::StructuralNearNullVectorSpace*      _nsp;

    MAST::FilterBase*                         _filter;
    
    MAST::MaterialPropertyCardBase*           _m_card;
    MAST::ElementPropertyCardBase*            _p_card;
    
    MAST::MeshFieldFunction*                  _density_function;
    libMesh::ExodusII_IO*                     _output;
    
    libMesh::FEType                           _fetype;
    libMesh::FEType                           _density_fetype;
    
    std::vector<MAST::Parameter*>             _params_for_sensitivity;
    std::map<std::string, MAST::Parameter*>   _parameters;
    std::set<MAST::FunctionBase*>             _field_functions;
    std::set<MAST::BoundaryConditionBase*>    _boundary_conditions;
    std::set<unsigned int>                    _dv_dof_ids;
    std::set<unsigned int>                    _dirichlet_bc_ids;

    std::vector<std::pair<unsigned int, MAST::Parameter*>>  _dv_params;

    
    //
    //  \section ex_6_system_discipline  System and Discipline
    //
    void _init_system_and_discipline() {
        
        //
        // make sure that the mesh has been initialized
        //
        libmesh_assert(_mesh);
        
        //
        // create the equation system
        //
        _eq_sys    = new  libMesh::EquationSystems(*_mesh);
        
        //
        // create the libmesh system and set the preferences for structural
        // eigenvalue problems
        //
        _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
        _sys->set_eigenproblem_type(libMesh::GHEP);
        
        //
        // initialize the system to the right set of variables
        //
        _sys_init       = new MAST::StructuralSystemInitialization(*_sys,
                                                                   _sys->name(),
                                                                   _fetype);
        _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
        
        //
        // Initialize the system for level set function.
        // A level set function is defined on a coarser mesh than the structural
        // mesh.
        // A level set function is assumed to be a first-order Lagrange finite element
        //
        _density_fetype      = libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE);
        _density_sys         = &(_eq_sys->add_system<libMesh::ExplicitSystem>("density"));
        _density_sys->add_variable("rho", _density_fetype);
    }

    
    void _init_eq_sys() {
        
        _eq_sys->init();
        _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
        _sys->set_exchange_A_and_B(true);
    }
    

    //
    //   variables added to the mesh
    //
    void _init_fetype() {
        
        // FEType to initialize the system. Get the order and type of element.
        std::string
        order_str   = _input("fe_order", "order of finite element shape basis functions",    "first"),
        family_str  = _input("fe_family",      "family of finite element shape functions", "lagrange");
        
        libMesh::Order
        o  = libMesh::Utility::string_to_enum<libMesh::Order>(order_str);
        libMesh::FEFamily
        fe = libMesh::Utility::string_to_enum<libMesh::FEFamily>(family_str);
        _fetype = libMesh::FEType(o, fe);
    }
        
    
    //
    //   \section ex_6_properties Properties
    //
    //
    //
    //   \subsection ex_6_material_properties Material Properties
    //

    void _init_material() {
        
        _rho_min  = _input("rho_min", "lower limit on density variable", 1.e-8);

        Real
        Eval      = _input("E", "modulus of elasticity", 72.e9),
        penalty   = _input("rho_penalty", "SIMP modulus of elasticity penalty", 4.),
        rhoval    = _input("rho", "material density", 2700.),
        nu_val    = _input("nu", "Poisson's ratio",  0.33),
        alpha_val = _input("alpha", "coefficient of thermal expansion", 1.5e-5),
        kval      = _input("k", "thermal conductivity",  1.e-2),
        cpval     = _input("cp", "thermal capacitance",  864.);
        
        
        MAST::Parameter
        *rho       = new MAST::Parameter("rho",      rhoval),
        *nu        = new MAST::Parameter("nu",       nu_val),
        *k         = new MAST::Parameter("k",          kval),
        *cp        = new MAST::Parameter("cp",        cpval),
        *alpha     = new MAST::Parameter("alpha_expansion", alpha_val);

        MAST::ConstantFieldFunction
        *rho_f   = new MAST::ConstantFieldFunction(  "rho",                 *rho),
        *nu_f    = new MAST::ConstantFieldFunction(   "nu",                  *nu),
        *k_f     = new MAST::ConstantFieldFunction( "k_th",                   *k),
        *cp_f    = new MAST::ConstantFieldFunction(   "cp",                  *cp),
        *alpha_f = new MAST::ConstantFieldFunction("alpha_expansion",     *alpha);

        _Ef      = new ElasticityFunction(Eval, _rho_min, penalty, *_density_function);
        
        _parameters[  rho->name()]     = rho;
        _parameters[   nu->name()]     = nu;
        _parameters[    k->name()]     = k;
        _parameters[   cp->name()]     = cp;
        _parameters[alpha->name()]     = alpha;
        _field_functions.insert(_Ef);
        _field_functions.insert(rho_f);
        _field_functions.insert(nu_f);
        _field_functions.insert(k_f);
        _field_functions.insert(cp_f);
        _field_functions.insert(alpha_f);

        _m_card = new MAST::IsotropicMaterialPropertyCard;
        _m_card->add(*_Ef);
        _m_card->add(*rho_f);
        _m_card->add(*nu_f);
        _m_card->add(*k_f);
        _m_card->add(*cp_f);
        _m_card->add(*alpha_f);
    }

    
    //
    //   \subsection  ex_6_section_properties Section Properties
    //

    void _init_section_property(){
        
        
        
        Real
        kappa_val = _input("kappa", "shear correction factor",  5./6.),
        th_v      =  _input("th", "thickness of 2D element",  0.001);
        
        MAST::Parameter
        *th       = new MAST::Parameter("th", th_v),
        *kappa    = new MAST::Parameter("kappa", kappa_val),
        *zero     = new MAST::Parameter("zero", 0.);
        
        MAST::ConstantFieldFunction
        *th_f     = new MAST::ConstantFieldFunction("h",       *th),
        *kappa_f  = new MAST::ConstantFieldFunction("kappa",  *kappa),
        *hoff_f   = new MAST::ConstantFieldFunction("off",   *zero);
        
        
        _parameters[th->name()]    = th;
        _parameters[kappa->name()] = kappa;
        _parameters[zero->name()]  = zero;
        _field_functions.insert(th_f);
        _field_functions.insert(kappa_f);
        _field_functions.insert(hoff_f);
        
        typename T::SectionPropertyCardType
        *p_card   = new typename T::SectionPropertyCardType;
        
        _p_card   = p_card;
        
        // set nonlinear strain if requested
        bool
        nonlinear = _input("if_nonlinear", "flag to turn on/off nonlinear strain", false);
        if (nonlinear) p_card->set_strain(MAST::NONLINEAR_STRAIN);
        
        p_card->add(*th_f);
        p_card->add(*kappa_f);
        p_card->add(*hoff_f);
        p_card->set_material(*_m_card);
        _discipline->set_property_for_subdomain(0, *p_card);
    }
    

    //
    //   \section  ex_6_initial_solution Initial Density field
    //
    //
    
    
    void initialize_solution() {
        
        //
        // initialize density field to a constant value of the specified
        // volume fraction
        //
        _vf    = _input("volume_fraction", "upper limit for the voluem fraction", 0.5);

        _density_sys->solution->zero();
        _density_sys->solution->add(_vf);
        _density_sys->solution->close();
    }
    
        
    //
    //   \subsection  ex_6_design_variable_init   Design Variables
    //
    //   initializes the design variable vector, called by the
    //   optimization interface.
    //
    void init_dvar(std::vector<Real>& x,
                   std::vector<Real>& xmin,
                   std::vector<Real>& xmax) {
        
        x.resize(_n_vars);
        xmin.resize(_n_vars);
        xmax.resize(_n_vars);
        
        std::fill(xmin.begin(), xmin.end(),    0.);
        std::fill(xmax.begin(), xmax.end(),    1.e0);

        //
        // now, check if the user asked to initialize dvs from a previous file
        //
        std::string
        nm    =  _input("restart_optimization_file", "filename with optimization history for restart", "");
        
        if (nm.length()) {
            
            unsigned int
            iter = _input("restart_optimization_iter", "restart iteration number from file", 0);
            this->initialize_dv_from_output_file(nm, iter, x);
        }
        else {
            
            for (unsigned int i=0; i<_n_vars; i++)
                x[i] = (*_dv_params[i].second)();
        }
    }


    //
    //  \subsection  ex_6_function_evaluation Function Evaluation
    //
    void evaluate(const std::vector<Real>& dvars,
                  Real& obj,
                  bool eval_obj_grad,
                  std::vector<Real>& obj_grad,
                  std::vector<Real>& fvals,
                  std::vector<bool>& eval_grads,
                  std::vector<Real>& grads) {
        
        libMesh::out << "New Evaluation" << std::endl;
        
        // copy DVs to level set function
        libMesh::NumericVector<Real>
        &base_phi = _density_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, dvars[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_density_sys->solution);

        // this will create a localized vector in _level_set_sys->curret_local_solution
        _density_sys->update();
        
        // create a localized vector for use in interpolation
        std::unique_ptr<libMesh::NumericVector<Real>>
        local_density_sol(libMesh::NumericVector<Real>::build(_sys->comm()).release());
        local_density_sol->init(_density_sys->n_dofs(),
                                _density_sys->n_local_dofs(),
                                _density_sys->get_dof_map().get_send_list(),
                                false,
                                libMesh::GHOSTED);
        _density_sys->solution->localize(*local_density_sol);

        _density_function->clear();
        _density_function->init(*local_density_sol, true);
        _sys->solution->zero();
        
        //////////////////////////////////////////////////////////////////////
        // check to see if the sensitivity of constraint is requested
        //////////////////////////////////////////////////////////////////////
        bool if_grad_sens = false;
        for (unsigned int i=0; i<eval_grads.size(); i++)
            if_grad_sens = (if_grad_sens || eval_grads[i]);

        // if sensitivity analysis is requested, then initialize the vectors
        std::vector<libMesh::NumericVector<Real>*> sens_vecs;
        if (eval_obj_grad || if_grad_sens)
            _initialize_sensitivity_data(sens_vecs);
        
        //*********************************************************************
        // DO NOT zero out the gradient vector, since GCMMA needs it for the  *
        // subproblem solution                                                *
        //*********************************************************************
        MAST::NonlinearImplicitAssembly                          nonlinear_assembly;
        MAST::StressAssembly                                     stress_assembly;
        MAST::StructuralNonlinearAssemblyElemOperations          nonlinear_elem_ops;
        
        /////////////////////////////////////////////////////////////////////
        // first constrain the indicator function and solve
        /////////////////////////////////////////////////////////////////////
        nonlinear_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        stress_assembly.set_discipline_and_system(*_discipline, *_sys_init);
        nonlinear_elem_ops.set_discipline_and_system(*_discipline, *_sys_init);
        
        MAST::StressStrainOutputBase                    stress;
        MAST::ComplianceOutput                          compliance;
        stress.set_discipline_and_system(*_discipline, *_sys_init);
        stress.set_participating_elements_to_all();
        stress.set_aggregation_coefficients(_p_val, 1., _vm_rho, _stress_lim) ;
        compliance.set_participating_elements_to_all();
        compliance.set_discipline_and_system(*_discipline, *_sys_init);

        //////////////////////////////////////////////////////////////////////
        // evaluate the stress constraint
        //////////////////////////////////////////////////////////////////////
        
        libMesh::out << "Static Solve" << std::endl;

        Real
        penalty          = _input("rho_penalty", "SIMP modulus of elasticity penalty", 4.),
        stress_penalty   = _input("stress_rho_penalty", "SIMP modulus of elasticity penalty for stress evaluation", 0.5);
        // set the elasticity penalty for solution
        _Ef->set_penalty_val(penalty);
        
        _sys->solve(nonlinear_elem_ops, nonlinear_assembly);
        SNESConvergedReason
        r = dynamic_cast<libMesh::PetscNonlinearSolver<Real>&>
        (*_sys->nonlinear_solver).get_converged_reason();
        
        // if the solver diverged due to linear solve, then there is a problem with
        // this geometry and we need to return with a high value set for the
        // constraints
        if (r == SNES_DIVERGED_LINEAR_SOLVE ||
            _sys->final_nonlinear_residual() > 1.e-1) {
            
            obj = 1.e11;
            for (unsigned int i=0; i<_n_ineq; i++)
                fvals[i] = 1.e11;
            return;
        }
        
        Real
        max_vm = 0.,
        vm_agg = 0.,
        vol    = 0.,
        comp   = 0.;

        // ask the system to update so that the localized solution is available for
        // further processing
        _sys->update();

        //////////////////////////////////////////////////////////////////////
        // evaluate the functions
        //////////////////////////////////////////////////////////////////////
        
        // evaluate the volume for used in the problem setup
        _evaluate_volume(*local_density_sol, sens_vecs, &vol, nullptr);
        libMesh::out << "volume: " << vol << std::endl;

        // evaluate the output based on specified problem type
        if (_problem == "compliance_volume") {
            
            nonlinear_assembly.calculate_output(*_sys->current_local_solution, false, compliance);
            comp      = compliance.output_total();
            obj       = _obj_scaling * comp;
            fvals[0]  = vol/_volume - _vf; // vol/vol0 - a <=
            libMesh::out << "compliance: " << comp << std::endl;
        }
        else if (_problem == "volume_stress") {
            
            // set the elasticity penalty for stress evaluation
            _Ef->set_penalty_val(stress_penalty);
            nonlinear_assembly.calculate_output(*_sys->current_local_solution, false, stress);
            max_vm    = stress.get_maximum_von_mises_stress();
            vm_agg    = stress.output_total();
            obj       = _obj_scaling * vol;
            fvals[0]  =  stress.output_total()/_stress_lim - 1.;  // g = sigma/sigma0-1 <= 0
            libMesh::out
            << "  max: "    << max_vm
            << "  constr: " << vm_agg
            << std::endl;
        }
        else
            libmesh_error();
        

        //////////////////////////////////////////////////////////////////////
        // evaluate the objective sensitivities, if requested
        //////////////////////////////////////////////////////////////////////
        if (eval_obj_grad) {
            
            if (_problem == "compliance_volume") {
                
                _evaluate_compliance_sensitivity(compliance,
                                                 nonlinear_elem_ops,
                                                 nonlinear_assembly,
                                                 obj_grad);
                
                for (unsigned int i=0; i<obj_grad.size(); i++)
                    obj_grad[i] *= _obj_scaling;
            }
            else if (_problem == "volume_stress") {
                
                _evaluate_volume(*local_density_sol, sens_vecs, nullptr, &obj_grad);
                for (unsigned int i=0; i<grads.size(); i++)
                    obj_grad[i] *= _obj_scaling;
            }
            else
                libmesh_error();
        }
        
        
        //////////////////////////////////////////////////////////////////////
        // evaluate the sensitivities for constraints
        //////////////////////////////////////////////////////////////////////
        if (if_grad_sens) {
            
            if (_problem == "compliance_volume") {
                
                _evaluate_volume(*local_density_sol, sens_vecs, nullptr, &grads);
                for (unsigned int i=0; i<grads.size(); i++)
                    grads[i] /= _volume;
            }
            else if (_problem == "volume_stress") {
                
                _evaluate_stress_sensitivity(penalty,
                                             stress_penalty,
                                             stress,
                                             nonlinear_elem_ops,
                                             nonlinear_assembly,
                                             grads);
            }
            else
                libmesh_error();
        }
        
        //
        // also the stress data for plotting
        //
        _Ef->set_penalty_val(stress_penalty);
        stress_assembly.update_stress_strain_data(stress, *_sys->solution);
        
        _density_function->clear();
        _clear_sensitivity_data(sens_vecs);
    }

    
    //
    // \subsection ex_6_sensitivity_vectors Initialize sensitivity data
    //
    void _initialize_sensitivity_data(std::vector<libMesh::NumericVector<Real>*>& dphi_vecs) {

        libmesh_assert_equal_to(dphi_vecs.size(), 0);
        
        dphi_vecs.resize(_n_vars, nullptr);
        
        // localized vectors are used for the level
        // set mesh function since it uses a different mesh than the analysis mesh
        // and the two can have different partitionings in the paralle environment.
        for (unsigned int i=0; i<_n_vars; i++) {

            libMesh::NumericVector<Real>
            *vec = nullptr;

            // non-zero value of the DV perturbation
            std::map<unsigned int, Real> nonzero_val;
            nonzero_val[_dv_params[i].first] = 1.;
            
            vec = libMesh::NumericVector<Real>::build(_sys->comm()).release();
            vec->init(_density_sys->n_dofs(),
                      _density_sys->n_local_dofs(),
                      _density_sys->get_dof_map().get_send_list(),
                      false,
                      libMesh::GHOSTED);
            _filter->compute_filtered_values(nonzero_val, *vec, false);

            dphi_vecs[i] = vec;
        }

        for ( unsigned int i=0; i<_n_vars; i++)
            dphi_vecs[i]->close();

        // we will use this ghosted vector to initialize the mesh function,
        // which is setup to reuse this vector, so we have to store it
        for ( unsigned int i=0; i<_n_vars; i++)
            _density_function->init_sens(*_dv_params[i].second, *dphi_vecs[i], true);
    }

    
    void _clear_sensitivity_data(std::vector<libMesh::NumericVector<Real>*>& dphi_vecs) {

        // delete the vectors that we do not need any more
        for (unsigned int i=0; i<dphi_vecs.size(); i++)
            delete dphi_vecs[i];
        dphi_vecs.clear();

        _density_function->clear();
    }
    

    //
    //  \subsection  ex_6_volume_sensitivity Sensitivity of Material Volume
    //
    void _evaluate_volume(const libMesh::NumericVector<Real>& sol,
                          const std::vector<libMesh::NumericVector<Real>*>& dsol_vecs,
                          Real               *volume,
                          std::vector<Real>  *grad) {
        
        unsigned int
        sys_num = _density_sys->number();
        
        if (volume) {

            *volume = 0.;
            
            libMesh::MeshBase::element_iterator
            it    =  _mesh->active_local_elements_begin(),
            end   =  _mesh->active_local_elements_end();
            
            Real
            rho = 0.;
            
            for ( ; it != end; it++) {
                
                const libMesh::Elem& e = **it;
                
                // compute the average element density value
                rho = 0.;
                for (unsigned int i=0; i<e.n_nodes(); i++) {
                    const libMesh::Node& n = *e.node_ptr(i);
                    rho += sol.el(n.dof_number(sys_num, 0, 0));
                }
                rho /= (1. * e.n_nodes());

                // use this density value to compute the volume
                *volume  +=  e.volume() * rho;
            }
            
            this->comm().sum(*volume);
        }
        
        
        if (grad) {
            
            std::fill(grad->begin(), grad->end(), 0.);
            ElementParameterDependence dep(*_filter);
            
            for (unsigned int i=0; i<_n_vars; i++) {
                
                libMesh::MeshBase::element_iterator
                it    =  _mesh->active_local_elements_begin(),
                end   =  _mesh->active_local_elements_end();
                
                Real
                drho = 0.;
                
                for ( ; it != end; it++) {
                    
                    const libMesh::Elem& e = **it;
                    
                    // do not compute if the element is not in the domain
                    // of influence of the parameter
                    if (!dep.if_elem_depends_on_parameter(e, *_dv_params[i].second))
                        continue;
                    
                    // compute the average element density value
                    drho = 0.;
                    for (unsigned int j=0; j<e.n_nodes(); j++) {
                        const libMesh::Node& n = *e.node_ptr(j);
                        drho += dsol_vecs[i]->el(n.dof_number(sys_num, 0, 0));
                    }
                    drho /= (1. * e.n_nodes());
                    
                    // use this density value to compute the volume
                    (*grad)[i]  +=  e.volume() * drho;
                }
            }

            this->comm().sum(*grad);
        }
    }
    
    
    
    //
    //  \subsection  ex_6_stress_sensitivity Sensitivity of Stress and Eigenvalues
    //
    void
    _evaluate_stress_sensitivity
    (const Real                    penalty,
     const Real                    stress_penalty,
     MAST::StressStrainOutputBase& stress,
     MAST::AssemblyElemOperations& nonlinear_elem_ops,
     MAST::NonlinearImplicitAssembly& nonlinear_assembly,
     std::vector<Real>& grads) {
        
        _sys->adjoint_solve(*_sys->current_local_solution,
                            false,
                            nonlinear_elem_ops,
                            stress,
                            nonlinear_assembly,
                            false);
        
        ElementParameterDependence dep(*_filter);
        nonlinear_assembly.attach_elem_parameter_dependence_object(dep);
        
        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        // first compute the sensitivity contribution from dot product of adjoint vector
        // and residual sensitivity
        std::vector<Real>
        g1(_n_vars, 0.),
        g2(_n_vars, 0.);
        std::vector<const MAST::FunctionBase*>
        p_vec(_n_vars, nullptr);
        for (unsigned int i=0; i<_n_vars; i++)
            p_vec[i] = _dv_params[i].second;

        //////////////////////////////////////////////////////////////////////
        // stress sensitivity
        //////////////////////////////////////////////////////////////////////
        _Ef->set_penalty_val(penalty);
        nonlinear_assembly.calculate_output_adjoint_sensitivity_multiple_parameters_no_direct
        (*_sys->current_local_solution,
         false,
         _sys->get_adjoint_solution(),
         p_vec,
         nonlinear_elem_ops,
         stress,
         g1);


        // we will skip the summation of sensitivity inside the stress object to minimize
        // communication cost. Instead, we will do it at the end for the constraint vector
        stress.set_skip_comm_sum(true);
        _Ef->set_penalty_val(stress_penalty);
        for (unsigned int i=0; i<_n_vars; i++) {
            
            // set the elasticity penalty for solution, which is needed for
            // computation of the residual sensitivity
            nonlinear_assembly.calculate_output_direct_sensitivity(*_sys->current_local_solution,
                                                                   false,
                                                                   nullptr,
                                                                   false,
                                                                   *_dv_params[i].second,
                                                                   stress);
            g2[i] = stress.output_sensitivity_total(*_dv_params[i].second);

            stress.clear_sensitivity_data();
        }
        stress.set_skip_comm_sum(false);

        // now sum the values across processors to sum the partial derivatives for
        // each parameter
        _sys->comm().sum(g2);
        
        // now compute contribution to the stress constraint
        for (unsigned int i=0; i<_n_vars; i++)
            grads[1*i+0] = 1./_stress_lim * (g1[i] + g2[i]);

        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    
    void
    _evaluate_compliance_sensitivity
    (MAST::ComplianceOutput&                  compliance,
     MAST::AssemblyElemOperations&            nonlinear_elem_ops,
     MAST::NonlinearImplicitAssembly&         nonlinear_assembly,
     std::vector<Real>& grads) {
        
        // Adjoint solution for compliance = - X
        
        ElementParameterDependence dep(*_filter);
        nonlinear_assembly.attach_elem_parameter_dependence_object(dep);

        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        // first compute the sensitivity contribution from dot product of adjoint vector
        // and residual sensitivity
        std::vector<Real>
        g1(_n_vars, 0.),
        g2(_n_vars, 0.);
        std::vector<const MAST::FunctionBase*>
        p_vec(_n_vars, nullptr);
        for (unsigned int i=0; i<_n_vars; i++)
            p_vec[i] = _dv_params[i].second;
        
        //////////////////////////////////////////////////////////////////////
        // compliance sensitivity
        //////////////////////////////////////////////////////////////////////
        // set the elasticity penalty for solution, which is needed for
        // computation of the residual sensitivity
        nonlinear_assembly.calculate_output_adjoint_sensitivity_multiple_parameters_no_direct
        (*_sys->current_local_solution,
         false,
         *_sys->current_local_solution,
         p_vec,
         nonlinear_elem_ops,
         compliance,
         g1);

        compliance.set_skip_comm_sum(true);
        for (unsigned int i=0; i<_n_vars; i++) {
            
            nonlinear_assembly.calculate_output_direct_sensitivity(*_sys->current_local_solution,
                                                                   false,
                                                                   nullptr,
                                                                   false,
                                                                   *_dv_params[i].second,
                                                                   compliance);
            g2[i] = compliance.output_sensitivity_total(*_dv_params[i].second);
        }
        
        compliance.set_skip_comm_sum(false);
        
        // now sum the values across processors to sum the partial derivatives for
        // each parameter
        _sys->comm().sum(g2);

        for (unsigned int i=0; i<_n_vars; i++)
            grads[i] = -g1[i] + g2[i];

        nonlinear_assembly.clear_elem_parameter_dependence_object();
    }

    void set_n_vars(const unsigned int n_vars) {_n_vars = n_vars;}

    //
    //  \subsection  ex_6_design_output  Output of Design Iterate
    //
    void output(unsigned int iter,
                const std::vector<Real>& x,
                Real obj,
                const std::vector<Real>& fval,
                bool if_write_to_optim_file) {
        
        libmesh_assert_equal_to(x.size(), _n_vars);
        
        Real
        sys_time     = _sys->time;
        
        std::string
        output_name  = _input("output_file_root", "prefix of output file names", "output"),
        modes_name   = output_name + "modes.exo";
        
        std::ostringstream oss;
        oss << "output_optim.e-s." << std::setfill('0') << std::setw(5) << iter ;
        
        //
        // copy DVs to level set function
        //
        libMesh::NumericVector<Real>
        &base_phi = _density_sys->get_vector("base_values");
        
        for (unsigned int i=0; i<_n_vars; i++)
            if (_dv_params[i].first >= base_phi.first_local_index() &&
                _dv_params[i].first <  base_phi.last_local_index())
                base_phi.set(_dv_params[i].first, x[i]);
        base_phi.close();
        _filter->compute_filtered_values(base_phi, *_density_sys->solution);
        // create a ghosted vector for use in interpolation
        std::unique_ptr<libMesh::NumericVector<Real>>
        local_density_sol(libMesh::NumericVector<Real>::build(_sys->comm()).release());
        local_density_sol->init(_density_sys->n_dofs(),
                                _density_sys->n_local_dofs(),
                                _density_sys->get_dof_map().get_send_list(),
                                false,
                                libMesh::GHOSTED);
        _density_sys->solution->localize(*local_density_sol);

        _density_function->init(*local_density_sol, true);
        
        std::vector<bool> eval_grads(this->n_ineq(), false);
        std::vector<Real> f(this->n_ineq(), 0.), grads;
        this->evaluate(x, obj, false, grads, f, eval_grads, grads);
        
        _sys->time = iter;
        _sys_init->get_stress_sys().time = iter;
        // "1" is the number of time-steps in the file, as opposed to the time-step number.
        libMesh::ExodusII_IO(*_mesh).write_timestep(oss.str(), *_eq_sys, 1, (1.*iter));
        
        _density_function->clear();
        
        //
        // set the value of time back to its original value
        //
        _sys->time = sys_time;
        
        //
        // increment the parameter values
        //
        unsigned int
        update_freq = _input("update_freq_optim_params", "number of iterations after which the optimization parameters are updated", 50),
        factor = iter/update_freq ;
        if (factor > 0 && iter%update_freq == 0) {
            
            Real
            p_val           = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0),
            vm_rho          = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0),
            constr_penalty  = _input("constraint_penalty", "constraint penalty in GCMMA",      50.),
            max_penalty     = _input("max_constraint_penalty", "maximum constraint penalty in GCMMA",      1.e7),
            initial_step    = _input("initial_rel_step", "initial relative step length in GCMMA",      0.5),
            min_step        = _input("minimum_rel_step", "minimum relative step length in GCMMA",      0.001);
            
            constr_penalty = std::min(constr_penalty*pow(10, factor), max_penalty);
            initial_step   = std::max(initial_step-0.01*factor, min_step);
            _p_val         = std::min(p_val+2*factor, 10.);
            _vm_rho        = std::min(vm_rho+factor*0.5, 2.);
            libMesh::out
            << "Updated values: c = " << constr_penalty
            << "  step = " << initial_step
            << "  p = " << _p_val
            << "  rho = " << _vm_rho << std::endl;
            
            _optimization_interface->set_real_parameter   ( "constraint_penalty",   constr_penalty);
            _optimization_interface->set_real_parameter   ("initial_rel_step",        initial_step);
        }

        MAST::FunctionEvaluation::output(iter, x, obj/_obj_scaling, f, if_write_to_optim_file);
    }

#if MAST_ENABLE_SNOPT == 1
    MAST::FunctionEvaluation::funobj
    get_objective_evaluation_function() {
    
        return _optim_obj;
    }

    MAST::FunctionEvaluation::funcon
    get_constraint_evaluation_function() {
    
        return _optim_con;
    }
#endif
    
    
    //
    // \section  ex_6_initialization  Initialization
    //
    //   \subsection  ex_6_constructor  Constructor
    //
    
    TopologyOptimizationSIMP(const libMesh::Parallel::Communicator& comm_in,
                                   MAST::Examples::GetPotWrapper& input):
    MAST::FunctionEvaluation             (comm_in),
    _initialized                         (false),
    _input                               (input),
    _problem                             (""),
    _volume                              (0.),
    _obj_scaling                         (0.),
    _stress_penalty                      (0.),
    _perimeter_penalty                   (0.),
    _stress_lim                          (0.),
    _p_val                               (0.),
    _vm_rho                              (0.),
    _vf                                  (0.),
    _rho_min                             (0.),
    _mesh                                (nullptr),
    _eq_sys                              (nullptr),
    _sys                                 (nullptr),
    _sys_init                            (nullptr),
    _discipline                          (nullptr),
    _nsp                                 (nullptr),
    _filter                              (nullptr),
    _m_card                              (nullptr),
    _p_card                              (nullptr),
    _output                              (nullptr) {
        
        libmesh_assert(!_initialized);
        
        //
        // call the initialization routines for each component
        //
        _mesh    =  new libMesh::SerialMesh(this->comm());

        _init_fetype();
        T::init_analysis_mesh(*this, *_mesh);
        _init_system_and_discipline();
        T::init_analysis_dirichlet_conditions(*this);
        _init_eq_sys();

        _nsp  = new MAST::StructuralNearNullVectorSpace;
        _sys->nonlinear_solver->nearnullspace_object = _nsp;

        // density function is used by elasticity modulus function. So, we
        // initialize this here
        _density_function        = new MAST::MeshFieldFunction(*_density_sys, "rho", libMesh::GHOSTED);

        _init_material();
        T::init_structural_loads(*this);
        T::init_thermoelastic_loads(*this);
        _init_section_property();
        _initialized = true;
        
        /////////////////////////////////////////////////
        // now initialize the design data.
        /////////////////////////////////////////////////
        
        //
        // initialize density field to a constant value of the specified
        // volume fraction
        //
        _vf    = _input("volume_fraction", "upper limit for the voluem fraction", 0.5);

        _density_sys->solution->zero();
        _density_sys->solution->add(_vf);
        _density_sys->solution->close();

        //
        // next, define a new parameter to define design variable for nodal level-set
        // function value
        //
        T::init_simp_dvs(*this);
        
        Real
        filter_radius          = _input("filter_radius", "radius of geometric filter for level set field", 0.015);
        _filter                = new MAST::FilterBase(*_density_sys, filter_radius, _dv_dof_ids);
        libMesh::NumericVector<Real>& vec = _density_sys->add_vector("base_values");
        vec = *_density_sys->solution;
        vec.close();


        _problem               = _input("problem_type", "{compliance_volume, volume_stress}", "compliance_volume");
        _volume                = T::reference_volume(*this);
        _obj_scaling           = 1./_volume;
        _stress_penalty        = _input("stress_penalty", "penalty value for stress_constraint", 0.);
        _perimeter_penalty     = _input("perimeter_penalty", "penalty value for perimeter in the objective function", 0.);
        _stress_lim            = _input("vm_stress_limit", "limit von-mises stress value", 2.e8);
        _p_val                 = _input("constraint_aggregation_p_val", "value of p in p-norm stress aggregation", 2.0);
        _vm_rho                = _input("constraint_aggregation_rho_val", "value of rho in p-norm stress aggregation", 2.0);
        _output                = new libMesh::ExodusII_IO(*_mesh);
        
        //
        // two inequality constraints: stress and eigenvalue.
        //
        _n_ineq = 1;
        
        std::string
        output_name = _input("output_file_root", "prefix of output file names", "output");
        output_name += "_optim_history.txt";
        this->set_output_file(output_name);
        
    }
    
    //
    //   \subsection  ex_6_destructor  Destructor
    //
    ~TopologyOptimizationSIMP() {
        
        {
            std::set<MAST::BoundaryConditionBase*>::iterator
            it   = _boundary_conditions.begin(),
            end  = _boundary_conditions.end();
            for ( ; it!=end; it++)
                delete *it;
        }
        
        {
            std::set<MAST::FunctionBase*>::iterator
            it   = _field_functions.begin(),
            end  = _field_functions.end();
            for ( ; it!=end; it++)
                delete *it;
        }
        
        {
            std::map<std::string, MAST::Parameter*>::iterator
            it   = _parameters.begin(),
            end  = _parameters.end();
            for ( ; it!=end; it++)
                delete it->second;
        }
        
        if (!_initialized)
            return;

        delete _nsp;
        
        delete _m_card;
        delete _p_card;
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _sys_init;
        
        delete _filter;
        delete _output;
        delete _density_function;

        for (unsigned int i=0; i<_dv_params.size(); i++)
            delete _dv_params[i].second;
    }
    

};


//
//   \subsection  ex_6_wrappers_snopt  Wrappers for SNOPT
//

MAST::FunctionEvaluation* _my_func_eval = nullptr;

#if MAST_ENABLE_SNOPT == 1

unsigned int
it_num = 0;

void
_optim_obj(int*    mode,
           int*    n,
           double* x,
           double* f,
           double* g,
           int*    nstate) {

    //
    // make sure that the global variable has been setup
    //
    libmesh_assert(_my_func_eval);

    //
    // initialize the local variables
    //
    Real
    obj = 0.;

    unsigned int
    n_vars  =  _my_func_eval->n_vars(),
    n_con   =  _my_func_eval->n_eq()+_my_func_eval->n_ineq();

    libmesh_assert_equal_to(*n, n_vars);

    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (0);

    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), false);
    
    //
    // copy the dvars
    //
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];


    _my_func_eval->_evaluate_wrapper(dvars,
                                     obj,
                                     *mode>0,       // request the derivatives of obj
                                     obj_grad,
                                     fvals,
                                     eval_grads,
                                     grads);

    //
    // now copy them back as necessary
    //
    *f  = obj;
    if (*mode > 0) {
        
        // output data to the file
        _my_func_eval->_output_wrapper(it_num, dvars, obj, fvals, true);
        it_num++;

        for (unsigned int i=0; i<n_vars; i++)
            g[i] = obj_grad[i];
    }

    if (obj > 1.e10) *mode = -1;
}






void
_optim_con(int*    mode,
           int*    ncnln,
           int*    n,
           int*    ldJ,
           int*    needc,
           double* x,
           double* c,
           double* cJac,
           int*    nstate) {

    //
    // make sure that the global variable has been setup
    //
    libmesh_assert(_my_func_eval);

    //
    // initialize the local variables
    //
    Real
    obj = 0.;

    unsigned int
    n_vars  =  _my_func_eval->n_vars(),
    n_con   =  _my_func_eval->n_eq()+_my_func_eval->n_ineq();

    libmesh_assert_equal_to(    *n, n_vars);
    libmesh_assert_equal_to(*ncnln, n_con);

    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (n_vars*n_con, 0.);

    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), *mode>0);

    //
    // copy the dvars
    //
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];


    _my_func_eval->_evaluate_wrapper(dvars,
                                     obj,
                                     false,       // request the derivatives of obj
                                     obj_grad,
                                     fvals,
                                     eval_grads,
                                     grads);

    //
    // now copy them back as necessary
    //
    // first the constraint functions
    //
    for (unsigned int i=0; i<n_con; i++)
        c[i] = fvals[i];

    if (*mode > 0) {
        //
        // next, the constraint gradients
        //
        for (unsigned int i=0; i<n_con*n_vars; i++)
            cJac[i] = grads[i];
    }
    
    if (obj > 1.e10) *mode = -1;
}
#endif

//
//   \subsection  ex_6_main Main function
//

int main(int argc, char* argv[]) {

    libMesh::LibMeshInit init(argc, argv);

    MAST::Examples::GetPotWrapper
    input(argc, argv, "input");

    std::unique_ptr<MAST::FunctionEvaluation>
    top_opt;
    
    std::string
    mesh = input("mesh", "inplane2d, bracket2d, truss2d, eyebar2d", "inplane2d");
    
    if (mesh == "inplane2d") {
        top_opt.reset
        (new TopologyOptimizationSIMP<MAST::Examples::Inplane2DModel>
         (init.comm(), input));
    }
    else if (mesh == "bracket2d") {
        top_opt.reset
        (new TopologyOptimizationSIMP<MAST::Examples::Bracket2DModel>
         (init.comm(), input));
    }
    else if (mesh == "eyebar2d") {
        top_opt.reset
        (new TopologyOptimizationSIMP<MAST::Examples::Eyebar2DModel>
         (init.comm(), input));
    }
    else if (mesh == "truss2d") {
        top_opt.reset
        (new TopologyOptimizationSIMP<MAST::Examples::Truss2DModel>
         (init.comm(), input));
    }
    else if (mesh == "bracket3d") {
        top_opt.reset
        (new TopologyOptimizationSIMP<MAST::Examples::Bracket3DModel>
         (init.comm(), input));
    }
    else
        libmesh_error();
    
    _my_func_eval = top_opt.get();
    
    //MAST::NLOptOptimizationInterface optimizer(NLOPT_LD_SLSQP);
    std::unique_ptr<MAST::OptimizationInterface>
    optimizer;
    
    std::string
    s          = input("optimizer", "optimizer to use in the example", "gcmma");

    if (s == "gcmma") {

        optimizer.reset(new MAST::GCMMAOptimizationInterface);
        
        unsigned int
        max_inner_iters        = input("max_inner_iters", "maximum inner iterations in GCMMA", 15);
        
        Real
        constr_penalty         = input("constraint_penalty", "constraint penalty in GCMMA", 50.),
        initial_rel_step       = input("initial_rel_step", "initial step size in GCMMA", 1.e-2),
        asymptote_reduction    = input("asymptote_reduction", "reduction of aymptote in GCMMA", 0.7),
        asymptote_expansion    = input("asymptote_expansion", "expansion of asymptote in GCMMA", 1.2);
        
        optimizer->set_real_parameter   ("constraint_penalty",  constr_penalty);
        optimizer->set_real_parameter   ("initial_rel_step",  initial_rel_step);
        optimizer->set_real_parameter   ("asymptote_reduction",  asymptote_reduction);
        optimizer->set_real_parameter   ("asymptote_expansion",  asymptote_expansion);
        optimizer->set_integer_parameter(   "max_inner_iters", max_inner_iters);
    }
    else if (s == "snopt") {
        
        optimizer.reset(new MAST::NPSOLOptimizationInterface);
    }
    else {
        
        libMesh::out
        << "Unrecognized optimizer specified: " << s << std::endl;
        libmesh_error();
    }
    
    if (optimizer.get()) {
        
        optimizer->attach_function_evaluation_object(*top_opt);

        bool
        verify_grads = input("verify_gradients", "If true, the gradients of objective and constraints will be verified without optimization", false);
        if (verify_grads) {
            
            std::vector<Real> xx1(top_opt->n_vars()), xx2(top_opt->n_vars());
            top_opt->init_dvar(xx1, xx2, xx2);
            top_opt->verify_gradients(xx1);
        }
        else
            optimizer->optimize();
    }
    
    // END_TRANSLATE
    return 0;
}
