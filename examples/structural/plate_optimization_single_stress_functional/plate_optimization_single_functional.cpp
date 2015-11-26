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

// C++ includes
#include <iostream>

// MAST includes
#include "examples/structural/plate_optimization_single_stress_functional/plate_optimization_single_functional.h"
#include "driver/driver_base.h"
#include "elasticity/stress_output_base.h"
#include "optimization/optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "elasticity/structural_nonlinear_assembly.h"

// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"


extern
MAST::FunctionEvaluation *__my_func_eval;



void
plate_stress_functional_optim_obj(int*    mode,
                                  int*    n,
                                  double* x,
                                  double* f,
                                  double* g,
                                  int*    nstate) {
    
    
    // make sure that the global variable has been setup
    libmesh_assert(__my_func_eval);
    
    // initialize the local variables
    Real
    obj = 0.;
    
    unsigned int
    n_vars  =  __my_func_eval->n_vars(),
    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
    
    libmesh_assert_equal_to(*n, n_vars);
    
    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (0);
    
    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), false);
    
    // copy the dvars
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];
    
    
    __my_func_eval->evaluate(dvars,
                             obj,
                             true,       // request the derivatives of obj
                             obj_grad,
                             fvals,
                             eval_grads,
                             grads);
    
    
    // now copy them back as necessary
    *f  = obj;
    for (unsigned int i=0; i<n_vars; i++)
        g[i] = obj_grad[i];
}






void
plate_stress_functional_optim_con(int*    mode,
                                  int*    ncnln,
                                  int*    n,
                                  int*    ldJ,
                                  int*    needc,
                                  double* x,
                                  double* c,
                                  double* cJac,
                                  int*    nstate) {
    
    
    // make sure that the global variable has been setup
    libmesh_assert(__my_func_eval);
    
    // initialize the local variables
    Real
    obj = 0.;
    
    unsigned int
    n_vars  =  __my_func_eval->n_vars(),
    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
    
    libmesh_assert_equal_to(    *n, n_vars);
    libmesh_assert_equal_to(*ncnln, n_con);
    
    std::vector<Real>
    dvars   (*n,    0.),
    obj_grad(*n,    0.),
    fvals   (n_con, 0.),
    grads   (n_vars*n_con, 0.);
    
    std::vector<bool>
    eval_grads(n_con);
    std::fill(eval_grads.begin(), eval_grads.end(), true);
    
    // copy the dvars
    for (unsigned int i=0; i<n_vars; i++)
        dvars[i] = x[i];
    
    
    __my_func_eval->evaluate(dvars,
                             obj,
                             true,       // request the derivatives of obj
                             obj_grad,
                             fvals,
                             eval_grads,
                             grads);
    
    
    // now copy them back as necessary
    
    // first the constraint functions
    for (unsigned int i=0; i<n_con; i++)
        c[i] = fvals[i];
    
    // next, the constraint gradients
    for (unsigned int i=0; i<n_con*n_vars; i++)
        cJac[i] = grads[i];
    
    
}





MAST::PlateBendingSingleStressFunctionalSizingOptimization::
PlateBendingSingleStressFunctionalSizingOptimization(std::ostream& output):
MAST::FunctionEvaluation(output),
_initialized(false),
_n_divs_x(0),
_n_divs_y(0),
_n_elems(0),
_n_stations_x(0) { }



void
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
init(GetPot& infile,
     libMesh::ElemType e_type,
     bool if_vk) {

    libmesh_assert(!_initialized);

    // number of elements
    _n_divs_x    = infile("n_divs_x", 16);
    _n_divs_y    = infile("n_divs_y", 16);
    _n_elems     = _n_divs_x*_n_divs_y;
    if (e_type == libMesh::TRI3)
        _n_elems    *= 2;

    // number of stations
    _n_stations_x = infile("n_stations", 8);
    
    
    // now setup the optimization data
    _n_vars                = _n_stations_x; // for thickness variable
    _n_eq                  = 0;
    _n_ineq                = 1;            // all stress constraints combined into single functional
    _max_iters             = 1000;
    
    
    
    // length of domain
    _length        = infile("length", 0.50);
    _width         = infile("width",  0.25);
    
    // limit stress
    _stress_limit  = infile("max_stress", 4.00e8);
    
    // create the mesh
    _mesh          = new libMesh::SerialMesh(_init->comm());
    
    // initialize the mesh with one element
    libMesh::MeshTools::Generation::build_square(*_mesh,
                                                 _n_divs_x, _n_divs_y,
                                                 0, _length,
                                                 0, _width,
                                                 e_type);
    _mesh->prepare_for_use();
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<libMesh::NonlinearImplicitSystem>("structural"));
    
    // FEType to initialize the system
    libMesh::FEType fetype (libMesh::FIRST, libMesh::LAGRANGE);
    
    // initialize the system to the right set of variables
    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
                                                               _sys->name(),
                                                               fetype);
    _discipline     = new MAST::StructuralDiscipline(*_eq_sys);
    
    
    // create and add the boundary condition and loads
    std::vector<unsigned int> constrained_vars(4);
    // not constraning ty, tz will keep it simply supported
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 5;  // tz
    
    // create and add the boundary condition and loads
    _dirichlet_bottom = new MAST::DirichletBoundaryCondition;
    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
    _dirichlet_top    = new MAST::DirichletBoundaryCondition;
    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
    
    _dirichlet_bottom->init (0, constrained_vars);
    _dirichlet_right->init  (1, constrained_vars);
    _dirichlet_top->init    (2, constrained_vars);
    _dirichlet_left->init   (3, constrained_vars);
    
    _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,    *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
    
    _discipline->init_system_dirichlet_bc(dynamic_cast<libMesh::System&>(*_sys));
    
    // initialize the equation system
    _eq_sys->init();
    
    
    // initialize the dv vector data
    const Real
    th_l                   = infile("thickness_lower", 0.0001),
    th_u                   = infile("thickness_upper",   0.2),
    th                     = infile("thickness",         0.2),
    dx                     = _length/(_n_stations_x-1);
    
    _dv_init.resize    (_n_vars);
    _dv_scaling.resize (_n_vars);
    _dv_low.resize     (_n_vars);
    
    // design variables for the thickness values
    for (unsigned int i=0; i<_n_vars; i++) {
        
        _dv_init[i]    =  infile("dv_init", th/th_u, i);
        _dv_low[i]     = th_l/th_u;
        _dv_scaling[i] =      th_u;
    }
    
    
    // create the thickness variables
    _th_station_parameters.resize(_n_vars);
    _th_station_functions.resize(_n_vars);
    
    std::map<Real, MAST::FieldFunction<Real>*> th_station_vals;
    
    for (unsigned int i=0; i<_n_stations_x; i++) {
        std::ostringstream oss;
        oss << "h_" << i;
        
        // now we need a parameter that defines the thickness at the
        // specified station and a constant function that defines the
        // field function at that location.
        MAST::Parameter* h               =
        new MAST::Parameter(oss.str(), infile("thickness", 0.002));
        
        MAST::ConstantFieldFunction* h_f =
        new MAST::ConstantFieldFunction("h", *h);
        
        // add this to the thickness map
        th_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                               (i*dx, h_f));
        
        // add the function to the parameter set
        _th_station_parameters[i]          = h;
        _th_station_functions[i]           = h_f;
        
        // tell the assembly system about the sensitvity parameter
        _discipline->add_parameter(*h);
    }
    
    // now create the h_y function and give it to the property card
    _th_f.reset(new MAST::MultilinearInterpolation("h", th_station_vals));
    
    
    // create the property functions and add them to the
    
    _E               = new MAST::Parameter(   "E", infile("E",    72.e9));
    _nu              = new MAST::Parameter(  "nu", infile("nu",    0.33));
    _kappa           = new MAST::Parameter("kappa",infile("kappa",5./6.));
    _rho             = new MAST::Parameter( "rho", infile("rho", 2700.0));
    _zero            = new MAST::Parameter("zero",                    0.);
    _press           = new MAST::Parameter(   "p", infile("press", 2.e6));
    
    
    _E_f             = new MAST::ConstantFieldFunction("E",            *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
    _rho_f           = new MAST::ConstantFieldFunction("rho",        *_rho);
    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
    
    // initialize the load
    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _p_load->add(*_press_f);
    _discipline->add_volume_load(0, *_p_load);
    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(  *_E_f);
    _m_card->add( *_nu_f);
    _m_card->add(*_rho_f);
    _m_card->add(*_kappa_f);
    
    // create the element property card
    _p_card         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card->add(*_th_f);
    _p_card->add(*_hoff_f);
    
    // tell the section property about the material property
    _p_card->set_material(*_m_card);
    if (if_vk) _p_card->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _discipline->set_property_for_subdomain(0, *_p_card);
    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _mesh->elements_begin(),
    e_end   = _mesh->elements_end();
    
    // points where stress is evaluated
    std::vector<libMesh::Point> pts;

    if (e_type == libMesh::QUAD4 ||
        e_type == libMesh::QUAD8 ||
        e_type == libMesh::QUAD9) {
        
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
        pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
    }
    else if (e_type == libMesh::TRI3 ||
             e_type == libMesh::TRI6) {
        
        pts.push_back(libMesh::Point(1./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(2./3., 1./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(2./3., 1./3.,-1.)); // lower skin
        pts.push_back(libMesh::Point(1./3., 2./3., 1.)); // upper skin
        pts.push_back(libMesh::Point(1./3., 2./3.,-1.)); // lower skin
    }
    else
        libmesh_assert(false); // should not get here


    _outputs = new MAST::StressStrainOutputBase;
    _outputs->set_points_for_evaluation(pts);
    _discipline->add_volume_output((*e_it)->subdomain_id(), *_outputs);
    
    // create the assembly object
    _assembly = new MAST::StructuralNonlinearAssembly;
    
    _assembly->attach_discipline_and_system(*_discipline, *_structural_sys);
    
    
    // create the function to calculate weight
    _weight = new MAST::PlateWeight(*_discipline);
    
    _initialized = true;
}




MAST::PlateBendingSingleStressFunctionalSizingOptimization::
~PlateBendingSingleStressFunctionalSizingOptimization() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card;
        
        delete _p_load;
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        
        delete _E_f;
        delete _nu_f;
        delete _kappa_f;
        delete _hoff_f;
        delete _press_f;
        delete _rho_f;
        
        delete _E;
        delete _nu;
        delete _kappa;
        delete _zero;
        delete _press;
        delete _rho;
        
        delete _weight;
        
        _assembly->clear_discipline_and_system();
        delete _assembly;
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _structural_sys;
        delete _outputs;
        
        
        // delete the h_y station functions
        {
            std::vector<MAST::ConstantFieldFunction*>::iterator
            it  = _th_station_functions.begin(),
            end = _th_station_functions.end();
            for (; it != end; it++)  delete *it;
        }
        
        
        // delete the h_y station parameters
        {
            std::vector<MAST::Parameter*>::iterator
            it  = _th_station_parameters.begin(),
            end = _th_station_parameters.end();
            for (; it != end; it++)  delete *it;
        }
    }
}



void
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
init_dvar(std::vector<Real>& x,
          std::vector<Real>& xmin,
          std::vector<Real>& xmax) {
    
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}



void
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
evaluate(const std::vector<Real>& dvars,
         Real& obj,
         bool eval_obj_grad,
         std::vector<Real>& obj_grad,
         std::vector<Real>& fvals,
         std::vector<bool>& eval_grads,
         std::vector<Real>& grads) {
    
    
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // set the parameter values equal to the DV value
    for (unsigned int i=0; i<_n_vars; i++)
        (*_th_station_parameters[i]) = dvars[i]*_dv_scaling[i];
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;
    
    // the optimization problem is defined as
    // min weight, subject to constraints on displacement and stresses
    Real
    wt      = 0.,
    pval    = 2.;
    
    
    // calculate weight
    (*_weight)(pt, 0., wt);
    
    
    //////////////////////////////////////////////////////////////////////
    // first zero the solution
    //////////////////////////////////////////////////////////////////////
    _sys->solution->zero();
    this->clear_stresss();
    
    
    
    //////////////////////////////////////////////////////////////////////
    // now solve using appropriate number of load steps this load step
    //////////////////////////////////////////////////////////////////////
    bool if_vk = (_p_card->strain_type() == MAST::VON_KARMAN_STRAIN);
    
    // set the number of load steps
    unsigned int
    n_steps = 1;
    if (if_vk) n_steps = 10;
    
    Real
    p0      = (*_press)();
    
    // now iterate over the load steps
    for (unsigned int i=0; i<n_steps; i++) {
        std::cout
        << "Load step: " << i << std::endl;
        
        (*_press)()  =  p0*(i+1.)/(1.*n_steps);
        _sys->solve();
    }
    
    // calculate the stresses
    _assembly->calculate_outputs(*(_sys->solution));
    
    
    //////////////////////////////////////////////////////////////////////
    // get the objective and constraints
    //////////////////////////////////////////////////////////////////////
    
    // now get the displacement constraint
    //pt(0) = 3.;
    //DenseRealVector disp_vec;
    //(*_disp_function)(pt, 0., disp_vec);
    // reference displacement value
    // w < w0 => w/w0 < 1. => w/w0 - 1. < 0
    //disp = disp_vec(0);
    
    // set the function and objective values
    obj = wt;
    
    // copy the element von Mises stress values as the functions
    fvals[0] =  -1. +
    _outputs->von_Mises_p_norm_functional_for_all_elems(pval)/_stress_limit;
    
    
    
    //////////////////////////////////////////////////////////////////
    //   evaluate sensitivity if needed
    //////////////////////////////////////////////////////////////////
    
    // sensitivity of the objective function
    if (eval_obj_grad) {
        
        Real w_sens = 0.;
        
        // set gradient of weight
        for (unsigned int i=0; i<_n_vars; i++) {
            _weight->derivative(MAST::PARTIAL_DERIVATIVE,
                                *_th_station_parameters[i],
                                pt,
                                0.,
                                w_sens);
            obj_grad[i] = w_sens*_dv_scaling[i];
        }
    }
    
    
    // now check if the sensitivity of objective function is requested
    bool if_sens = false;
    
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_sens = (if_sens || eval_grads[i]);
    
    if (if_sens) {
        
        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////
        
        // we are going to choose to use one parametric sensitivity at a time
        for (unsigned int i=0; i<_n_vars; i++) {
            
            libMesh::ParameterVector params;
            params.resize(1);
            params[0]  = _th_station_parameters[i]->ptr();
            
            // iterate over each dv and calculate the sensitivity
            _sys->add_sensitivity_solution(0).zero();
            this->clear_stresss();
            
            // sensitivity analysis
            _sys->sensitivity_solve(params);
            
            // evaluate sensitivity of the outputs
            _assembly->calculate_output_sensitivity(params,
                                                    true,    // true for total sensitivity
                                                    *(_sys->solution));
            
            // copy the sensitivity values in the output
            grads[i] = _dv_scaling[i]/_stress_limit *
            _outputs->von_Mises_p_norm_functional_sensitivity_for_all_elems
            (pval, _th_station_parameters[i]);
        }
    }
    
    
    // write the evaluation output
    //this->output(0, dvars, obj, fvals, false);
}





void
MAST::PlateBendingSingleStressFunctionalSizingOptimization::clear_stresss() {
    
    _outputs->clear(false);
}





void
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
output(unsigned int iter,
       const std::vector<Real>& x,
       Real obj,
       const std::vector<Real>& fval,
       bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
        
    // write the solution for visualization
    libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                        *_eq_sys);
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
}




MAST::FunctionEvaluation::funobj
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
get_objective_evaluation_function() {
    
    return plate_stress_functional_optim_obj;
}



MAST::FunctionEvaluation::funcon
MAST::PlateBendingSingleStressFunctionalSizingOptimization::
get_constraint_evaluation_function() {
    
    return plate_stress_functional_optim_con;
}

