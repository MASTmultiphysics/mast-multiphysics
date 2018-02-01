///*
// * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
// * Copyright (C) 2013-2018  Manav Bhatia
// *
// * This library is free software; you can redistribute it and/or
// * modify it under the terms of the GNU Lesser General Public
// * License as published by the Free Software Foundation; either
// * version 2.1 of the License, or (at your option) any later version.
// *
// * This library is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// * Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with this library; if not, write to the Free Software
// * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// */
//
//// C++ includes
//#include <iostream>
//
//// MAST includes
//#include "examples/structural/topology_optim_2D/topology_optim_2D.h"
//#include "elasticity/stress_output_base.h"
//#include "optimization/optimization_interface.h"
//#include "optimization/function_evaluation.h"
//#include "base/nonlinear_implicit_assembly.h"
//#include "elasticity/structural_nonlinear_assembly.h"
//#include "base/real_output_function.h"
//#include "base/nonlinear_system.h"
//
//
//// libMesh includes
//#include "libmesh/parallel_mesh.h"
//#include "libmesh/mesh_generation.h"
//#include "libmesh/exodusII_io.h"
//#include "libmesh/numeric_vector.h"
//#include "libmesh/sparse_matrix.h"
//#include "libmesh/getpot.h"
//#include "libmesh/string_to_enum.h"
//
//
//extern
//libMesh::LibMeshInit     *__init;
//extern
//MAST::FunctionEvaluation *__my_func_eval;
//
//
//
//void
//topology_optim_2D_optim_obj(int*    mode,
//                            int*    n,
//                            double* x,
//                            double* f,
//                            double* g,
//                            int*    nstate) {
//    
//    
//    // make sure that the global variable has been setup
//    libmesh_assert(__my_func_eval);
//    
//    // initialize the local variables
//    Real
//    obj = 0.;
//    
//    unsigned int
//    n_vars  =  __my_func_eval->n_vars(),
//    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
//    
//    libmesh_assert_equal_to(*n, n_vars);
//    
//    std::vector<Real>
//    dvars   (*n,    0.),
//    obj_grad(*n,    0.),
//    fvals   (n_con, 0.),
//    grads   (0);
//    
//    std::vector<bool>
//    eval_grads(n_con);
//    std::fill(eval_grads.begin(), eval_grads.end(), false);
//    
//    // copy the dvars
//    for (unsigned int i=0; i<n_vars; i++)
//        dvars[i] = x[i];
//    
//    
//    __my_func_eval->evaluate(dvars,
//                             obj,
//                             true,       // request the derivatives of obj
//                             obj_grad,
//                             fvals,
//                             eval_grads,
//                             grads);
//    
//    
//    // now copy them back as necessary
//    *f  = obj;
//    for (unsigned int i=0; i<n_vars; i++)
//        g[i] = obj_grad[i];
//}
//
//
//
//
//
//
//void
//topology_optim_2D_optim_con(int*    mode,
//                            int*    ncnln,
//                            int*    n,
//                            int*    ldJ,
//                            int*    needc,
//                            double* x,
//                            double* c,
//                            double* cJac,
//                            int*    nstate) {
//    
//    
//    // make sure that the global variable has been setup
//    libmesh_assert(__my_func_eval);
//    
//    // initialize the local variables
//    Real
//    obj = 0.;
//    
//    unsigned int
//    n_vars  =  __my_func_eval->n_vars(),
//    n_con   =  __my_func_eval->n_eq()+__my_func_eval->n_ineq();
//    
//    libmesh_assert_equal_to(    *n, n_vars);
//    libmesh_assert_equal_to(*ncnln, n_con);
//    
//    std::vector<Real>
//    dvars   (*n,    0.),
//    obj_grad(*n,    0.),
//    fvals   (n_con, 0.),
//    grads   (n_vars*n_con, 0.);
//    
//    std::vector<bool>
//    eval_grads(n_con);
//    std::fill(eval_grads.begin(), eval_grads.end(), true);
//    
//    // copy the dvars
//    for (unsigned int i=0; i<n_vars; i++)
//        dvars[i] = x[i];
//    
//    
//    __my_func_eval->evaluate(dvars,
//                             obj,
//                             true,       // request the derivatives of obj
//                             obj_grad,
//                             fvals,
//                             eval_grads,
//                             grads);
//    
//    
//    // now copy them back as necessary
//    
//    // first the constraint functions
//    for (unsigned int i=0; i<n_con; i++)
//        c[i] = fvals[i];
//    
//    // next, the constraint gradients
//    for (unsigned int i=0; i<n_con*n_vars; i++)
//        cJac[i] = grads[i];
//    
//    
//}
//
//
//
//
//MAST::YoungsModulus::YoungsModulus(const std::string& nm,
//                                   MAST::FieldFunction<Real> &rho,
//                                   const Real base_modulus,
//                                   const Real penalty):
//MAST::FieldFunction<Real>(nm),
//_rho(rho),
//_base_modulus(base_modulus),
//_penalty(penalty)  {
//    
//    _functions.insert(&rho);
//}
//
//
//
//MAST::YoungsModulus::~YoungsModulus() {
//
//}
//
//
//
//void
//MAST::YoungsModulus::operator() (const libMesh::Point& p,
//                                 Real t,
//                                 Real& v) const {
//    
//    _rho(p, t, v);
//    v = _base_modulus*pow(v, _penalty);
//}
//
//
//void
//MAST::YoungsModulus::derivative(   const MAST::FunctionBase& f,
//                                const libMesh::Point& p,
//                                Real t,
//                                Real& v) const {
//    Real
//    r = 0., dr = 0.;
//    
//    _rho(p, t, r);
//    _rho.derivative( f, p, t, dr);
//    v = _base_modulus*pow(r, _penalty-1.)*_penalty*dr;
//}
//
//
//
//
//
//
//
//MAST::TopologyOptimization2D::
//TopologyOptimization2D(const libMesh::Parallel::Communicator& comm):
//MAST::FunctionEvaluation(comm),
//_initialized(false),
//_penalty(0.),
//_volume_fraction(0.),
//_n_divs_x(0),
//_n_divs_y(0),
//_n_elems(0) { }
//
//
//
//
//void
//MAST::TopologyOptimization2D::init(GetPot& infile,
//                                   libMesh::ElemType e_type,
//                                   bool if_vk) {
//    
//    libmesh_assert(!_initialized);
//    
//    // number of elements
//    _n_divs_x    = infile("n_divs_x", 16);
//    _n_divs_y    = infile("n_divs_y", 16);
//    _n_elems     = _n_divs_x*_n_divs_y;
//    if (e_type == libMesh::TRI3)
//        _n_elems    *= 2;
//    
//    
//    // now setup the optimization data
//    _n_vars                = _n_elems; // for thickness variable
//    _n_eq                  = 0;
//    _n_ineq                = 1;        // volume fraction constraint
//    _max_iters             = 1000;
//    
//    
//    
//    // length of domain
//    _length        = infile("length", 0.50);
//    _width         = infile("width",  0.25);
//    
//    // topology optimization parameters
//    _penalty          = infile(        "penalty",  3.);
//    _volume_fraction  = infile("volume_fraction", 0.3);
//    
//    // create the mesh
//    _mesh          = new libMesh::SerialMesh(__init->comm());
//    
//    // initialize the mesh with one element
//    libMesh::MeshTools::Generation::build_square(*_mesh,
//                                                 _n_divs_x, _n_divs_y,
//                                                 0, _length,
//                                                 0, _width,
//                                                 e_type);
//    _mesh->write("mesh.exo");
//    // create the equation system
//    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
//    
//    // create the libmesh system
//    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
//    _rho_sys   = &(_eq_sys->add_system<libMesh::ExplicitSystem>("density_vars"));
//    
//    // FEType to initialize the system
//    libMesh::FEType fetype     (libMesh::FIRST,    libMesh::LAGRANGE);
//    libMesh::FEType fetype_rho (libMesh::CONSTANT, libMesh::MONOMIAL);
//    
//    
//    // initialize the system to the right set of variables
//    _structural_sys = new MAST::StructuralSystemInitialization(*_sys,
//                                                               _sys->name(),
//                                                               fetype);
//    _discipline     = new MAST::PhysicsDisciplineBase(*_eq_sys);
//
//    // density variable for plotting
//    _rho_sys->add_variable("rho", fetype_rho);
//    
//    // create and add the boundary condition and loads
//    _dirichlet_right  = new MAST::DirichletBoundaryCondition;
//    _dirichlet_left   = new MAST::DirichletBoundaryCondition;
//    
//    _dirichlet_right->init  (1, _structural_sys->vars());
//    _dirichlet_left->init   (3, _structural_sys->vars());
//    
//    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
//    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
//    
//    _discipline->init_system_dirichlet_bc(*_sys);
//    
//    // initialize the equation system
//    _eq_sys->init();
//    
//
//    // create the property functions and add them to the
//    _nu              = new MAST::Parameter(  "nu", infile("nu",    0.33));
//    _kappa           = new MAST::Parameter("kappa",infile("kappa",5./6.));
//    _zero            = new MAST::Parameter("zero",                    0.);
//    _press           = new MAST::Parameter(   "p", infile("press", 2.e6));
//    _th              = new MAST::Parameter(  "th", infile("th",   2.e-3));
//    
//    _nu_f            = new MAST::ConstantFieldFunction("nu",          *_nu);
//    _kappa_f         = new MAST::ConstantFieldFunction("kappa",    *_kappa);
//    _hoff_f          = new MAST::ConstantFieldFunction("off",       *_zero);
//    _press_f         = new MAST::ConstantFieldFunction("pressure", *_press);
//    _th_f            = new MAST::ConstantFieldFunction("h",           *_th);
//    
//    // initialize the load
//    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
//    _p_load->add(*_press_f);
//    _discipline->add_side_load(2, *_p_load); // pressure on upper surface
//
//    
//    // create the compliance output object
//    _output = new MAST::RealOutputFunction(MAST::STRUCTURAL_COMPLIANCE);
//    _discipline->add_volume_output(0, *_output);
//    
//    
//    // resize the elem vector
//    _elems.resize(_n_elems);
//
//    // element iterators to define property cards
//    libMesh::MeshBase::element_iterator
//    el_it  = _mesh->local_elements_begin(),
//    el_end = _mesh->local_elements_end();
//    
//    unsigned int
//    counter = 0;
//    
//    bool
//    insert_success = false;
//    
//    for (; el_it != el_end; el_it++) {
//        
//        libMesh::Elem* el = *el_it;
//
//        // set the subdomain ID of the element so that the
//        // property card can be associated with this ID
//        el->subdomain_id()      = counter;
//        
//        std::ostringstream oss;
//        oss << "rho_" << counter;
//
//        
//        // create the density parameters for the elements
//        MAST::Parameter*
//        rho    = new MAST::Parameter(oss.str(), 1.);
//        MAST::ConstantFieldFunction*
//        rho_f  = new MAST::ConstantFieldFunction("rho", *rho);
//        MAST::YoungsModulus*
//        E_f    = new MAST::YoungsModulus("E",
//                                         *rho_f,
//                                         70.e9,
//                                         _penalty);
//        
//        
//        // create the material property card
//        MAST::IsotropicMaterialPropertyCard*
//        m_card         = new MAST::IsotropicMaterialPropertyCard;
//        
//        // add the material properties to the card
//        m_card->add(     *E_f);
//        m_card->add(   *_nu_f);
//        m_card->add(   *rho_f);
//        m_card->add(*_kappa_f);
//        
//        // create the element property card
//        MAST::Solid2DSectionElementPropertyCard*
//        p_card         = new MAST::Solid2DSectionElementPropertyCard;
//        
//        // add the section properties to the card
//        p_card->add(*_th_f);
//        p_card->add(*_hoff_f);
//        
//        // tell the section property about the material property
//        p_card->set_material(*m_card);
//        
//        _discipline->set_property_for_subdomain(counter, *p_card);
//        
//        
//        // now add the property and card pointers to the map
//        _elems[counter]         =  el;
//        
//        insert_success =
//        _elem_rho.insert(std::pair<const libMesh::Elem*, MAST::Parameter*>
//                         (el, rho)).second;
//        
//        libmesh_assert(insert_success);
//
//        insert_success =
//        _elem_rho_f.insert(std::pair<const libMesh::Elem*, MAST::ConstantFieldFunction*>
//                           (el, rho_f)).second;
//        
//        libmesh_assert(insert_success);
//
//        
//        insert_success =
//        _elem_E_f.insert(std::pair<const libMesh::Elem*, MAST::YoungsModulus*>
//                         (el, E_f)).second;
//        
//        libmesh_assert(insert_success);
//
//
//        insert_success =
//        _elem_m_card.insert(std::pair<const libMesh::Elem*, MAST::IsotropicMaterialPropertyCard*>
//                            (el, m_card)).second;
//        
//        libmesh_assert(insert_success);
//
//        
//        insert_success =
//        _elem_p_card.insert(std::pair<const libMesh::Elem*, MAST::Solid2DSectionElementPropertyCard*>
//                            (el, p_card)).second;
//        
//        libmesh_assert(insert_success);
//
//        counter++;
//    }
//
//    
//    // create the assembly object
//    _assembly = new MAST::NonlinearImplicitAssembly;
//    _elem_ops = new MAST::StructuralNonlinearAssemblyElemOperations;
//    
//    _assembly->attach_discipline_and_system(*_elem_ops,
//                                            *_discipline,
//                                            *_structural_sys);
//    
//    
//    _initialized = true;
//}
//
//
//
//
//MAST::TopologyOptimization2D::~TopologyOptimization2D() {
//    
//    
//    if (_initialized) {
//
//        delete _p_load;
//        delete _dirichlet_right;
//        delete _dirichlet_left;
//        
//        delete _nu_f;
//        delete _kappa_f;
//        delete _hoff_f;
//        delete _press_f;
//        delete _th_f;
//        
//        delete _nu;
//        delete _kappa;
//        delete _zero;
//        delete _press;
//        delete _th;
//        
//        _assembly->clear_discipline_and_system();
//        delete _assembly;
//        delete _elem_ops;
//
//        delete _eq_sys;
//        delete _mesh;
//        
//        delete _discipline;
//        delete _structural_sys;
//        
//        delete _output;
//        
//        // delete the element data
//        {
//            std::map<const libMesh::Elem*, MAST::Parameter*>::iterator
//            it   = _elem_rho.begin(),
//            end  = _elem_rho.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//
//        {
//            std::map<const libMesh::Elem*, MAST::ConstantFieldFunction*>::iterator
//            it   = _elem_rho_f.begin(),
//            end  = _elem_rho_f.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::YoungsModulus*>::iterator
//            it   = _elem_E_f.begin(),
//            end  = _elem_E_f.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::IsotropicMaterialPropertyCard*>::iterator
//            it   = _elem_m_card.begin(),
//            end  = _elem_m_card.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//        
//        {
//            std::map<const libMesh::Elem*, MAST::Solid2DSectionElementPropertyCard*>::iterator
//            it   = _elem_p_card.begin(),
//            end  = _elem_p_card.end();
//            
//            for (; it != end; it++)
//                delete it->second;
//        }
//
//    }
//}
//
//
//
//void
//MAST::TopologyOptimization2D::init_dvar(std::vector<Real>& x,
//                                                          std::vector<Real>& xmin,
//                                                          std::vector<Real>& xmax) {
//    // one DV for each element
//    x.resize(_n_vars);
//    xmin.resize(_n_vars);
//    xmax.resize(_n_vars);
//    
//    std::fill(   x.begin(),    x.end(),   0.3);
//    std::fill(xmin.begin(), xmin.end(), 1.e-5);
//    std::fill(xmax.begin(), xmax.end(),    1.);
//}
//
//
//
//void
//MAST::TopologyOptimization2D::evaluate(const std::vector<Real>& dvars,
//                                                         Real& obj,
//                                                         bool eval_obj_grad,
//                                                         std::vector<Real>& obj_grad,
//                                                         std::vector<Real>& fvals,
//                                                         std::vector<bool>& eval_grads,
//                                                         std::vector<Real>& grads) {
//    
//    
//    libmesh_assert_equal_to(dvars.size(), _n_vars);
//    
//    // set the parameter values equal to the DV value
//    for (unsigned int i=0; i<_n_vars; i++) {
//        const libMesh::Elem* el = _elems[i];
//        MAST::Parameter* par = _elem_rho[el];
//        (*par)() = dvars[i];
//    }
//    
//    // DO NOT zero out the gradient vector, since GCMMA needs it for the
//    // subproblem solution
//    
//    libMesh::Point pt; // dummy point object
//    
//    libMesh::out << "New Eval" << std::endl;
//    
//    // the optimization problem is defined as
//    // min weight, subject to constraints on displacement and stresses
//    
//    
//    
//    //////////////////////////////////////////////////////////////////////
//    // first zero the solution
//    //////////////////////////////////////////////////////////////////////
//    _sys->solution->zero();
//    _output->clear();
//    _sys->solve();
//
//
//    //////////////////////////////////////////////////////////////////////
//    // get the objective and constraints
//    //////////////////////////////////////////////////////////////////////
//    
//    // now get the compliance
//    _assembly->calculate_outputs(*_sys->solution);
//    obj = _output->get_value(); // negative, since J = -K
//    
//    
//    fvals[0] = 0.;
//    Real
//    vol       = 0.,
//    total_vol = 0.;
//    for (unsigned int i=0; i<_n_elems; i++) {
//        vol = _elems[i]->volume();
//        fvals[0]  += dvars[i] * vol; // constraint:  xi vi - V <= 0
//        total_vol += vol;
//    }
//    fvals[0] /= (total_vol * _volume_fraction);
//    fvals[0] -=              1.;
//    
//    
//    
//    //////////////////////////////////////////////////////////////////
//    //   evaluate sensitivity if needed
//    //////////////////////////////////////////////////////////////////
//    
//    // sensitivity of the objective function
//    if (eval_obj_grad) {
//        
//        // we are going to choose to use one parametric sensitivity at a time
//        for (unsigned int i=0; i<_n_elems; i++) {
//            
//            MAST::Parameter* f = _elem_rho[_elems[i]];
//            
//            libMesh::ParameterVector params;
//            params.resize(1);
//            params[0]  = f->ptr();
//            
//            // iterate over each dv and calculate the sensitivity
//            _sys->add_sensitivity_solution(0).zero();
//            
//            // sensitivity analysis
//            _sys->sensitivity_solve(params);
//            _assembly->calculate_output_sensitivity(params,
//                                                    true,
//                                                    *_sys->solution);
//            
//            obj_grad[i] = _output->get_sensitivity(f);
//        }
//    }
//    
//    // now check if the sensitivity of constraint function is requested
//    if (eval_grads[0]) {
//        
//        //////////////////////////////////////////////////////////////////
//        // indices used by GCMMA follow this rule:
//        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
//        //////////////////////////////////////////////////////////////////
//        
//        for (unsigned int i=0; i<_n_elems; i++)
//            grads[i] = _elems[i]->volume()/(_volume_fraction * total_vol);
//    }
//    
//    
//    // write the evaluation output
//    //this->output(0, dvars, obj, fvals, false);
//}
//
//
//
//
//
//
//void
//MAST::TopologyOptimization2D::output(unsigned int iter,
//                                                       const std::vector<Real>& x,
//                                                       Real obj,
//                                                       const std::vector<Real>& fval,
//                                                       bool if_write_to_optim_file) const {
//    
//    libmesh_assert_equal_to(x.size(), _n_vars);
//    
//    // set the desity value in the auxiliary system for output
//    for (unsigned int i=0; i<_n_vars; i++) {
//        const libMesh::Elem* el = _elems[i];
//        _rho_sys->solution->set(el->dof_number(_rho_sys->number(), 0, 0), x[i]);
//    }
//    
//    
//    // write the solution for visualization
//    std::set<std::string> nm;
//    nm.insert(_sys->name());
//    nm.insert(_rho_sys->name());
//    libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
//                                                        *_eq_sys,
//                                                        &nm);
//    
//    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
//}
//
//
//
//
//MAST::FunctionEvaluation::funobj
//MAST::TopologyOptimization2D::get_objective_evaluation_function() {
//    
//    return topology_optim_2D_optim_obj;
//}
//
//
//
//MAST::FunctionEvaluation::funcon
//MAST::TopologyOptimization2D::get_constraint_evaluation_function() {
//    
//    return topology_optim_2D_optim_con;
//}
//
