/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2016  Manav Bhatia
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
#include "examples/structural/stiffened_plate_optimization_thermally_stressed_piston_theory_flutter/stiffened_plate_thermally_stressed_piston_theory_flutter_optimization.h"
#include "examples/structural/stiffened_plate_optimization/stiffened_plate_optimization_base.h"
#include "examples/structural/beam_bending/beam_bending.h"
#include "driver/driver_base.h"
#include "base/nonlinear_system.h"
#include "elasticity/stress_output_base.h"
#include "optimization/optimization_interface.h"
#include "optimization/function_evaluation.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "elasticity/piston_theory_boundary_condition.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
#include "elasticity/structural_fluid_interaction_assembly.h"
#include "elasticity/structural_near_null_vector_space.h"
#include "aeroelasticity/time_domain_flutter_solver.h"
#include "aeroelasticity/time_domain_flutter_root.h"


// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_solver.h"

extern
libMesh::LibMeshInit     *__init;
extern
MAST::FunctionEvaluation *__my_func_eval;



void
stiffened_plate_thermally_stressed_piston_theory_flutter_optim_obj(int*    mode,
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
stiffened_plate_thermally_stressed_piston_theory_flutter_optim_con(int*    mode,
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





MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
StiffenedPlateThermallyStressedPistonTheorySizingOptimization
(const libMesh::Parallel::Communicator& comm):
MAST::FunctionEvaluation(comm),
_initialized(false),
_n_eig(0),
_n_divs_x(0),
_n_divs_between_stiff(0),
_n_stiff(0),
_n_plate_elems(0),
_n_elems_per_stiff(0),
_n_elems(0),
_n_dv_stations_x(0),
_n_load_steps(0),
_mesh(NULL),
_eq_sys(NULL),
_sys(NULL),
_structural_sys(NULL),
_discipline(NULL),
_nonlinear_assembly(NULL),
_modal_assembly(NULL),
_fsi_assembly(NULL),
_E(NULL),
_alpha(NULL),
_nu(NULL),
_kappa(NULL),
_rho(NULL),
_temp(NULL),
_zero(NULL),
_velocity(NULL),
_mach(NULL),
_rho_air(NULL),
_gamma_air(NULL),
_E_f(NULL),
_nu_f(NULL),
_alpha_f(NULL),
_kappa_f(NULL),
_rho_f(NULL),
_temp_f(NULL),
_thzoff_stiff_f(NULL),
_ref_temp_f(NULL),
_velocity_f(NULL),
_mach_f(NULL),
_rho_air_f(NULL),
_gamma_air_f(NULL),
_hoff_plate_f(NULL),
_weight(NULL),
_m_card(NULL),
_p_card_plate(NULL),
_dirichlet_left(NULL),
_dirichlet_right(NULL),
_dirichlet_bottom(NULL),
_dirichlet_top(NULL),
_T_load(NULL),
_piston_bc(NULL),
_flutter_solver(NULL)
{ }




void
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
init(GetPot& infile,
     libMesh::ElemType e_type,
     bool if_vk) {
    
    libmesh_assert(!_initialized);
    
    // number of eigenvalues
    _n_eig = 20;
    
    // number of elements
    _n_divs_x                = infile("n_divs_x",             40);
    _n_divs_between_stiff    = infile("n_divs_between_stiff", 10);
    _n_stiff                 = infile("n_stiffeners",          3);
    
    _n_plate_elems           = _n_divs_x*(_n_stiff+1)*_n_divs_between_stiff;
    if (e_type == libMesh::TRI3)
        _n_plate_elems    *= 2;
    
    _n_elems_per_stiff       = _n_divs_x;
    _n_elems                 = _n_plate_elems + _n_stiff * _n_elems_per_stiff;
    
    
    // number of stations
    _n_dv_stations_x         = infile("n_stations", 4);
    
    // number of load steps
    _n_load_steps            = infile("n_load_steps", 20);
    
    
    // now setup the optimization data
    _n_vars                = _n_dv_stations_x + 2*_n_dv_stations_x * _n_stiff; // for thickness variable
    _n_eq                  = 0;
    _n_ineq                = _n_eig + 1 + _n_elems; // constraint that each eigenvalue > 0 flutter constraint + one element stress functional per elem
    _max_iters             = 1000;
    
    
    
    // length of domain
    _length        = infile("length", 0.50);
    _width         = infile("width",  0.25);
    
    // limit stress
    _stress_limit  = infile("max_stress", 4.00e8);
    
    // create the mesh
    _mesh          = new libMesh::SerialMesh(__init->comm());
    
    // initialize the mesh with one element
    MAST::StiffenedPanelMesh panel_mesh;
    panel_mesh.init(_n_stiff,
                    _n_divs_x,
                    _n_divs_between_stiff,
                    _length,
                    _width,
                    *_mesh,
                    e_type,
                    true);
    
    // create the equation system
    _eq_sys    = new  libMesh::EquationSystems(*_mesh);
    
    // create the libmesh system
    _sys       = &(_eq_sys->add_system<MAST::NonlinearSystem>("structural"));
    _sys->set_eigenproblem_type(libMesh::GHEP);
    
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
    
    //_dirichlet_bottom->init (0, constrained_vars);
    //_dirichlet_right->init  (1, constrained_vars);
    //_dirichlet_top->init    (2, constrained_vars);
    //_dirichlet_left->init   (3, constrained_vars);
    
    _dirichlet_bottom->init (0, _structural_sys->vars());
    _dirichlet_right->init  (1, _structural_sys->vars());
    _dirichlet_top->init    (2, _structural_sys->vars());
    _dirichlet_left->init   (3, _structural_sys->vars());
    
    _discipline->add_dirichlet_bc(0, *_dirichlet_bottom);
    _discipline->add_dirichlet_bc(1,  *_dirichlet_right);
    _discipline->add_dirichlet_bc(2,    *_dirichlet_top);
    _discipline->add_dirichlet_bc(3,   *_dirichlet_left);
    
    _discipline->init_system_dirichlet_bc(*_sys);
    
    // initialize the equation system
    _eq_sys->init();
    _sys->initialize_condensed_dofs(*_discipline);
    _sys->eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    _sys->set_exchange_A_and_B(true);
    _sys->set_n_requested_eigenvalues(_n_eig);
    
    // initialize the dv vector data
    const Real
    th_l                   = infile("thickness_lower",  0.001),
    th_u                   = infile("thickness_upper",    0.2),
    th                     = infile("thickness",          0.2),
    dx                     = _length/(_n_dv_stations_x-1);
    
    _dv_init.resize    (_n_vars);
    _dv_scaling.resize (_n_vars);
    _dv_low.resize     (_n_vars);
    _problem_parameters.resize(_n_vars);
    
    // design variables for the thickness values
    for (unsigned int i=0; i<_n_vars; i++) {
        
        _dv_init[i]    =  infile("dv_init", th/th_u, i);
        _dv_low[i]     =  th_l/th_u;
        _dv_scaling[i] =      th_u;
    }
    
    
    // create the thickness variables
    _th_station_parameters_plate.resize(_n_dv_stations_x);
    _th_station_functions_plate.resize(_n_dv_stations_x);
    
    std::map<Real, MAST::FieldFunction<Real>*>
    thy_station_vals,
    thz_station_vals;
    
    for (unsigned int i=0; i<_n_dv_stations_x; i++) {
        std::ostringstream oss;
        oss << "h_" << i;
        
        // now we need a parameter that defines the thickness at the
        // specified station and a constant function that defines the
        // field function at that location.
        MAST::Parameter* h               =
        new MAST::Parameter(oss.str(), infile("thickness", 0.002));
        
        MAST::ConstantFieldFunction* h_f =
        new MAST::ConstantFieldFunction(oss.str(), *h);
        
        // add this to the thickness map
        thy_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                (i*dx, h_f));
        
        // add the function to the parameter set
        _th_station_parameters_plate[i]          = h;
        _th_station_functions_plate[i]           = h_f;
        
        // tell the assembly system about the sensitvity parameter
        _discipline->add_parameter(*h);
        _problem_parameters[i] = h;
    }
    
    // now create the h_y function and give it to the property card
    _th_plate_f.reset(new MAST::MultilinearInterpolation("h", thy_station_vals));
    thy_station_vals.clear();
    
    
    // create the property functions and add them to the card
    
    _E               = new MAST::Parameter(   "E", infile("E",         72.e9));
    _nu              = new MAST::Parameter(  "nu", infile("nu",         0.33));
    _kappa           = new MAST::Parameter("kappa",infile("kappa",     5./6.));
    _alpha           = new MAST::Parameter("alpha",infile("alpha",    2.5e-5));
    _rho             = new MAST::Parameter( "rho", infile("rho",      2700.0));
    _zero            = new MAST::Parameter("zero",                         0.);
    _temp            = new MAST::Parameter( "temperature",infile("temp", 60.));
    _p_cav           = new MAST::Parameter("p_cav",infile("p_cav",        0.));
    _velocity        = new MAST::Parameter("V"   ,                         0.);
    _mach            = new MAST::Parameter("mach", infile("mach",        4.5));
    _rho_air         = new MAST::Parameter("rho" , infile("rho_f",      1.05));
    _gamma_air       = new MAST::Parameter("gamma", infile("gamma",      1.4));
    
    
    _E_f             = new MAST::ConstantFieldFunction("E",                   *_E);
    _nu_f            = new MAST::ConstantFieldFunction("nu",                 *_nu);
    _kappa_f         = new MAST::ConstantFieldFunction("kappa",           *_kappa);
    _alpha_f         = new MAST::ConstantFieldFunction("alpha_expansion", *_alpha);
    _rho_f           = new MAST::ConstantFieldFunction("rho",               *_rho);
    _temp_f          = new MAST::ConstantFieldFunction("temperature",      *_temp);
    _p_cav_f         = new MAST::ConstantFieldFunction("pressure",        *_p_cav);
    _ref_temp_f      = new MAST::ConstantFieldFunction("ref_temperature",  *_zero);
    _thzoff_stiff_f  = new MAST::ConstantFieldFunction("hz_off",           *_zero);
    _hoff_plate_f    = new MAST::SectionOffset("off",
                                               *_th_plate_f,
                                               0.);
    _velocity_f      = new MAST::ConstantFieldFunction("V",      *_velocity);
    _mach_f          = new MAST::ConstantFieldFunction("mach",       *_mach);
    _rho_air_f       = new MAST::ConstantFieldFunction("rho",     *_rho_air);
    _gamma_air_f     = new MAST::ConstantFieldFunction("gamma", *_gamma_air);
    
    // velocity constraint for flutter analysis
    _V0_flutter      =  infile("V0_flutter",     410.);
    _n_V_divs_flutter=  infile("n_V_divs",        10);
    
    // initialize the load
    _T_load          = new MAST::BoundaryConditionBase(MAST::TEMPERATURE);
    _T_load->add(*_temp_f);
    _T_load->add(*_ref_temp_f);
    _discipline->add_volume_load(0, *_T_load);          // for the panel
    for (unsigned int i=0; i<_n_stiff; i++)
        _discipline->add_volume_load(i+1, *_T_load);    // for the stiffeners

    // pressure load
    _p_load          = new MAST::BoundaryConditionBase(MAST::SURFACE_PRESSURE);
    _p_load->add(*_p_cav_f);
    _discipline->add_volume_load(0, *_p_load);          // for the panel

    
    // create the material property card
    _m_card         = new MAST::IsotropicMaterialPropertyCard;
    
    // add the material properties to the card
    _m_card->add(  *_E_f);
    _m_card->add( *_nu_f);
    _m_card->add(*_rho_f);
    _m_card->add(*_kappa_f);
    _m_card->add(*_alpha_f);
    
    // create the element property card
    _p_card_plate         = new MAST::Solid2DSectionElementPropertyCard;
    
    // add the section properties to the card
    _p_card_plate->add(*_th_plate_f);
    _p_card_plate->add(*_hoff_plate_f);
    
    // tell the section property about the material property
    _p_card_plate->set_material(*_m_card);
    if (if_vk) _p_card_plate->set_strain(MAST::VON_KARMAN_STRAIN);
    
    _discipline->set_property_for_subdomain(0, *_p_card_plate);
    
    // now initialize the piston theory boundary conditions
    RealVectorX  vel = RealVectorX::Zero(3);
    vel(0)           = 1.;  // flow along the x-axis
    _piston_bc       = new MAST::PistonTheoryBoundaryCondition(1,     // order
                                                               vel);  // vel vector
    _piston_bc->add(*_velocity_f);
    _piston_bc->add(*_mach_f);
    _piston_bc->add(*_rho_air_f);
    _piston_bc->add(*_gamma_air_f);
    _discipline->add_volume_load(0, *_piston_bc);
    _discipline->add_parameter(*_velocity);
    
    // now add the property cards for each stiffener
    // element orientation
    libMesh::Point orientation;
    orientation(2) = 1.;
    
    // property card per stiffener
    _p_card_stiff.resize(_n_stiff);
    
    // thickness per stiffener station
    _thy_station_parameters_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thy_station_functions_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thz_station_parameters_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thz_station_functions_stiff.resize(_n_dv_stations_x*_n_stiff);
    _thy_stiff_f.resize(_n_stiff);
    _thz_stiff_f.resize(_n_stiff);
    _hyoff_stiff_f.resize(_n_stiff);
    
    for (unsigned int i=0; i<_n_stiff; i++) {
        
        // this map is used to store the thickness parameter along length
        thy_station_vals.clear();
        thz_station_vals.clear();
        
        // first define the thickness station parameters and the thickness
        // field function
        for (unsigned int j=0; j<_n_dv_stations_x; j++) {
            std::ostringstream ossy, ossz;
            ossy << "h_y_" << j << "_stiff_" << i;
            ossz << "h_y_" << j << "_stiff_" << i;
            
            // now we need a parameter that defines the thickness at the
            // specified station and a constant function that defines the
            // field function at that location.
            MAST::Parameter
            *h_y  = new MAST::Parameter(ossy.str(), infile("thickness", 0.002)),
            *h_z  = new MAST::Parameter(ossz.str(), infile("thickness", 0.002));
            
            MAST::ConstantFieldFunction
            *h_y_f = new MAST::ConstantFieldFunction(ossy.str(), *h_y),
            *h_z_f = new MAST::ConstantFieldFunction(ossy.str(), *h_z);
            
            // add this to the thickness map
            thy_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                    (j*dx, h_y_f));
            thz_station_vals.insert(std::pair<Real, MAST::FieldFunction<Real>*>
                                    (j*dx, h_z_f));
            
            // add the function to the parameter set
            _thy_station_parameters_stiff[i*_n_dv_stations_x+j]          = h_y;
            _thy_station_functions_stiff [i*_n_dv_stations_x+j]          = h_y_f;
            _thz_station_parameters_stiff[i*_n_dv_stations_x+j]          = h_z;
            _thz_station_functions_stiff [i*_n_dv_stations_x+j]          = h_z_f;
            
            
            // tell the assembly system about the sensitvity parameter
            _discipline->add_parameter(*h_y);
            _discipline->add_parameter(*h_z);
            _problem_parameters[(2*i+1)*_n_dv_stations_x+j] = h_y;
            _problem_parameters[(2*i+2)*_n_dv_stations_x+j] = h_z;
        }
        
        // now create the h_y function and give it to the property card
        _thy_stiff_f[i]   = new MAST::MultilinearInterpolation("hy", thy_station_vals);
        _thz_stiff_f[i]   = new MAST::MultilinearInterpolation("hz", thz_station_vals);
        _hyoff_stiff_f[i] = new MAST::SectionOffset("hy_off",
                                                    *_thy_stiff_f[i],
                                                    -1.);
        
        _p_card_stiff[i]  = new MAST::Solid1DSectionElementPropertyCard;
        
        
        // add the section properties to the card
        _p_card_stiff[i]->add(*_thy_stiff_f[i]);
        _p_card_stiff[i]->add(*_thz_stiff_f[i]);
        _p_card_stiff[i]->add(*_hyoff_stiff_f[i]);
        _p_card_stiff[i]->add(*_thzoff_stiff_f);
        _p_card_stiff[i]->y_vector() = orientation;
        
        // tell the section property about the material property
        _p_card_stiff[i]->set_material(*_m_card);
        if (if_vk) _p_card_stiff[i]->set_strain(MAST::VON_KARMAN_STRAIN);
        
        _p_card_stiff[i]->init();
        
        // the domain ID of the stiffener is 1 plus the stiff number
        _discipline->set_property_for_subdomain(i+1, *_p_card_stiff[i]);
    }
    
    // initialize the null space object and assign it to the structural
    // module
    _nsp = new MAST::StructuralNearNullVectorSpace;
    _sys->nonlinear_solver->nearnullspace_object = _nsp;

    
    
    // create the output objects, one for each element
    libMesh::MeshBase::const_element_iterator
    e_it    = _mesh->local_elements_begin(),
    e_end   = _mesh->local_elements_end();

    // points where stress is evaluated
    std::vector<libMesh::Point> pts;
    
    for ( ; e_it != e_end; e_it++) {
        
        pts.clear();
        if ((*e_it)->type() == libMesh::QUAD4 ||
            (*e_it)->type() == libMesh::QUAD8 ||
            (*e_it)->type() == libMesh::QUAD9) {
            
            pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3), -1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
            pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3), 1.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3),  1/sqrt(3),-1.)); // lower skin
        }
        else if ((*e_it)->type() == libMesh::TRI3 ||
                 (*e_it)->type() == libMesh::TRI6) {
            
            pts.push_back(libMesh::Point(1./3., 1./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(1./3., 1./3.,-1.)); // lower skin
            pts.push_back(libMesh::Point(2./3., 1./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(2./3., 1./3.,-1.)); // lower skin
            pts.push_back(libMesh::Point(1./3., 2./3., 1.)); // upper skin
            pts.push_back(libMesh::Point(1./3., 2./3.,-1.)); // lower skin
        }
        else if ((*e_it)->type() == libMesh::EDGE2 ||
                 (*e_it)->type() == libMesh::EDGE3) {
            
            pts.push_back(libMesh::Point(-1/sqrt(3), 1., 0.)); // upper skin
            pts.push_back(libMesh::Point(-1/sqrt(3),-1., 0.)); // lower skin
            pts.push_back(libMesh::Point( 1/sqrt(3), 1., 0.)); // upper skin
            pts.push_back(libMesh::Point( 1/sqrt(3),-1., 0.)); // lower skin
        }
        else
            libmesh_assert(false); // should not get here
        
        
        MAST::StressStrainOutputBase * output = new MAST::StressStrainOutputBase;
        
        // tell the object to evaluate the data for this object only
        std::set<const libMesh::Elem*> e_set;
        e_set.insert(*e_it);
        output->set_elements_in_domain(e_set);
        output->set_points_for_evaluation(pts);
        output->set_volume_loads(_discipline->volume_loads());
        _outputs.push_back(output);
        
        _discipline->add_volume_output((*e_it)->subdomain_id(), *output);
    }
    
    // make sure that the number of output objects here is the same as the
    // number of local elems in the mesh
    libmesh_assert_equal_to(_outputs.size(), _mesh->n_local_elem());
    
    // create the assembly object
    _nonlinear_assembly = new MAST::StructuralNonlinearAssembly;
    _modal_assembly     = new MAST::StructuralModalEigenproblemAssembly;
    _fsi_assembly       = new MAST::StructuralFluidInteractionAssembly;
    
    // flutter solver
    _flutter_solver  = new MAST::TimeDomainFlutterSolver;
    std::string nm("flutter_output.txt");
    if (__init->comm().rank() == 0)
        _flutter_solver->set_output_file(nm);
    
    
    
    // create the function to calculate weight
    _weight = new MAST::StiffenedPlateWeight(*_discipline);
    
    _initialized = true;
}




MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
~StiffenedPlateThermallyStressedPistonTheorySizingOptimization() {
    
    if (_initialized) {
        
        delete _m_card;
        delete _p_card_plate;
        for (unsigned int i=0; i<_n_stiff; i++) delete _p_card_stiff[i];
        
        delete _T_load;
        delete _p_load;
        
        delete _dirichlet_bottom;
        delete _dirichlet_right;
        delete _dirichlet_top;
        delete _dirichlet_left;
        
        delete _E_f;
        delete _alpha_f;
        delete _nu_f;
        delete _kappa_f;
        delete _rho_f;
        delete _hoff_plate_f;
        delete _temp_f;
        delete _p_cav_f;
        delete _ref_temp_f;
        delete _thzoff_stiff_f;
        for (unsigned int i=0; i<_n_stiff; i++) delete _hyoff_stiff_f[i];
        for (unsigned int i=0; i<_n_stiff; i++) delete _thy_stiff_f[i];
        for (unsigned int i=0; i<_n_stiff; i++) delete _thz_stiff_f[i];
        delete _velocity_f;
        delete _mach_f;
        delete _rho_air_f;
        delete _gamma_air_f;
        
        
        delete _E;
        delete _nu;
        delete _alpha;
        delete _kappa;
        delete _rho;
        delete _zero;
        delete _temp;
        delete _p_cav;
        delete _velocity;
        delete _mach;
        delete _rho_air;
        delete _gamma_air;
        
        
        delete _weight;
        
        delete _nonlinear_assembly;
        delete _modal_assembly;
        delete _fsi_assembly;
        
        // delete the basis vectors
        if (_basis.size())
            for (unsigned int i=0; i<_basis.size(); i++)
                delete _basis[i];
        
        delete _eq_sys;
        delete _mesh;
        
        delete _discipline;
        delete _structural_sys;
        
        delete _flutter_solver;
        delete _piston_bc;
        
        delete _nsp;
        
        
        // iterate over the output quantities and delete them
        {
            std::vector<MAST::StressStrainOutputBase*>::iterator
            it   =   _outputs.begin(),
            end  =   _outputs.end();
            for ( ; it != end; it++) delete *it;
            
            _outputs.clear();
        }
        
        
        // delete the h_y station functions
        {
            std::vector<MAST::ConstantFieldFunction*>::iterator
            it  = _th_station_functions_plate.begin(),
            end = _th_station_functions_plate.end();
            for (; it != end; it++)  delete *it;
        }
        
        
        {
            std::vector<MAST::ConstantFieldFunction*>::iterator
            it  = _thy_station_functions_stiff.begin(),
            end = _thy_station_functions_stiff.end();
            for (; it != end; it++)  delete *it;
        }

        {
            std::vector<MAST::ConstantFieldFunction*>::iterator
            it  = _thz_station_functions_stiff.begin(),
            end = _thz_station_functions_stiff.end();
            for (; it != end; it++)  delete *it;
        }

        
        // delete the h_y station parameters
        {
            std::vector<MAST::Parameter*>::iterator
            it  = _th_station_parameters_plate.begin(),
            end = _th_station_parameters_plate.end();
            for (; it != end; it++)  delete *it;
        }
        
        
        {
            std::vector<MAST::Parameter*>::iterator
            it  = _thy_station_parameters_stiff.begin(),
            end = _thy_station_parameters_stiff.end();
            for (; it != end; it++)  delete *it;
        }

        {
            std::vector<MAST::Parameter*>::iterator
            it  = _thz_station_parameters_stiff.begin(),
            end = _thz_station_parameters_stiff.end();
            for (; it != end; it++)  delete *it;
        }
    }
}



void
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
init_dvar(std::vector<Real>& x,
          std::vector<Real>& xmin,
          std::vector<Real>& xmax) {
    
    // one DV for each element
    x       = _dv_init;
    xmin    = _dv_low;
    xmax.resize(_n_vars);
    std::fill(xmax.begin(), xmax.end(), 1.);
}


// define the class to provide the interface for nonlinear steady-state
// solution that can be used by the flutter solver to solve for different
// flight velocities
namespace MAST {
    
    class StiffenedPlateSteadySolverInterface:
    public MAST::FlutterSolverBase::SteadySolver {
        
    public:
        StiffenedPlateSteadySolverInterface(MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization& obj,
                              bool if_output,
                              bool if_clear_vector_on_exit,
                              unsigned int n_steps):
        MAST::FlutterSolverBase::SteadySolver(),
        _obj(obj),
        _if_write_output(if_output),
        _n_steps(n_steps),
        _if_only_aero_load_steps(false),
        _if_clear_vector_on_exit(if_clear_vector_on_exit),
        inc (0) {
            
            _obj._sys->add_vector("base_solution");
            writer = new libMesh::ExodusII_IO(_obj._sys->get_mesh());
        }
        
        virtual ~StiffenedPlateSteadySolverInterface() {
            
            if (_if_clear_vector_on_exit)
                _obj._sys->remove_vector("base_solution");
            delete writer;
        }
        
        
        /*!
         *  solves for the steady state solution, and @returns
         *  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solve() {
            
            libMesh::NumericVector<Real>& sol = _obj._sys->get_vector("base_solution");
            *_obj._sys->solution = sol;
            
            // now do the solve
            libmesh_assert(_obj._initialized);
            
            bool if_vk = (_obj._p_card_plate->strain_type() == MAST::VON_KARMAN_STRAIN);
            
            ///////////////////////////////////////////////////////////////
            // first, solve the quasi-steady problem
            ///////////////////////////////////////////////////////////////
            // set the number of load steps
            unsigned int
            n_steps = 1;
            if (if_vk) n_steps = _n_steps;
            
            Real
            T0      = (*_obj._temp)(),
            V0      = (*_obj._velocity)(),
            p0      = (*_obj._p_cav)();
            
            // zero the solution before solving
            _obj.clear_stresss();
            
            _obj._nonlinear_assembly->attach_discipline_and_system(*_obj._discipline,
                                                                   *_obj._structural_sys);
            
            // now iterate over the load steps
            for (unsigned int i=0; i<n_steps; i++) {
                
                // modify aero component
                (*_obj._velocity)()  =  V0*(i+1.)/(1.*n_steps);
                
                // modify the thermal load if specified by the user
                if (!_if_only_aero_load_steps) {
                    
                    (*_obj._temp)()      =  T0*(i+1.)/(1.*n_steps);
                    (*_obj._p_cav)()     =  p0*(i+1.)/(1.*n_steps);
                }

                libMesh::out
                << "Load step: " << i
                << "  : T = " << (*_obj._temp)()
                << "  : p = " << (*_obj._p_cav)()
                << "  : V = " << (*_obj._velocity)()
                << std::endl;

                _obj._sys->solve();
                _obj._discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
                writer->write_timestep("output_all_steps.exo",
                                       *_obj._eq_sys,
                                       inc+1,
                                       inc+1);
                inc++;
                
            }
            
            // copy the solution to the base solution vector
            sol = *_obj._sys->solution;
            
            _obj._nonlinear_assembly->clear_discipline_and_system();
            
            return sol;
        }
        
        
        /*!
         *   sets the number of steps to be used for nonlinaer steady analysis.
         *   The default is 25 for noninear and 1 for linear.
         */
        void set_n_load_steps(unsigned int n) {
            _n_steps = n;
        }
        
        
        /*!
         *   Tells the solver to increment only the aero loads during the
         *   nonlinear stepping. \p false by default.
         */
        void set_modify_only_aero_load(bool f) {
            _if_only_aero_load_steps = f;
        }
        
        
        /*!
         * @returns  a non-const-reference to the solution.
         */
        virtual libMesh::NumericVector<Real>&
        solution() { return _obj._sys->get_vector("base_solution"); }

        
        
        /*!
         * @returns  a const-reference to the solution.
         */
        virtual const libMesh::NumericVector<Real>&
        solution() const { return _obj._sys->get_vector("base_solution"); }
        
        
    protected:
        
        /*!
         *   pointer to the object that hold all the solution data
         */
        MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization& _obj;
        
        /*!
         *   flag to toggle output
         */
        bool _if_write_output;
        
        /*!
         *   number of nonliear load increment steps
         */
        unsigned int _n_steps;
        
        /*!
         *   only aero load is modified during nonlinear iterations
         */
        bool _if_only_aero_load_steps;
        
        
        /*!
         *   deletes the solution vector from system when the class is
         *   destructed unless this flag is false.
         */
        bool _if_clear_vector_on_exit;
        libMesh::ExodusII_IO* writer;
        unsigned int inc;
    };
}





void
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
evaluate(const std::vector<Real>& dvars,
         Real& obj,
         bool eval_obj_grad,
         std::vector<Real>& obj_grad,
         std::vector<Real>& fvals,
         std::vector<bool>& eval_grads,
         std::vector<Real>& grads) {
    
    libmesh_assert(_initialized);
    libmesh_assert_equal_to(dvars.size(), _n_vars);
    
    // while this function tries to ensure that the gradient and function values
    // returned to the optimizer on all ranks is the same, we may run the
    // possibility of slight differences in the input DVs due to unforeseen
    // reasons. Hence, we will borrow the DV values from rank 0, just to
    // make sure that all processors are evaluating the functions at the
    // same exact values.
    std::vector<Real>
    my_dvars(dvars.begin(), dvars.end());
    __init->comm().broadcast(my_dvars);
    
    
    // set the parameter values equal to the DV value
    // first the plate thickness values
    for (unsigned int i=0; i<_n_vars; i++)
        (*_problem_parameters[i]) = my_dvars[i]*_dv_scaling[i];
    
    
    // DO NOT zero out the gradient vector, since GCMMA needs it for the
    // subproblem solution
    // zero the function evaluations
    std::fill(fvals.begin(), fvals.end(), 0.);
    

    
    libMesh::Point pt; // dummy point object
    
    libMesh::out << "New Eval" << std::endl;
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "th     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << (*_problem_parameters[i])() << std::endl;
    
    
    // the optimization problem is defined as
    // min weight, subject to constraints on displacement and stresses
    Real
    wt      = 0.,
    pval    = 2.;
    
    bool
    if_write_output = true;
    
    // calculate weight
    (*_weight)(pt, 0., wt);
    
    
    //////////////////////////////////////////////////////////////////////
    // first zero the solution
    //////////////////////////////////////////////////////////////////////
    _sys->solution->zero();
    this->clear_stresss();
    // solve for the steady state at zero velocity
    (*_velocity) = 0.;
    
    
    //////////////////////////////////////////////////////////////////////
    // steady state solution
    //////////////////////////////////////////////////////////////////////
    MAST::StiffenedPlateSteadySolverInterface steady_solve(*this,
                                                           if_write_output,
                                                           false,
                                                           _n_load_steps);

    steady_solve.solve();
    // use this solution as the base solution later if no flutter is found.
    libMesh::NumericVector<Real>&
    steady_sol_wo_aero = _sys->add_vector("steady_sol_wo_aero");
    steady_sol_wo_aero = *_sys->solution;
    
    
    // now that we have the solution under the influence of thermal loads,
    // we will use a small number of load steps to get to the steady-state
    // solution
    steady_solve.set_n_load_steps(50);
    steady_solve.set_modify_only_aero_load(true);
    
    
    //////////////////////////////////////////////////////////////////////
    // perform the modal and flutter analysis
    //////////////////////////////////////////////////////////////////////
    // modal analysis is about the base state exclusing aerodynamic loads.
    // So, they will act as generalized coordinates that will not provide
    // diagonal reduced order mass/stiffness operator. The eigenvalues
    // will also be independent of velocity.
    _modal_assembly->attach_discipline_and_system(*_discipline, *_structural_sys);
    _modal_assembly->set_base_solution(steady_sol_wo_aero);
    _sys->eigenproblem_solve();
    _modal_assembly->clear_discipline_and_system();
    
    unsigned int
    nconv = std::min(_sys->get_n_converged_eigenvalues(),
                     _sys->get_n_requested_eigenvalues());
    if (_basis.size() > 0)
        libmesh_assert(_basis.size() == nconv);
    else {
        _basis.resize(nconv);
        for (unsigned int i=0; i<_basis.size(); i++)
            _basis[i] = NULL;
    }
    
    
    // vector of eigenvalues
    std::vector<Real> eig_vals(nconv);
    
    bool if_all_eig_positive = true;
    
    libMesh::ExodusII_IO*
    writer = nullptr;
    
    if (if_write_output)
        writer = new libMesh::ExodusII_IO(*_mesh);

    for (unsigned int i=0; i<nconv; i++) {
        
        // create a vector to store the basis
        if (_basis[i] == NULL)
            _basis[i] = _sys->solution->zero_clone().release();
        
        // now write the eigenvalue
        Real
        re = 0.,
        im = 0.;
        _sys->get_eigenpair(i, re, im, *_basis[i]);
        
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << re << std::endl;
        
        eig_vals[i]          = re;
        if_all_eig_positive  = (if_all_eig_positive && (re>0.))?true:false;
        
        if (if_write_output) {
            
            // copy the solution for output
            (*_sys->solution) = *_basis[i];
            
            // We write the file in the ExodusII format.
            std::set<std::string> nm;
            nm.insert(_sys->name());
            writer->write_timestep("modes.exo",
                                   *_eq_sys,
                                   i+1, i);
        }
    }
    
    
    //////////////////////////////////////////////////////////////////////
    // perform the flutter analysis
    //////////////////////////////////////////////////////////////////////
    std::pair<bool, MAST::FlutterRootBase*> sol(false, nullptr);
    
    // perform flutter analysis only if all modal eigenvalues are positive
    if (if_all_eig_positive) {
        
        // the flutter solver calculates the stability eigenvalues about
        // an equilibrium state that includes the aerodynamic loads. The
        // stability analysis starts with V=0, which corresponds to the
        // equilibrium state calculated for the thermal loads. However,
        // a new equilibrium state is calculated for each velocity. So,
        // the flutter root, if found will depend on the equilibrium state
        // which, in turn, depends on the velocity.
        _fsi_assembly->attach_discipline_and_system(*_discipline,
                                                    *_structural_sys);
        _fsi_assembly->set_base_solution(steady_solve.solution());
        _flutter_solver->clear_solutions();
        _flutter_solver->attach_assembly(*_fsi_assembly);
        _flutter_solver->attach_steady_solver(steady_solve);
        _flutter_solver->initialize(*_velocity,
                                    0.0e3,                // lower V
                                    2*_V0_flutter,        // upper V
                                    _n_V_divs_flutter,    // number of divisions
                                    _basis);              // basis vectors
        
        sol = _flutter_solver->analyze_and_find_critical_root_without_tracking(1.e-3, 20);
        _flutter_solver->print_sorted_roots();
        _fsi_assembly->clear_discipline_and_system();
        _flutter_solver->clear_assembly_object();
    }
    else
        libMesh::out
        << "** Negative frequency: skipping flutter analysis **"
        << std::endl;
    
    
    // if the flutter root was not found, then we will use the steady sol wo
    // aero as the base solution for all the following computations
    // Otherwise, the equilibrium state is dependent on the velocity.
    if (!sol.second) {
        
        // velocity should be set to zero for all residual calculations
        (*_velocity) = 0.;
        steady_solve.solution() = steady_sol_wo_aero;
        *_sys->solution         = steady_sol_wo_aero;
    }
    
    
    
    // now calculate the stress output based on the velocity output
    _nonlinear_assembly->attach_discipline_and_system(*_discipline, *_structural_sys);
    _nonlinear_assembly->calculate_outputs(steady_solve.solution());
    _nonlinear_assembly->clear_discipline_and_system();
    
    
    if (if_write_output) {
        
        libMesh::out << "Writing output to : output.exo" << std::endl;
        
        std::set<std::string> nm;
        nm.insert(_sys->name());
        // write the solution for visualization
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                            *_eq_sys,
                                                            &nm);
        
        _discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
    }
    
    
    
    
    if (sol.second && if_write_output) {
        // now write the flutter mode to an output file.
        // Flutter mode Y = sum_i (X_i * (xi_re + xi_im)_i)
        // using the right eigenvector of the system.
        // where i is the structural mode
        //
        // The time domain simulation assumes the temporal solution to be
        // X(t) = (Y_re + i Y_im) exp(p t)
        //      = (Y_re + i Y_im) exp(p_re t) * (cos(p_im t) + i sin(p_im t))
        //      = exp(p_re t) (Z_re + i Z_im ),
        // where Z_re = Y_re cos(p_im t) - Y_im sin(p_im t), and
        //       Z_im = Y_re sin(p_im t) + Y_im cos(p_im t).
        //
        // We write the simulation of the mode over a period of oscillation
        //
        
        
        // first calculate the real and imaginary vectors
        std::auto_ptr<libMesh::NumericVector<Real> >
        re(_sys->solution->zero_clone().release()),
        im(_sys->solution->zero_clone().release());
        
        
        // first the real part
        _sys->solution->zero();
        for (unsigned int i=0; i<_basis.size(); i++) {
            re->add(sol.second->eig_vec_right(i).real(), *_basis[i]);
            im->add(sol.second->eig_vec_right(i).imag(), *_basis[i]);
        }
        re->close();
        im->close();
        
        // now open the output processor for writing
        libMesh::ExodusII_IO flutter_mode_output(*_mesh);
        
        // use N steps in a time-period
        Real
        t_sys = _sys->time,
        pi    = acos(-1.);
        unsigned int
        N_divs = 100;
        
        
        for (unsigned int i=0; i<=N_divs; i++) {
            _sys->time   =  2.*pi*(i*1.)/(N_divs*1.);
            
            _sys->solution->zero();
            _sys->solution->add( cos(_sys->time), *re);
            _sys->solution->add(-sin(_sys->time), *im);
            _sys->solution->close();
            flutter_mode_output.write_timestep("flutter_mode.exo",
                                               *_eq_sys,
                                               i+1,
                                               _sys->time);
        }
        
        // reset the system time
        _sys->time = t_sys;
    }
    
    //////////////////////////////////////////////////////////////////////
    // get the objective and constraints
    //////////////////////////////////////////////////////////////////////
    
    // set the function and objective values
    obj = wt;
    
    // parallel sum of the weight
    this->comm().sum(obj);
    
    
    // set the eigenvalue constraints  -eig <= 0. scale
    // by an arbitrary 1/1.e7 factor
    for (unsigned int i=0; i<nconv; i++)
        fvals[i] = -eig_vals[i]/1.e7;
    
    
    // we need to make sure that the stress functional in a parallel environment
    // correspond to unique spots in the constraint function vector. Below,
    // the output functionals are created for each local element. These will
    // be mapped to unique spots by identifying the number of elements on each
    // subdomain
    std::vector<unsigned int>
    beginning_elem_id(__init->comm().size(), 0);
    for (unsigned int i=1; i<__init->comm().size(); i++)
        beginning_elem_id[i] = _mesh->n_elem_on_proc(i-1);
    
    // now use this info to identify the beginning elem id
    for (unsigned int i=2; i<__init->comm().size(); i++)
        beginning_elem_id[i] += beginning_elem_id[i-1];

    const unsigned int
    my_id0 = beginning_elem_id[__init->comm().rank()];
    
    // copy the element von Mises stress values as the functions
    for (unsigned int i=0; i<_outputs.size(); i++)
        fvals[_n_eig+1+i+my_id0] =  -1. +
        _outputs[i]->von_Mises_p_norm_functional_for_all_elems(pval)/_stress_limit;
    
    // now sum the value of the stress constraints, so that all ranks
    // have the same values
    for (unsigned int i=0; i<_n_elems; i++)
        this->comm().sum(fvals[_n_eig+1+i]);
    
    // copy the flutter velocity to the constraint vector
    //     Vf        >= V0
    // or, V0/Vf     <= 1
    // or, V0/Vf - 1 <= 0
    //
    if (!if_all_eig_positive)
        fvals[_n_eig+0]  =  100.;
    else if (sol.second)
        fvals[_n_eig+0]  =  _V0_flutter/sol.second->V - 1.;
    else
        fvals[_n_eig+0]  =  -100.;
    
    //////////////////////////////////////////////////////////////////
    //   evaluate sensitivity if needed
    //////////////////////////////////////////////////////////////////
    
    // sensitivity of the objective function
    if (eval_obj_grad) {
        
        Real w_sens = 0.;
        
        // set gradient of weight
        for (unsigned int i=0; i<_n_vars; i++) {
            
            _weight->derivative(MAST::PARTIAL_DERIVATIVE,
                                *_problem_parameters[i],
                                pt,
                                0.,
                                w_sens);
            obj_grad[i] = w_sens*_dv_scaling[i];
        }
        
        // parallel sum
        this->comm().sum(obj_grad);
    }
    
    
    // now check if the sensitivity of objective function is requested
    bool if_sens = false;
    
    for (unsigned int i=0; i<eval_grads.size(); i++)
        if_sens = (if_sens || eval_grads[i]);
    
    if (if_sens) {
        
        // first initialize the gradient vector to zero
        std::fill(grads.begin(), grads.end(), 0.);
        
        //////////////////////////////////////////////////////////////////
        // indices used by GCMMA follow this rule:
        // grad_k = dfi/dxj  ,  where k = j*NFunc + i
        //////////////////////////////////////////////////////////////////

        libMesh::ParameterVector params;
        params.resize(1);
        
        // copy the solution to be used for sensitivity
        // If a flutter solution was found, then this depends on velocity.
        // Otherwise, it is independent of velocity
        *_sys->solution = steady_solve.solution();

        // first do sensitivity analysis wrt velocity, which is necessary for
        // flutter sensitivity.
        libMesh::NumericVector<Real>& dXdV = _sys->add_vector("sol_V_sens");
        std::vector<Real> dsigma_dV(_outputs.size(), 0.);
        
        if (sol.second) {
            
            params[0] = _velocity->ptr();
            _nonlinear_assembly->attach_discipline_and_system(*_discipline,
                                                              *_structural_sys);
            _sys->sensitivity_solve(params);
            dXdV = _sys->get_sensitivity_solution(0);
            
            // This accounts for the dependence of equilibrium state on velocity
            // dsigma/dp   =  par sigma/par p + par sigma/ par V dV/dp
            this->clear_stresss();
            _nonlinear_assembly->calculate_output_sensitivity(params,
                                                              true,
                                                              *(_sys->solution));
            for (unsigned int j=0; j<_outputs.size(); j++)
                dsigma_dV[j] = _outputs[j]->von_Mises_p_norm_functional_sensitivity_for_all_elems
                (pval, _velocity);
            
            _nonlinear_assembly->clear_discipline_and_system();
        }
        else
            dXdV.zero();
        
        
        // we are going to choose to use one parametric sensitivity at a time
        for (unsigned int i=0; i<_n_vars; i++) {
            
            params[0]  = _problem_parameters[i]->ptr();
            
            // copy the solution to be used for sensitivity
            // If a flutter solution was found, then this depends on velocity.
            // Otherwise, it is independent of velocity
            *_sys->solution = steady_solve.solution();

            // iterate over each dv and calculate the sensitivity
            libMesh::NumericVector<Real>& dXdp = _sys->add_sensitivity_solution(0);
            dXdp.zero();
            this->clear_stresss();
            
            // sensitivity analysis
            _nonlinear_assembly->attach_discipline_and_system(*_discipline,
                                                              *_structural_sys);
            _sys->sensitivity_solve(params);
            
            // evaluate sensitivity of the outputs
            _nonlinear_assembly->calculate_output_sensitivity(params,
                                                              true,    // true for total sensitivity
                                                              *(_sys->solution));
            
            _nonlinear_assembly->clear_discipline_and_system();

            // copy the sensitivity values in the output. This accounts for the
            // sensitivity of state wrt parameter. However, if a flutter root
            // was found, the state depends on velocity, which depends on the
            // parameter. Hence, the total sensitivity of stress constraint
            // would need to include the latter component, which was added
            // above.
            for (unsigned int j=0; j<_outputs.size(); j++)
                grads[(i*_n_ineq) + (j+_n_eig+1+my_id0)] = _dv_scaling[i]/_stress_limit *
                _outputs[j]->von_Mises_p_norm_functional_sensitivity_for_all_elems
                (pval, _problem_parameters[i]);
            
            // if all eigenvalues are positive, calculate at the sensitivity of
            // flutter velocity
            if (if_all_eig_positive && sol.second) {
                //
                //           g = V0/Vf - 1 <= 0
                //   Hence, sensitivity is
                //   -V0/Vf^2  dVf
                //
                _fsi_assembly->set_base_solution(steady_solve.solution());
                _fsi_assembly->set_base_solution(dXdp, true);
                _fsi_assembly->attach_discipline_and_system(*_discipline, *_structural_sys);
                _flutter_solver->attach_assembly(*_fsi_assembly);
                _flutter_solver->calculate_sensitivity(*sol.second,
                                                       params,
                                                       0,
                                                       &dXdp,
                                                       &dXdV);
                _fsi_assembly->clear_discipline_and_system();
                _flutter_solver->clear_assembly_object();
                
                grads[(i*_n_ineq) + (_n_eig+0)]  =  -(_dv_scaling[i] *
                                                      _V0_flutter/pow(sol.second->V, 2) *
                                                      sol.second->V_sens);
                
                // add the partial derivative of stress due to velocity sensitivity
                for (unsigned int j=0; j<_outputs.size(); j++)
                    grads[(i*_n_ineq) + (j+_n_eig+1+my_id0)] += _dv_scaling[i]/_stress_limit *
                    dsigma_dV[j] * sol.second->V_sens;
                
                // The natural frequencies were calculated about a base solution
                // that did not include the aerodynamic loads. Hence, we will
                // use a different base soltuion here that what was used for
                // flutter solution. However, a new sensitivity solution is required
                // for this state, since the previous dXdp was calculated for
                // X including aero loads.
                dXdp.zero();
                
                // sensitivity analysis
                (*_velocity) = 0.;
                *_sys->solution = steady_sol_wo_aero;
                _nonlinear_assembly->attach_discipline_and_system(*_discipline,
                                                                  *_structural_sys);
                _sys->sensitivity_solve(params);
                
                _nonlinear_assembly->clear_discipline_and_system();
            }
            else
                // if no root was found, then set the sensitivity to a zero value
                grads[(i*_n_ineq) + (_n_eig+0)]  =  0.;
            
            
            // now, sum the sensitivity of the stress function gradients
            // so that all processors have the same values
            for (unsigned int j=0; j<_n_elems; j++)
                __init->comm().sum(grads[(i*_n_ineq) + (j+_n_eig+1)]);
            
            
            // calculate the sensitivity of the eigenvalues
            std::vector<Real> eig_sens(nconv);
            _modal_assembly->set_base_solution(steady_sol_wo_aero);
            _modal_assembly->set_base_solution(dXdp, true);
            _modal_assembly->attach_discipline_and_system(*_discipline, *_structural_sys);
            // this should not be necessary, but currently the eigenproblem sensitivity
            // depends on availability of matrices before sensitivity
            _sys->assemble_eigensystem();
            _sys->eigenproblem_sensitivity_solve(params, eig_sens);
            _modal_assembly->clear_discipline_and_system();
            
            for (unsigned int j=0; j<nconv; j++)
                grads[(i*_n_ineq) + j] = -_dv_scaling[i]*eig_sens[j]/1.e7;
        }
    }
}





void
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::clear_stresss() {
    
    // iterate over the output quantities and delete them
    std::vector<MAST::StressStrainOutputBase*>::iterator
    it   =   _outputs.begin(),
    end  =   _outputs.end();
    
    for ( ; it != end; it++)
        (*it)->clear(false);
}





void
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
output(unsigned int iter,
       const std::vector<Real>& x,
       Real obj,
       const std::vector<Real>& fval,
       bool if_write_to_optim_file) const {
    
    libmesh_assert_equal_to(x.size(), _n_vars);
    
    // write the DVs in the physical dimension
    for (unsigned int i=0; i<_n_vars; i++)
        libMesh::out
        << "th     [ " << std::setw(10) << i << " ] = "
        << std::setw(20) << (*_problem_parameters[i])() << std::endl;
    
    
    // write the solution for visualization
    std::set<std::string> nm;
    nm.insert(_sys->name());
    libMesh::ExodusII_IO(*_mesh).write_equation_systems("output.exo",
                                                        *_eq_sys,
                                                        &nm);
    _discipline->plot_stress_strain_data<libMesh::ExodusII_IO>("stress_output.exo");
    
    MAST::FunctionEvaluation::output(iter, x, obj, fval, if_write_to_optim_file);
}




MAST::FunctionEvaluation::funobj
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
get_objective_evaluation_function() {
    
    return stiffened_plate_thermally_stressed_piston_theory_flutter_optim_obj;
}



MAST::FunctionEvaluation::funcon
MAST::StiffenedPlateThermallyStressedPistonTheorySizingOptimization::
get_constraint_evaluation_function() {
    
    return stiffened_plate_thermally_stressed_piston_theory_flutter_optim_con;
}

