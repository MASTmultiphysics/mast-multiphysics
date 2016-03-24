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
#include "examples/fluid/panel_small_disturbance_frequency_domain_analysis_2D/panel_small_disturbance_frequency_domain_analysis_2d.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "base/nonlinear_system.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "solver/complex_solver_base.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "boundary_condition/rigid_surface_motion.h"
#include "aeroelasticity/frequency_function.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"



extern libMesh::LibMeshInit* __init;



MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis::
PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis() {
    
    
    // initialize the libMesh object
    _mesh              = new libMesh::ParallelMesh(__init->comm());
    _eq_sys            = new libMesh::EquationSystems(*_mesh);
    
    // add the system to be used for analysis
    _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    _sys->set_init_B_matrix();
    
    // initialize the flow conditions
    GetPot infile("input.in");
    
    
    const unsigned int
    dim                 = 2,
    nx_divs             = 3,
    ny_divs             = 1,
    panel_bc_id         = 10,
    symmetry_bc_id      = 11;
    
    libMesh::ElemType
    elem_type           =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "QUAD4"));
    
    libMesh::FEFamily
    fe_type             =
    libMesh::Utility::string_to_enum<libMesh::FEFamily>(infile("fe_family", "LAGRANGE"));
    
    libMesh::Order
    fe_order            =
    libMesh::Utility::string_to_enum<libMesh::Order>(infile("fe_order", "FIRST"));
    
    std::vector<Real>
    x_div_loc        (nx_divs+1),
    x_relative_dx    (nx_divs+1),
    y_div_loc        (ny_divs+1),
    y_relative_dx    (ny_divs+1);
    
    std::vector<unsigned int>
    x_divs           (nx_divs),
    y_divs           (ny_divs);
    
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MeshInitializer::CoordinateDivisions);
    
    std::vector<MeshInitializer::CoordinateDivisions*>
    divs(dim);
    
    
    // now read in the values: x-coord
    if (nx_divs > 0) {
        
        for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
            
            x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
            x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
            
            if (i_div < nx_divs) //  this is only till nx_divs
                x_divs[i_div]       = infile( "x_div_nelem", 0, i_div);
        }
        
        divs[0] = x_coord_divs.get();
        x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    }
    
    
    // now read in the values: y-coord
    if ((dim > 1) && (ny_divs > 0)) {
        
        for (unsigned int i_div=0; i_div<ny_divs+1; i_div++) {
            
            y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
            y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
            
            if (i_div < ny_divs) //  this is only till ny_divs
                y_divs[i_div]    = infile( "y_div_nelem",  0, i_div);
        }
        
        divs[1] = y_coord_divs.get();
        y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    }
    
    
    
    
    // initialize the mesh
    MAST::PanelMesh2D().init(0.,               // t/c
                             false,            // if cos bump
                             0,                // n max bumps
                             panel_bc_id,
                             symmetry_bc_id,
                             divs,
                             *_mesh,
                             elem_type);
    
    _discipline        = new MAST::ConservativeFluidDiscipline(*_eq_sys);
    _fluid_sys         = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                         _sys->name(),
                                                                         libMesh::FEType(fe_order, fe_type),
                                                                         dim);
    
    
    // initialize the equation system for analysis
    _eq_sys->init();
    
    // print the information
    _eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD),
    _symm_wall     = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    _slip_wall     = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
    
    _flight_cond    =  new MAST::FlightCondition;
    for (unsigned int i=0; i<3; i++) {
        
        _flight_cond->body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
        _flight_cond->body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
        _flight_cond->body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
        _flight_cond->body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
        _flight_cond->body_angular_rates(i) = infile("body_angular_rates", 0., i);
    }
    
    _flight_cond->ref_chord       = infile("ref_c",    1.);
    _flight_cond->altitude        = infile( "alt",     0.);
    _flight_cond->mach            = infile("mach",     .5);
    _flight_cond->gas_property.cp = infile(  "cp",  1003.);
    _flight_cond->gas_property.cv = infile(  "cv",   716.);
    _flight_cond->gas_property.T  = infile("temp",   300.);
    _flight_cond->gas_property.rho= infile( "rho",   1.05);
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _discipline->set_flight_condition(*_flight_cond);
    
    // define parameters
    _omega             = new MAST::Parameter("omega",     100.);
    _velocity          = new MAST::Parameter("velocity",  _flight_cond->velocity_magnitude);
    _b_ref             = new MAST::Parameter("b_ref",       1.);
    
    
    // now define the constant field functions based on this
    _omega_f           = new MAST::ConstantFieldFunction("omega",       *_omega);
    _velocity_f        = new MAST::ConstantFieldFunction("velocity", *_velocity);
    _b_ref_f           = new MAST::ConstantFieldFunction("b_ref",       *_b_ref);
    
    // initialize the frequency function
    _freq_function     = new MAST::FrequencyFunction("freq",
                                                     *_omega_f,
                                                     *_velocity_f,
                                                     *_b_ref_f);
    
    // initialize the motion object
    _motion            = new MAST::RigidSurfaceMotion;
    _motion->init(*_freq_function,                 // frequency function
                  _flight_cond->body_yaw_axis,     // plunge vector
                  _flight_cond->body_pitch_axis,   // pitch axis
                  RealVectorX::Zero(3),            // hinge location
                  1.,                              // plunge amplitude
                  0.,                              // pitch amplitude
                  0.);                             // pitch phase lead
    
    // tell the physics about boundary conditions
    _slip_wall->add(*_motion);
    _discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _discipline->add_side_load(              i, *_far_field);
}






MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis::
~PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis() {
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _fluid_sys;
    
    delete _far_field;
    delete _symm_wall;
    delete _slip_wall;
    
    delete _flight_cond;
    
    delete _omega;
    delete _velocity;
    delete _b_ref;
    
    delete _omega_f;
    delete _velocity_f;
    delete _b_ref_f;
    
    delete _freq_function;
    
    delete _motion;
}





MAST::Parameter*
MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis::
get_parameter(const std::string &nm) {
    
    MAST::Parameter *rval = NULL;
    
    // look through the vector of parameters to see if the name is available
    std::vector<MAST::Parameter*>::iterator
    it   =  _params_for_sensitivity.begin(),
    end  =  _params_for_sensitivity.end();
    
    bool
    found = false;
    
    for ( ; it != end; it++) {
        
        if (nm == (*it)->name()) {
            rval    = *it;
            found   = true;
        }
    }
    
    // if the param was not found, then print the message
    if (!found) {
        std::cout
        << std::endl
        << "Parameter not found by name: " << nm << std::endl
        << "Valid names are: "
        << std::endl;
        for (it = _params_for_sensitivity.begin(); it != end; it++)
            std::cout << "   " << (*it)->name() << std::endl;
        std::cout << std::endl;
    }
    
    return rval;
}





const libMesh::NumericVector<Real>&
MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis::
solve(bool if_write_output) {
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(4);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_e();
    
    // create the vector for storing the base solution.
    // we will swap this out with the system solution, initialize and
    // then swap it back.
    libMesh::NumericVector<Real>& base_sol =
    _sys->add_vector("fluid_base_solution");
    _sys->solution->swap(base_sol);
    _fluid_sys->initialize_solution(s);
    _sys->solution->swap(base_sol);
    
    // create the nonlinear assembly object
    MAST::FrequencyDomainLinearizedComplexAssembly   assembly;
    
    // complex solver for solution of small-disturbance system of eqs.
    MAST::ComplexSolverBase                          solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(*_discipline,
                                          solver,
                                          *_fluid_sys);
    assembly.set_base_solution(base_sol);
    assembly.set_frequency_function(*_freq_function);
    
    solver.solve_block_matrix();
    
    if (if_write_output) {
        
        libMesh::NumericVector<Real>
        &real_sol = solver.real_solution(),
        &imag_sol = solver.imag_solution();
        
        
        // first, write the real part
        std::cout
        << "Writing real output to : real_output.exo" << std::endl;
        
        _sys->solution->swap(real_sol);
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("real_output.exo",
                                                            *_eq_sys);
        _sys->solution->swap(real_sol);
        
        
        // next, write the imag part
        std::cout
        << "Writing imag output to : imag_output.exo" << std::endl;
        
        _sys->solution->swap(imag_sol);
        libMesh::ExodusII_IO(*_mesh).write_equation_systems("imag_output.exo",
                                                            *_eq_sys);
        _sys->solution->swap(imag_sol);
        
        
        // now write the mode to an output file.
        // mode Y = sum_i (X_i * (xi_re + xi_im)_i)
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
            _sys->solution->add( cos(_sys->time), real_sol);
            _sys->solution->add(-sin(_sys->time), imag_sol);
            _sys->solution->close();
            flutter_mode_output.write_timestep("complex_sol_transient.exo",
                                               *_eq_sys,
                                               i+1,
                                               _sys->time);
        }
        
        // reset the system time
        _sys->time = t_sys;
    }
    
    
    assembly.clear_discipline_and_system();
    _sys->remove_vector("fluid_base_solution");
    
    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::PanelInviscidSmallDisturbanceFrequencyDomain2DAnalysis::
sensitivity_solve(MAST::Parameter& p, bool if_write_output) {
    
    /*_discipline->add_parameter(p);
     
     // create the nonlinear assembly object
     MAST::StructuralNonlinearAssembly   assembly;
     
     assembly.attach_discipline_and_system(*_discipline, *_structural_sys);
     
     libMesh::NonlinearImplicitSystem&      nonlin_sys   =
     dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
     
     libMesh::ParameterVector params;
     params.resize(1);
     params[0]  =  p.ptr();
     
     // zero the solution before solving
     nonlin_sys.add_sensitivity_solution(0).zero();
     this->clear_stresss();
     
     nonlin_sys.sensitivity_solve(params);
     
     // evaluate sensitivity of the outputs
     assembly.calculate_output_sensitivity(params,
     true,    // true for total sensitivity
     *(_sys->solution));
     
     
     assembly.clear_discipline_and_system();
     _discipline->remove_parameter(p);
     
     // write the solution for visualization
     if (if_write_output) {
     
     std::ostringstream oss1, oss2;
     oss1 << "output_" << p.name() << ".exo";
     oss2 << "output_" << p.name() << ".exo";
     
     std::cout
     << "Writing sensitivity output to : " << oss1.str()
     << "  and stress/strain sensitivity to : " << oss2.str()
     << std::endl;
     
     
     _sys->solution->swap(_sys->get_sensitivity_solution(0));
     
     // write the solution for visualization
     libMesh::ExodusII_IO(*_mesh).write_equation_systems(oss1.str(),
     *_eq_sys);
     _discipline->plot_stress_strain_data<libMesh::ExodusII_IO>(oss2.str(), &p);
     
     _sys->solution->swap(_sys->get_sensitivity_solution(0));
     }
     */
    return _sys->get_sensitivity_solution(0);
}

