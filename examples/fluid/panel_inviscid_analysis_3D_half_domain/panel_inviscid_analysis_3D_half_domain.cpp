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
#include <iostream>



// MAST includes
#include "examples/fluid/panel_inviscid_analysis_3D_half_domain/panel_inviscid_analysis_3D_half_domain.h"
#include "examples/fluid/meshing/panel_mesh_3D_half_domain.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "base/transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/elem.h"



extern libMesh::LibMeshInit* __init;



MAST::PanelInviscidAnalysis3DHalfDomain::PanelInviscidAnalysis3DHalfDomain() {


    // initialize the libMesh object
    _mesh              = new libMesh::ParallelMesh(__init->comm());
    _eq_sys            = new libMesh::EquationSystems(*_mesh);
    
    // add the system to be used for analysis
    _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    
    
    // initialize the flow conditions
    GetPot infile("input.in");

    
    const unsigned int
    dim                 = 3,
    nx_divs             = 3,
    ny_divs             = 2,
    nz_divs             = 1,
    panel_bc_id         = 10,
    symmetry_bc_id      = 11,
    n_max_bumps_x       = infile("n_max_bumps_x",    1),
    n_max_bumps_y       = infile("n_max_bumps_y",    1);
    
    const bool
    if_cos_bump         = infile("if_cos_bump",  false);

    const Real
    t_by_c              = infile("t_by_c",         0.05);
    
    libMesh::ElemType
    elem_type           =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(infile("elem_type", "HEX8"));

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
    y_relative_dx    (ny_divs+1),
    z_div_loc        (nz_divs+1),
    z_relative_dx    (nz_divs+1);
    
    std::vector<unsigned int>
    x_divs           (nx_divs),
    y_divs           (ny_divs),
    z_divs           (nz_divs);
    
    std::unique_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MeshInitializer::CoordinateDivisions),
    z_coord_divs    (new MeshInitializer::CoordinateDivisions);
    
    std::vector<MeshInitializer::CoordinateDivisions*>
    divs(dim);
    
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<nx_divs+1; i_div++) {
        
        x_div_loc[i_div]        = infile("x_div_loc",   0., i_div);
        x_relative_dx[i_div]    = infile( "x_rel_dx",   0., i_div);
        
        if (i_div < nx_divs) //  this is only till nx_divs
            x_divs[i_div]       = infile( "x_div_nelem", 0, i_div);
    }
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    
    
    // now read in the values: y-coord
    const Real
    y0  = infile("y_div_loc", 0., 1),
    y1  = infile("y_div_loc", 0., 2),
    dy0 = infile("y_rel_dx",  0., 1),
    dy1 = infile("y_rel_dx",  0., 2);
    
    y_div_loc[0]         = 0.5*(y0+y1);
    y_div_loc[1]         = y1;
    y_div_loc[2]         = infile("y_div_loc", 0., 3);
    
    y_relative_dx[0] = .5*(dy0+dy1);
    y_relative_dx[1] = dy1;
    y_relative_dx[2] = infile( "y_rel_dx", 0., 3);
    
    // halve the num y-divs for panel
    y_divs[0]    = infile( "y_div_nelem",  0, 1)/2;
    y_divs[1]    = infile( "y_div_nelem",  0, 2);
    
    
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    
    // now read in the values: z-coord
    for (unsigned int i_div=0; i_div<nz_divs+1; i_div++) {
        
        z_div_loc[i_div]     = infile("z_div_loc", 0., i_div);
        z_relative_dx[i_div] = infile( "z_rel_dx", 0., i_div);
        
        if (i_div < nz_divs) //  this is only till ny_divs
            z_divs[i_div]    = infile( "z_div_nelem",  0, i_div);
    }
    
    divs[2] = z_coord_divs.get();
    z_coord_divs->init(nz_divs, z_div_loc, z_relative_dx, z_divs);
    
    
    // initialize the mesh
    MAST::PanelMesh3DHalfDomain().init(t_by_c,               // t/c
                                       if_cos_bump,          // if cos bump
                                       n_max_bumps_x,        // n max bumps in x
                                       n_max_bumps_y,        // n max bumps in y
                                       panel_bc_id,
                                       symmetry_bc_id,
                                       divs,
                                       *_mesh,
                                       elem_type);
    
    
    
    _mesh->print_info();
    
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
    _far_field    = new MAST::BoundaryConditionBase(    MAST::FAR_FIELD),
    _slip_wall    = new MAST::BoundaryConditionBase(    MAST::SLIP_WALL);
    _symm_wall    = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    
    
    // tell the physics about these conditions
    _discipline->add_side_load(              1, *_symm_wall);
    _discipline->add_side_load(    panel_bc_id, *_slip_wall);
    _discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    //if (dim == 2)
    for (unsigned int i=2; i<=5; i++)
        _discipline->add_side_load(              i, *_far_field);
    
    // time step control
    _max_time_steps    =   infile("max_time_steps", 1000);
    _time_step_size    =   infile("initial_dt",    1.e-2);
    
    
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
}






MAST::PanelInviscidAnalysis3DHalfDomain::~PanelInviscidAnalysis3DHalfDomain() {
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _fluid_sys;
    
    delete _far_field;
    delete _slip_wall;
    delete _symm_wall;
    
    delete _flight_cond;
}



const libMesh::NumericVector<Real>&
MAST::PanelInviscidAnalysis3DHalfDomain::solve(bool if_write_output) {
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(5);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_u3();
    s(4) = _flight_cond->rho_e();
    _fluid_sys->initialize_solution(s);
    
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                                 assembly;
    MAST::ConservativeFluidTransientAssemblyElemOperations  elem_ops;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(elem_ops,
                                          *_discipline,
                                          solver,
                                          *_fluid_sys);
    
    MAST::NonlinearSystem& nonlin_sys = _fluid_sys->system();

    
    // file to write the solution for visualization
    libMesh::Nemesis_IO exodus_writer(*_mesh);
    
    // time solver parameters
    unsigned int
    t_step            = 0,
    n_iters_change_dt = 4,
    iter_count_dt     = 0;
    
    Real
    tval       = 0.,
    vel_0      = 0.,
    vel_1      = 1.e+12,
    p          = 0.5,
    factor     = 0.,
    min_factor = 1.5;
    
    solver.dt            = _time_step_size;
    solver.beta          = 1.0;
    
    // set the previous state to be same as the current state to account for
    // zero velocity as the initial condition
    solver.solution(1).zero();
    solver.solution(1).add(1., solver.solution());
    solver.solution(1).close();
    
    if (if_write_output)
        libMesh::out << "Writing output to : output.exo" << std::endl;

    // loop over time steps
    while ((t_step <= _max_time_steps) &&
           (vel_1  >=  1.e-8)) {
        
        // change dt if the iteration count has increased to threshold
        if (iter_count_dt == n_iters_change_dt) {
            
            libMesh::out
            << "Changing dt:  old dt = " << solver.dt
            << "    new dt = ";
            
            factor        = std::pow(vel_0/vel_1, p);
            factor        = std::max(factor, min_factor);
            solver.dt    *= factor;
            
            libMesh::out << solver.dt << std::endl;
            
            iter_count_dt = 0;
            vel_0         = vel_1;
        }
        else
            iter_count_dt++;
        
        libMesh::out
        << "Time step: " << t_step
        << " :  t = " << tval
        << " :  xdot-L2 = " << vel_1
        << std::endl;
        
        // write the time-step
        if (if_write_output) {
            
            exodus_writer.write_timestep("output.exo",
                                         *_eq_sys,
                                         t_step+1,
                                         nonlin_sys.time);
        }
        
        solver.solve();
        
        solver.advance_time_step();
        
        // get the velocity L2 norm
        vel_1 = solver.velocity().l2_norm();
        if (t_step == 0) vel_0 = vel_1;
        
        tval  += solver.dt;
        t_step++;
    }
    
    assembly.clear_discipline_and_system();

    return *(_sys->solution);
}





const libMesh::NumericVector<Real>&
MAST::PanelInviscidAnalysis3DHalfDomain::
sensitivity_solve(MAST::Parameter& p,
                  bool if_write_output) {
    
    libmesh_assert(false); // to be implemented

    return _sys->get_sensitivity_solution(0);
}

