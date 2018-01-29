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
#include "examples/fluid/ramp_laminar_analysis_2D/ramp_viscous_analysis_2d.h"
#include "examples/fluid/meshing/ramp_mesh_2D.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "base/transient_assembly.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "base/nonlinear_system.h"


// libMesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/elem.h"



extern libMesh::LibMeshInit* __init;



MAST::RampLaminarAnalysis2D::RampLaminarAnalysis2D() {


    // initialize the libMesh object
    _mesh              = new libMesh::ParallelMesh(__init->comm());
    _eq_sys            = new libMesh::EquationSystems(*_mesh);
    
    // add the system to be used for analysis
    _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    
    
    // initialize the flow conditions
    GetPot infile("input.in");

    
    const unsigned int
    dim                 = 2,
    nx_divs             = infile("nx_divs",          2),
    ny_divs             = infile("ny_divs",          1),
    ramp_bc_id          = infile("ramp_bc_id",      10),
    symmetry_bc_id      = infile("symmetry_bc_id",  11);
    
    
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
    
    std::unique_ptr<MeshInitializer::CoordinateDivisions>
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
    MAST::RampMesh2D().init(0.,
                            ramp_bc_id,
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
    
    // create and add the boundary condition and loads
    _dirichlet = new MAST::DirichletBoundaryCondition;
    std::vector<unsigned int> constrained_vars(2);
    constrained_vars[0] = _fluid_sys->vars()[1];
    constrained_vars[1] = _fluid_sys->vars()[2];
    _dirichlet->init (ramp_bc_id, constrained_vars);
    _discipline->add_dirichlet_bc(ramp_bc_id, *_dirichlet);
    _discipline->init_system_dirichlet_bc(*_sys);

    
    // initialize the equation system for analysis
    _eq_sys->init();
    
    // print the information
    _eq_sys->print_info();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field    = new MAST::BoundaryConditionBase(    MAST::FAR_FIELD),
    _slip_wall    = new MAST::BoundaryConditionBase( MAST::NO_SLIP_WALL);
    _symm_wall    = new MAST::BoundaryConditionBase(MAST::SYMMETRY_WALL);
    
    
    // tell the physics about these conditions
    _discipline->add_side_load(     ramp_bc_id, *_slip_wall);
    _discipline->add_side_load( symmetry_bc_id, *_symm_wall);
    // all boundaries except the bottom are far-field
    for (unsigned int i=1; i<=3; i++)
        _discipline->add_side_load(              i, *_far_field);
    
        
    // time step control
    _max_time_steps    =   infile("max_time_steps", 1000);
    _time_step_size    =   infile("initial_dt",    1.e-4);
    
    
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
    _flight_cond->gas_property.if_viscous= infile( "if_viscous",   true);
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _discipline->set_flight_condition(*_flight_cond);
}






MAST::RampLaminarAnalysis2D::~RampLaminarAnalysis2D() {
    
    delete _eq_sys;
    delete _mesh;
    delete _dirichlet;
    
    delete _discipline;
    delete _fluid_sys;
    
    delete _far_field;
    delete _slip_wall;
    delete _symm_wall;
    
    delete _flight_cond;
}




const libMesh::NumericVector<Real>&
MAST::RampLaminarAnalysis2D::solve(bool if_write_output) {
    
    // initialize the solution
    RealVectorX s = RealVectorX::Zero(4);
    s(0) = _flight_cond->rho();
    s(1) = _flight_cond->rho_u1();
    s(2) = _flight_cond->rho_u2();
    s(3) = _flight_cond->rho_e();
    _fluid_sys->initialize_solution(s);
    
    // create the nonlinear assembly object
    MAST::TransientAssembly                                  assembly;
    MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(elem_ops,
                                          *_discipline,
                                          solver,
                                          *_fluid_sys);
    
    MAST::NonlinearSystem&  nonlin_sys   =
    dynamic_cast<MAST::NonlinearSystem&>(_fluid_sys->system());

    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(*_mesh);
    
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
MAST::RampLaminarAnalysis2D::sensitivity_solve(MAST::Parameter& p,
                                          bool if_write_output) {
    
    libmesh_assert(false); // to be implemented

    return _sys->get_sensitivity_solution(0);
}

