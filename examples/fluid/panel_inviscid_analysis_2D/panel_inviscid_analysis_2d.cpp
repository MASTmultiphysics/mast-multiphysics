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


// MAST includes
#include "examples/fluid/panel_inviscid_analysis_2D/panel_inviscid_analysis_2d.h"
#include "examples/fluid/meshing/panel_mesh_2D.h"
#include "examples/base/input_wrapper.h"


// libMesh includes
#include "libmesh/parallel_mesh.h"
#include "libmesh/string_to_enum.h"


MAST::Examples::PanelAnalysis2D::
PanelAnalysis2D(const libMesh::Parallel::Communicator& comm_in):
MAST::Examples::FluidExampleBase(comm_in) {
    
}


MAST::Examples::PanelAnalysis2D::~PanelAnalysis2D() {
    
}


void
MAST::Examples::PanelAnalysis2D::_init_mesh() {

    _mesh              = new libMesh::ParallelMesh(this->comm());
    const unsigned int
    dim                 = 2,
    nx_divs             = 3,
    ny_divs             = 1,
    n_max_bumps_x       = (*_input)(_prefix+"n_max_bumps_x", "number of half sin-waves on panel", 1),
    n_elems_panel       = (*_input)(_prefix+"n_elems_panel", "number of elements on the panel",  10),
    n_elems_ff          = (*_input)(_prefix+"n_elems_le", "number of elements from panel to far-field boundary",  20);

    const bool
    if_cos_bump         = (*_input)("if_cos_bump", "if the panel curvature is defined by a cos function", false);

    const Real
    c                   = (*_input)("chord", "panel chord", 0.3),
    l_by_c              = (*_input)("l_by_c", "far-field distance to panel chord ratio", 5.),
    h_ff_by_h_le        = (*_input)("h_far_field_by_h_le", "relative element size at far-field boundary to element size at leading edge", 20.),
    t_by_c              = (*_input)("t_by_c", "panel thickness-to-chord ratio", 0.05);
    
    std::string
    t = (*_input)(_prefix+"elem_type", "type of geometric element in the mesh", "quad4");
    
    libMesh::ElemType
    e_type = libMesh::Utility::string_to_enum<libMesh::ElemType>(t);
    
    // if high order FE is used, libMesh requires atleast a second order
    // geometric element.
    if (_fetype.order > 1 && e_type == libMesh::QUAD4)
        e_type = libMesh::QUAD9;
    else if (_fetype.order > 1 && e_type == libMesh::TRI3)
        e_type = libMesh::TRI6;

    std::vector<Real>
    x_div_loc        = {-l_by_c*c, 0., c, (1.+l_by_c)*c},
    x_relative_dx    = {h_ff_by_h_le, 1., 1., h_ff_by_h_le},
    y_div_loc        = {0., l_by_c*c},
    y_relative_dx    = {1., h_ff_by_h_le};

    std::vector<unsigned int>
    x_divs           = {n_elems_ff, n_elems_panel, n_elems_ff},
    y_divs           = {n_elems_ff};

    std::unique_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs    (new MeshInitializer::CoordinateDivisions),
    y_coord_divs    (new MeshInitializer::CoordinateDivisions);

    std::vector<MeshInitializer::CoordinateDivisions*>
    divs(dim);


    divs[0] = x_coord_divs.get();
    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);

    // initialize the mesh
    MAST::PanelMesh2D().init(t_by_c,
                             if_cos_bump,
                             n_max_bumps_x,
                             divs,
                             *_mesh,
                             e_type);
}



void
MAST::Examples::PanelAnalysis2D::_init_loads() {
    
    std::vector<unsigned int>
    slip      =  {4},
    no_slip,
    symm      =  {5},
    far_field =  {1, 2, 3};
    _init_boundary_conditions(slip, no_slip, symm, far_field);
}


//
//const libMesh::NumericVector<Real>&
//MAST::PanelInviscidAnalysis2D::solve(bool if_write_output) {
//    
//    // initialize the solution
//    RealVectorX s = RealVectorX::Zero(4);
//    s(0) = _flight_cond->rho();
//    s(1) = _flight_cond->rho_u1();
//    s(2) = _flight_cond->rho_u2();
//    s(3) = _flight_cond->rho_e();
//    _fluid_sys->initialize_solution(s);
//    
//    // create the nonlinear assembly object
//    MAST::TransientAssembly                                  assembly;
//    MAST::ConservativeFluidTransientAssemblyElemOperations   elem_ops;
//
//    // Transient solver for time integration
//    MAST::FirstOrderNewmarkTransientSolver  solver;
//    
//    // now solve the system
//    assembly.set_discipline_and_system(*_discipline,
//                                       *_fluid_sys);
//
//    MAST::NonlinearSystem&  nonlin_sys = _fluid_sys->system();
//
//    
//    // file to write the solution for visualization
//    libMesh::ExodusII_IO exodus_writer(*_mesh);
//    
//    // time solver parameters
//    unsigned int
//    t_step            = 0,
//    n_iters_change_dt = 4,
//    iter_count_dt     = 0;
//    
//    Real
//    tval       = 0.,
//    vel_0      = 0.,
//    vel_1      = 1.e+12,
//    p          = 0.5,
//    factor     = 0.,
//    min_factor = 1.5;
//    
//    solver.dt            = _time_step_size;
//    solver.beta          = 1.0;
//    
//    // set the previous state to be same as the current state to account for
//    // zero velocity as the initial condition
//    solver.solution(1).zero();
//    solver.solution(1).add(1., solver.solution());
//    solver.solution(1).close();
//    
//    
//    if (if_write_output)
//        libMesh::out << "Writing output to : output.exo" << std::endl;
//
//    // loop over time steps
//    while ((t_step <= _max_time_steps) &&
//           (vel_1  >=  1.e-8)) {
//        
//        // change dt if the iteration count has increased to threshold
//        if (iter_count_dt == n_iters_change_dt) {
//            
//            libMesh::out
//            << "Changing dt:  old dt = " << solver.dt
//            << "    new dt = ";
//            
//            factor        = std::pow(vel_0/vel_1, p);
//            factor        = std::max(factor, min_factor);
//            solver.dt    *= factor;
//            
//            libMesh::out << solver.dt << std::endl;
//            
//            iter_count_dt = 0;
//            vel_0         = vel_1;
//        }
//        else
//            iter_count_dt++;
//        
//        libMesh::out
//        << "Time step: " << t_step
//        << " :  t = " << tval
//        << " :  xdot-L2 = " << vel_1
//        << std::endl;
//
//        // write the time-step
//        if (if_write_output) {
//
//            exodus_writer.write_timestep("output.exo",
//                                         *_eq_sys,
//                                         t_step+1,
//                                         nonlin_sys.time);
//        }
//        
//        solver.solve();
//        
//        solver.advance_time_step();
//
//        // get the velocity L2 norm
//        vel_1 = solver.velocity().l2_norm();
//        if (t_step == 0) vel_0 = vel_1;
//
//        tval  += solver.dt;
//        t_step++;
//    }
//    
//    assembly.clear_discipline_and_system();
//
//    return *(_sys->solution);
//}
