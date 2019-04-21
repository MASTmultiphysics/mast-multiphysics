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

// C++ includes
#include <iostream>



// MAST includes
#include "tests/fluid/build_conservative_fluid_elem.h"
#include "examples/base/rigid_surface_motion.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/frequency_domain_linearized_complex_assembly.h"
#include "solver/complex_solver_base.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "base/boundary_condition_base.h"
#include "aeroelasticity/frequency_function.h"
#include "fluid/frequency_domain_linearized_conservative_fluid_elem.h"
#include "base/nonlinear_system.h"

// libMesh includes
#include "libmesh/mesh_generation.h"


// MAST includes


extern libMesh::LibMeshInit* __init;



MAST::BuildConservativeFluidElem::BuildConservativeFluidElem() {
    
    
    // initialize the libMesh object
    _mesh              = new libMesh::ParallelMesh(__init->comm());
    _eq_sys            = new libMesh::EquationSystems(*_mesh);
    
    // add the system to be used for analysis
    _sys = &(_eq_sys->add_system<MAST::NonlinearSystem>("fluid"));
    
    // initialize the mesh
    unsigned int
    dim       = 2;

    libMesh::MeshTools::Generation::build_square(*_mesh, 1, 1);
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    _discipline        = new MAST::ConservativeFluidDiscipline(*_eq_sys);
    _fluid_sys         = new MAST::ConservativeFluidSystemInitialization(*_sys,
                                                                         _sys->name(),
                                                                         fe_type,
                                                                         dim);
    
    
    // initialize the equation system for analysis
    _eq_sys->init();
    
    // create the oundary conditions for slip-wall and far-field
    _far_field     = new MAST::BoundaryConditionBase(MAST::FAR_FIELD);
    _slip_wall     = new MAST::BoundaryConditionBase(MAST::SLIP_WALL);
    
    // tell the physics about these conditions
    _discipline->add_side_load(1, *_far_field);
    
    
    
    _flight_cond    =  new MAST::FlightCondition;
    _flight_cond->body_roll_axis(0)     = 1.;
    _flight_cond->body_pitch_axis(2)    = 1.;
    _flight_cond->body_yaw_axis(1)      = 1.;
    _flight_cond->body_euler_angles     = RealVectorX::Zero(3);
    _flight_cond->body_angular_rates    = RealVectorX::Zero(3);
    
    _flight_cond->ref_chord       = 1.;
    _flight_cond->altitude        = 0.;
    _flight_cond->mach            = .5;
    _flight_cond->gas_property.cp = 1003.;
    _flight_cond->gas_property.cv = 716.;
    _flight_cond->gas_property.T  = 300.;
    _flight_cond->gas_property.rho= 1.05;
    
    _flight_cond->init();
    
    // tell the discipline about the fluid values
    _discipline->set_flight_condition(*_flight_cond);
    
    // define parameters
    _omega             = new MAST::Parameter("omega",       0.);
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
                  0.,                              // plunge amplitude
                  1.,                              // pitch amplitude
                  0.);                             // pitch phase lead
    _displacement = new MAST::RigidSurfaceDisplacement(*_motion);
    _normal_rot   = new MAST::RigidSurfaceNormalRotation(*_motion);

    _slip_wall->add(*_displacement);
    _slip_wall->add(*_normal_rot);
    _discipline->add_side_load(0, *_slip_wall);
    
    
    // initialize the solution
    _base_sol    = RealVectorX::Zero(4),
    
    _base_sol(0) = _flight_cond->rho();
    _base_sol(1) = _flight_cond->rho_u1();
    _base_sol(2) = _flight_cond->rho_u2();
    _base_sol(3) = _flight_cond->rho_e();

}






MAST::BuildConservativeFluidElem::~BuildConservativeFluidElem() {
    
    delete _eq_sys;
    delete _mesh;
    
    delete _discipline;
    delete _fluid_sys;
    
    delete _far_field;
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
    delete _displacement;
    delete _normal_rot;
}





