
// C++ includes
#include <iostream>


// MAST includes
#include "fluid/conservative_fluid_discipline.h"
#include "fluid/conservative_fluid_system_initialization.h"
#include "fluid/conservative_fluid_transient_assembly.h"
#include "base/physics_discipline_base.h"
#include "solver/first_order_newmark_transient_solver.h"
#include "fluid/flight_condition.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/getpot.h"

int main_5(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::System& sys =
    eq_sys.add_system<libMesh::NonlinearImplicitSystem>("fluid");
    
    // initialize the mesh
    unsigned int
    dim       = 2;

    libMesh::GmshIO(mesh).read("/Users/manav/Documents/acads/Projects/gmsh_models/naca0012/naca0012_mesh0.msh");
    mesh.prepare_for_use();
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    MAST::ConservativeFluidDiscipline            fluid_discipline(eq_sys);
    MAST::ConservativeFluidSystemInitialization  fluid_sys(sys,
                                                           sys.name(),
                                                           fe_type,
                                                           dim);
    
    
    // add the Dirichlet boundary condition on left boundary
    MAST::DirichletBoundaryCondition   dirichlet;
    dirichlet.init(0, fluid_sys.vars());
    
    // add the boundary condition to the physics
    fluid_discipline.add_dirichlet_bc(0, dirichlet);
    
    // tell the system about the constraints. This needs to happen
    // before the equation system initialization
    fluid_discipline.init_system_dirichlet_bc(sys);
    
    // initialize the equation system for analysis
    eq_sys.init();
    
    // print the information
    mesh.print_info();
    eq_sys.print_info();
    
    // create the oundary conditions for slip-wall and far-field
    MAST::BoundaryConditionBase
    far_field(MAST::FAR_FIELD),
    slip_wall(MAST::SLIP_WALL);
    
    // tell the physics about these conditions
    fluid_discipline.add_side_load(0, slip_wall);
    fluid_discipline.add_side_load(1, far_field);
    
    
    // initialize the flow conditions
    GetPot infile("infile.txt");
    MAST::FlightCondition flight_cond;
    for (unsigned int i=0; i<3; i++) {
        
        flight_cond.body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
        flight_cond.body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
        flight_cond.body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
        flight_cond.body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
        flight_cond.body_angular_rates(i) = infile("body_angular_rates", 0., i);
    }
    flight_cond.ref_chord       = infile("ref_c",    1.);
    flight_cond.altitude        = infile( "alt",     0.);
    flight_cond.mach            = infile("mach",     .5);
    flight_cond.gas_property.cp = infile(  "cp",  1003.);
    flight_cond.gas_property.cv = infile(  "cv",   716.);
    flight_cond.gas_property.T  = infile("temp",    60.);
    flight_cond.gas_property.rho= infile( "rho",   1.05);
    flight_cond.init();
    
    // tell the discipline about the fluid values
    fluid_discipline.set_flight_condition(flight_cond);
    
    // create the nonlinear assembly object
    MAST::ConservativeFluidTransientAssembly   assembly;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;

    // now solve the system
    assembly.attach_discipline_and_system(fluid_discipline,
                                          solver,
                                          fluid_sys);
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(mesh);
    
    // time solver parameters
    unsigned int t_step  = 0;
    unsigned int n_steps = 100;
    solver.dt            = 1.0e1;
    solver.beta          = 1.0;
    solver.set_initial_condition(250.);
    
    // loop over time steps
    while (t_step < n_steps) {
        
        // write the time-step
        exodus_writer.write_timestep("out.exo",
                                     eq_sys,
                                     t_step+1,
                                     sys.time);
        
        solver.solve();
        
        
        solver.advance_time_step();
        
        t_step++;
    }
    
    assembly.clear_discipline_and_system();
    
    return 0;
}

