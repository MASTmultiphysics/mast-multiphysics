
// C++ includes
#include <iostream>


// MAST includes
#include "heat_conduction/heat_conduction_discipline.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_transient_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "property_cards/element_property_card_3D.h"
#include "driver/driver_base.h"
#include "solver/first_order_newmark_transient_solver.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"


int main(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::System& sys = eq_sys.add_system<libMesh::NonlinearImplicitSystem>("thermal");
    
    // initialize the mesh
    unsigned int
    nx        = 5,
    ny        = 5,
    nz        = 5;
    const Real
    width     = 10.,
    height    = 10.,
    depth     = 10.;
    libMesh::MeshTools::Generation::build_cube(mesh,
                                               nx, ny, nz,
                                               0., width,
                                               0., height,
                                               0., depth);
    mesh.prepare_for_use();
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    MAST::HeatConductionDiscipline            heat_cond(eq_sys);
    MAST::HeatConductionSystemInitialization  heat_cond_sys(sys,
                                                            sys.name(),
                                                            fe_type);
    
    
    // add the Dirichlet boundary condition on left boundary
    MAST::DirichletBoundaryCondition   dirichlet;
    dirichlet.init(3, heat_cond_sys.vars());
    
    // add the boundary condition to the physics
    heat_cond.add_dirichlet_bc(3, dirichlet);
    
    // tell the system about the constraints. This needs to happen
    // before the equation system initialization
    heat_cond.init_system_dirichlet_bc(sys);
    
    // initialize the equation system for analysis
    eq_sys.init();
    
    // print the information
    mesh.print_info();
    eq_sys.print_info();
    
    // add a flux load on the right boundary
    // this parameter defines the constant value of the flux
    MAST::Parameter q ("q", 5.0);
    
    // define the field function for boudnary condition
    MAST::ConstantFieldFunction q_flux("heat_flux", q);
    
    // create a flux boundary condition based on this
    MAST::BoundaryConditionBase flux_load(MAST::HEAT_FLUX);
    flux_load.add(q_flux);
    
    // tell the physics about this load
    heat_cond.add_side_load(0, flux_load);
    
    
    // initialize the material properties. First the parameters, which
    // are used to define the properties to be constants
    MAST::Parameter  k (   "k", 2.0);   // thermal conductance
    MAST::Parameter cp (  "cp", 2.0);   // thermal capacitance
    MAST::Parameter rho ("rho", 2.0);   // material density
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction   k_th( "k_th",   k);
    MAST::ConstantFieldFunction  cp_th(   "cp",  cp);
    MAST::ConstantFieldFunction rho_th(  "rho", rho);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials;
    materials.add(  k_th);
    materials.add( cp_th);
    materials.add(rho_th);
    
    
    // initialize the element section property
    MAST::Parameter h("h", 0.02);
    
    // define a constant field function for thickness
    MAST::ConstantFieldFunction  h_val("h", h);
    
    // add to a section property card
    MAST::ElementPropertyCard3D section_property;
    section_property.add(h_val);
    
    
    // tell the section property to use the material property card
    section_property.set_material(materials);
    
    
    // tell the conduction physics about the section properties
    heat_cond.set_property_for_subdomain(0, section_property);
    
    // create the nonlinear assembly object
    MAST::HeatConductionTransientAssembly   assembly;
    
    // Transient solver for time integration
    MAST::FirstOrderNewmarkTransientSolver  solver;
    
    // now solve the system
    assembly.attach_discipline_and_system(heat_cond, solver, heat_cond_sys);
    
    // file to write the solution for visualization
    libMesh::ExodusII_IO exodus_writer(mesh);
    
    // time solver parameters
    unsigned int t_step  = 0;
    unsigned int n_steps = 100;
    solver.dt            = 1.0e1;
    solver.beta          = 1.;
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

