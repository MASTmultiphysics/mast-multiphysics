
// C++ includes
#include <iostream>


// MAST includes
#include "heat_conduction/heat_conduction_discipline.h"
#include "heat_conduction/heat_conduction_system_initialization.h"
#include "heat_conduction/heat_conduction_nonlinear_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "driver/driver_base.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"


int main_7(int argc, const char * argv[]) {

    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::System& sys = eq_sys.add_system<libMesh::NonlinearImplicitSystem>("thermal");
    
    // initialize the mesh
    unsigned int
    nx        = 10,
    ny        = 10;
    const Real
    width     = 10.,
    height    = 10.;
    libMesh::MeshTools::Generation::build_square(mesh,
                                                 nx, ny,
                                                 0., width,
                                                 0., height);
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
    MAST::Parameter k  (  "k", 2.0);   // thermal conductance
    
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction  k_th("k_th", k);

    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials;
    materials.add(k_th);
    
    
    // initialize the element section property
    MAST::Parameter h("h", 0.02);
    
    // define a constant field function for thickness
    MAST::ConstantFieldFunction  h_val("h", h);
    
    // add to a section property card
    MAST::Solid2DSectionElementPropertyCard section_property;
    section_property.add(h_val);
    
    
    // tell the section property to use the material property card
    section_property.set_material(materials);
    
    
    // tell the conduction physics about the section properties
    heat_cond.set_property_for_subdomain(0, section_property);
    /*
    // create the nonlinear assembly object
    MAST::HeatConductionNonlinearAssembly   assembly;
    
    // now solve the system
    MAST::Driver::nonlinear_solution(heat_cond, heat_cond_sys, assembly);
     */
    // write the solution for visualization
    libMesh::ExodusII_IO(mesh).write_equation_systems("mesh.exo", eq_sys);
    
    return 0;
}

