
// C++ includes
#include <iostream>


// MAST includes
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/physics_discipline_base.h"
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "boundary_condition/dirichlet_boundary_condition.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "driver/driver_base.h"


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
    libMesh::System& sys =
    eq_sys.add_system<libMesh::NonlinearImplicitSystem>("structural");
    
    // initialize the mesh
    unsigned int
    nx        = 10;
    const Real
    length     = 10.;
    libMesh::MeshTools::Generation::build_line(mesh,
                                               nx,
                                               0., length);
    mesh.prepare_for_use();
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    MAST::StructuralDiscipline            structural_discipline(eq_sys);
    MAST::StructuralSystemInitialization  structural_sys(sys,
                                                         sys.name(),
                                                         fe_type);
    
    
    // add the Dirichlet boundary condition on left boundary
    MAST::DirichletBoundaryCondition   dirichlet;
    dirichlet.init(0, structural_sys.vars());
    
    // add the boundary condition to the physics
    structural_discipline.add_dirichlet_bc(0, dirichlet);
    
    // tell the system about the constraints. This needs to happen
    // before the equation system initialization
    structural_discipline.init_system_dirichlet_bc(sys);
    
    // initialize the equation system for analysis
    eq_sys.init();
    
    // print the information
    mesh.print_info();
    eq_sys.print_info();
    
    // add a flux load on the right boundary
    // this parameter defines the constant value of the flux
    MAST::Parameter p ("p", 5.0);
    
    // define the field function for boudnary condition
    MAST::ConstantFieldFunction press("pressure", p);
    
    // create a flux boundary condition based on this
    MAST::BoundaryConditionBase flux_load(MAST::SURFACE_PRESSURE);
    flux_load.add(press);
    
    // tell the physics about this load
    structural_discipline.add_side_load(1, flux_load);
    
    
    // initialize the material properties. First the parameters, which
    // are used to define the properties to be constants
    MAST::Parameter
    E  (  "E",  2.0),   // Young's modulus
    nu (  "nu",0.33);   // Poisson's ratio
    
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction
    E_f ("E",  E),
    nu_f("nu", nu);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials;
    materials.add( E_f);
    materials.add(nu_f);
    
    
    // initialize the element section property
    MAST::Parameter
    hy   ("hy", 0.02),
    hz   ("hz", 0.04),
    zero ("zero", 0.);
    
    // define a constant field function for thickness
    MAST::ConstantFieldFunction
    hy_f    ("hy",       hy),
    hz_f    ("hz",       hz),
    hyoff_f ("hy_off", zero),
    hzoff_f ("hz_off", zero);
    
    // add to a section property card
    MAST::Solid1DSectionElementPropertyCard section_property;
    section_property.add   (hy_f);
    section_property.add   (hz_f);
    section_property.add(hyoff_f);
    section_property.add(hzoff_f);
    section_property.y_vector()(1) = 1.;  // y-vector of the section is along y-axis
    
    // tell the section property to use the material property card
    section_property.set_material(materials);
    section_property.init();
    
    // tell the conduction physics about the section properties
    structural_discipline.set_property_for_subdomain(0, section_property);
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;
    
    // now solve the system
    MAST::Driver::nonlinear_solution(structural_discipline, structural_sys, assembly);
    
    sys.solution->print();
    
    // write the solution for visualization
    libMesh::ExodusII_IO(mesh).write_equation_systems("mesh.exo", eq_sys);
    
    return 0;
}

