
// C++ includes
#include <iostream>


// MAST includes
#include "elasticity/structural_discipline.h"
#include "elasticity/structural_system_initialization.h"
#include "elasticity/structural_modal_eigenproblem_assembly.h"
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
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"


int main(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::CondensedEigenSystem& sys =
    eq_sys.add_system<libMesh::CondensedEigenSystem>("structural");
    sys.set_eigenproblem_type(libMesh::GHEP);
    sys.eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    
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
    MAST::DirichletBoundaryCondition   dirichlet1, dirichlet2;
    // constrain to be a pinned-pinned beam
    std::vector<unsigned int> constrained_vars(4);
    constrained_vars[0] = 0;  // u
    constrained_vars[1] = 1;  // v
    constrained_vars[2] = 2;  // w
    constrained_vars[3] = 3;  // tx
    dirichlet1.init(0, constrained_vars);
    dirichlet2.init(1, constrained_vars);
    
    // add the boundary condition to the physics
    structural_discipline.add_dirichlet_bc(0, dirichlet1);
    structural_discipline.add_dirichlet_bc(1, dirichlet2);
    
    
    // initialize the equation system for analysis
    eq_sys.init();

    // print the information
    mesh.print_info();
    eq_sys.print_info();
    
    // add a flux load on the right boundary
    
    // initialize the material properties. First the parameters, which
    // are used to define the properties to be constants
    MAST::Parameter
    E  (  "E", 72.0e9),   // Young's modulus
    nu (  "nu",  0.33),   // Poisson's ratio
    rho( "rho", 2700.);    // Density
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction
    E_f  (  "E",  E),
    nu_f ( "nu", nu),
    rho_f("rho",rho);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials;
    materials.add(  E_f);
    materials.add( nu_f);
    materials.add(rho_f);
    
    
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
    MAST::StructuralModalEigenproblemAssembly   assembly;
    assembly.set_exchange_A_and_B_matrices(true);
    structural_discipline.init_system_dirichlet_bc(sys);
    const unsigned int n_eig = 4;
    eq_sys.parameters.set<unsigned int>("eigenpairs")    = n_eig;
    eq_sys.parameters.set<unsigned int>("basis vectors") = n_eig*3;
    
    assembly.attach_discipline_and_system(structural_discipline, structural_sys);
    sys.solve();
    assembly.clear_discipline_and_system();

    // Get the number of converged eigen pairs.
    unsigned int nconv = std::min(sys.get_n_converged(),
                                  n_eig);
    
    for (unsigned int i=0; i<nconv; i++)
    {
        std::ostringstream file_name;
        
        // We write the file in the ExodusII format.
        file_name << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << i
        << ".exo";
        
        // now write the eigenvlaues
        std::pair<Real, Real> val = sys.get_eigenpair(i);
        std::complex<Real> eigval = std::complex<Real>(val.first, val.second);
        eigval = 1./eigval;
        libMesh::out
        << std::setw(35) << std::fixed << std::setprecision(15)
        << eigval.real() << std::endl;
        
        // We write the file in the ExodusII format.
        libMesh::ExodusII_IO(mesh).write_equation_systems(file_name.str(),
                                                          eq_sys);
    }
    
    return 0;
}

