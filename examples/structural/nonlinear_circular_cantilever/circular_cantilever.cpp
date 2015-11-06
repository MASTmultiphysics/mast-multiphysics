
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
#include "property_cards/isotropic_element_property_card_3D.h"
#include "driver/driver_base.h"


// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"


class StructuralSystemInitializationUVW:
public MAST::SystemInitialization  {
    
public:
    StructuralSystemInitializationUVW(libMesh::System& sys,
                                      const std::string& prefix,
                                      const libMesh::FEType& fe_type):
    MAST::SystemInitialization(sys, prefix) {
        
        _vars.resize(3);
        
        std::string nm = prefix + "_ux";
        _vars[0] = sys.add_variable(nm, fe_type);
        
        nm = prefix + "_uy";
        _vars[1] = sys.add_variable(nm, fe_type);
        
        nm = prefix + "_uz";
        _vars[2] = sys.add_variable(nm, fe_type);
    }
    
    
    virtual ~StructuralSystemInitializationUVW() { }
    
protected:
    
};


int main_1(int argc, const char * argv[]) {
    
    // initialize the libMesh object
    libMesh::LibMeshInit              init(argc, argv);
    libMesh::ParallelMesh             mesh(init.comm());
    libMesh::EquationSystems          eq_sys(mesh);
    
    // add the system to be used for analysis
    libMesh::System& sys =
    eq_sys.add_system<libMesh::NonlinearImplicitSystem>("structural");
    
    // initialize the mesh
    const unsigned int
    nx           = 10, // elems along length
    ny           = 4,  // elems in thickness direction
    nz           = 4,  // elems in width direction
    n_load_steps = 10;
    libMesh::MeshTools::Generation::build_cube(mesh,
                                               nx, ny, nz,
                                               0.,  1.,
                                               -.5, .5,
                                               -.5, .5);
    mesh.prepare_for_use();
    
    // modify the mesh to a circular arc of 45 degrees
    const Real
    r = 100.,
    pi    = acos(-1.),
    theta = pi/180.*45.;
    
    libMesh::MeshBase::node_iterator
    node_it  = mesh.nodes_begin(),
    node_end = mesh.nodes_end();
    
    Real eta = 0.;
    std::set<libMesh::Node*> end_nodes;
    
    for ( ; node_it != node_end; node_it++) {
        libMesh::Point& pt = **node_it;
        eta = pt(0);
        pt(0) = (r+pt(1))*cos(pi/2.-theta*eta);
        pt(1) = (r+pt(1))*sin(pi/2.-theta*eta);
        
        // also prepare the set of nodes on the right side for
        // application of point force
        if (eta == 1.)
            end_nodes.insert(*node_it);
    }
    
    
    // variable type
    libMesh::FEType fe_type(libMesh::FIRST,
                            libMesh::LAGRANGE);
    
    MAST::StructuralDiscipline        structural_discipline(eq_sys);
    StructuralSystemInitializationUVW structural_sys(sys,
                                                     sys.name(),
                                                     fe_type);
    
    
    // add the Dirichlet boundary condition on left boundary
    MAST::DirichletBoundaryCondition   dirichlet;
    dirichlet.init(4, structural_sys.vars());
    
    // add the boundary condition to the physics
    structural_discipline.add_dirichlet_bc(4, dirichlet);
    
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
    E  (  "E",  10.0e6),   // Young's modulus
    nu (  "nu",   0.0);   // Poisson's ratio
    
    
    // define material property based on the parameter
    MAST::ConstantFieldFunction
    E_f ("E",  E),
    nu_f("nu", nu);
    
    
    // add them to a material property card
    MAST::IsotropicMaterialPropertyCard       materials;
    materials.add( E_f);
    materials.add(nu_f);
    
    
    // add to a section property card
    MAST::IsotropicElementPropertyCard3D section_property;
    
    // tell the section property to use the material property card
    section_property.set_material(materials);
    
    // tell the conduction physics about the section properties
    structural_discipline.set_property_for_subdomain(0, section_property);
    
    // create the nonlinear assembly object
    MAST::StructuralNonlinearAssembly   assembly;

    // write the solution for visualization
    libMesh::ExodusII_IO out(mesh);

    // now solve the system
    for (unsigned int i=0; i<n_load_steps; i++) {

        p = 600.*(i+1)/n_load_steps;
        
        assembly.attach_discipline_and_system(structural_discipline, structural_sys);
        
        libMesh::NonlinearImplicitSystem&      nonlin_sys   =
        dynamic_cast<libMesh::NonlinearImplicitSystem&>(assembly.system());
        
        nonlin_sys.solve();
        
        assembly.clear_discipline_and_system();
        
        sys.solution->print();
        
        out.write_timestep("out.exo", eq_sys, i, (i+1.)/n_load_steps);
    }
    
    
    return 0;
}

