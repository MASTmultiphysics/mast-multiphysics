// C++ Stanard Includes
#include <math.h>

// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/face_quad4.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"
#include "elasticity/structural_element_2d.h"
#include "elasticity/structural_system_initialization.h"
#include "base/physics_discipline_base.h"
#include "base/nonlinear_implicit_assembly.h"
#include "elasticity/structural_nonlinear_assembly.h"
#include "base/nonlinear_system.h"
#include "elasticity/structural_element_base.h"
#include "mesh/geom_elem.h"

// Custom includes
#include "test_helpers.h"

#define pi 3.14159265358979323846

extern libMesh::LibMeshInit* p_global_init;


TEST_CASE("structural_element_2d_base_tests",
          "[2D],[structural],[base]")
{
    const int n_elems = 1;
    const int n_nodes = 4;
    
    // Point Coordinates
    const std::vector<Real> x0 = {-1.0, 1.0, 1.0, -1.0};
    const std::vector<Real> y0 = {-1.0, -1.0, 1.0, 1.0};
    const std::vector<Real> z0 = {0.0, 0.0, 0.0, 0.0};
    
    RealMatrixX temp = RealMatrixX::Zero(3,4);
    temp << -1.0,  1.0, 1.0, -1.0, 
            -1.0, -1.0, 1.0,  1.0, 
             0.0,  0.0, 0.0,  0.0;
    const RealMatrixX X = temp;
    
    std::vector<Real> x = x0;
    std::vector<Real> y = y0;
    std::vector<Real> z = z0;
    
    /**
     *  First create the mesh with the one element we are testing.
     */
    // Setup the mesh properties
    libMesh::ReplicatedMesh mesh(p_global_init->comm());
    mesh.set_mesh_dimension(2);
    mesh.set_spatial_dimension(2);
    mesh.reserve_elem(n_elems);
    mesh.reserve_nodes(n_nodes);
    
    // Add nodes to the mesh
    for (uint i=0; i<n_nodes; i++)
    {
        mesh.add_point(libMesh::Point(x[i], y[i], z[i]));
    }
    
    // Add the element to the mesh
    libMesh::Elem *reference_elem = new libMesh::Quad4;
    reference_elem->set_id(0);    
    reference_elem->subdomain_id() = 0;
    reference_elem = mesh.add_elem(reference_elem);
    for (int i=0; i<n_nodes; i++)
    {
        reference_elem->set_node(i) = mesh.node_ptr(i);
    }
    
    // Prepare the mesh for use
    mesh.prepare_for_use();
    //mesh.print_info();
    
    const Real elem_volume = reference_elem->volume();
    // Calculate true volume using 2D shoelace formula
    Real true_volume = 0.0;
    for (uint i=0; i<n_nodes-1; i++)
    {
        true_volume += x[i]*y[i+1] - x[i+1]*y[i];
    }
    true_volume += x[n_nodes-1]*y[0] - x[0]*y[n_nodes-1];
    true_volume = std::abs(true_volume);
    true_volume *= 0.5;
    
    // Ensure the libMesh element has the expected volume
    REQUIRE( elem_volume == true_volume );
            
    /**
     *  Setup the material and section properties to be used in the element
     */
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    material.add(alpha_f);
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    section.add(kappa_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    
    /**
     *  Now we setup the structural system we will be solving.
     */
    libMesh::EquationSystems equation_systems(mesh);
    
    MAST::NonlinearSystem& system = equation_systems.add_system<MAST::NonlinearSystem>("structural");
    
    libMesh::FEType fetype(libMesh::FIRST, libMesh::LAGRANGE);
    
    MAST::StructuralSystemInitialization structural_system(system, 
                                                           system.name(), 
                                                           fetype);
    
    MAST::PhysicsDisciplineBase discipline(equation_systems);
    
    discipline.set_property_for_subdomain(0, section);
    
    equation_systems.init();
    //equation_systems.print_info();
    
    MAST::NonlinearImplicitAssembly assembly;
    assembly.set_discipline_and_system(discipline, structural_system);
    
    // Create the MAST element from the libMesh reference element
    MAST::GeomElem geom_elem;
    geom_elem.init(*reference_elem, structural_system);
    std::unique_ptr<MAST::StructuralElementBase> elem_base = build_structural_element(structural_system, geom_elem, section);
    
    // Cast the base structural element as a 2D structural element
    MAST::StructuralElement2D* elem = (dynamic_cast<MAST::StructuralElement2D*>(elem_base.get()));
    
    SECTION("number_strain_components")
    {
        REQUIRE(elem->n_direct_strain_components() == 3);
        REQUIRE(elem->n_von_karman_strain_components() == 2);
    }
    
    SECTION("no_incompatible_modes")
    {
        REQUIRE_FALSE( elem->if_incompatible_modes() );
    }
    
    SECTION("return_section_property")
    {
        const MAST::ElementPropertyCardBase& elem_section = elem->elem_property();
        CHECK( elem_section.if_isotropic() );
    }
    
    SECTION("set_get_local_solution")
    {
        const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices (reference_elem, dof_indices);
        uint n_dofs = uint(dof_indices.size());
        
        RealVectorX elem_solution = 5.3*RealVectorX::Ones(n_dofs);
        elem->set_solution(elem_solution);
        
        const RealVectorX& local_solution = elem->local_solution();
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(elem_solution);
        std::vector<double> truth = eigen_matrix_to_std_vector(local_solution);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("set_get_local_solution_sensitivity")
    {
        const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices (reference_elem, dof_indices);
        uint n_dofs = uint(dof_indices.size());
        
        RealVectorX elem_solution_sens = 3.1*RealVectorX::Ones(n_dofs);
        elem->set_solution(elem_solution_sens, true);
        
        const RealVectorX& local_solution_sens = elem->local_solution(true);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(elem_solution_sens);
        std::vector<double> truth = eigen_matrix_to_std_vector(local_solution_sens);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("element shape can be transformed")
    {
        const Real V0 = reference_elem->volume();
        
        // Stretch in x-direction
        transform_element(mesh, X, 0.0, 0.0, 0.0, 3.1, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == 12.4);
        
        // Stretch in y-direction
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 3.1, 0.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == 12.4);
        
        // Rotation about z-axis
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 60.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Rotation about y-axis
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 30.0, 0.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Rotation about x-axis
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 20.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Shifted in x-direction
        transform_element(mesh, X, 10.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Shifted in y-direction
        transform_element(mesh, X, 0.0, 7.5, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Shifted in z-direction
        transform_element(mesh, X, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE(reference_elem->volume() == V0);
        
        // Shear in x
        transform_element(mesh, X, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0, 5.2, 0.0);
        REQUIRE(reference_elem->volume() == Approx(V0));
        
        // Shear in y
        transform_element(mesh, X, 0.0, 0.0, 4.2, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -6.4);
        REQUIRE(reference_elem->volume() == Approx(V0));
    }
}
