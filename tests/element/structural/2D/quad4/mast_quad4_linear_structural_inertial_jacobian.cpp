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

TEST_CASE("quad4_linear_structural_inertial_consistent", 
          "[quad],[quad4],[dynamic],[2D],[element]")
{
    const int n_elems = 1;
    const int n_nodes = 4;
    
    // Point Coordinates
    RealMatrixX temp = RealMatrixX::Zero(3,4);
    temp << -1.0,  1.0, 1.0, -1.0, 
            -1.0, -1.0, 1.0,  1.0, 
             0.0,  0.0, 0.0,  0.0;
    const RealMatrixX X = temp;
    
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
        mesh.add_point(libMesh::Point(X(0,i), X(1,i), X(2,i)));
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
    Real true_volume = get_shoelace_area(X);
    
    // Ensure the libMesh element has the expected volume
    REQUIRE( elem_volume == true_volume );
            
    /**
     *  Setup the material and section properties to be used in the element
     */
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    material.add(rho_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    section.add(kappa_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    
    // Set the strain type to linear for the section
    section.set_strain(MAST::LINEAR_STRAIN);
    
    // Set the bending operator
    section.set_bending_model(MAST::MINDLIN);
    
    // Set the mass matrix to be lumped as opposed to the default, consistent
    section.set_diagonal_mass_matrix(false);
    REQUIRE_FALSE( section.if_diagonal_mass_matrix() );
    
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
    
    // Get element DOFs
    const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    dof_map.dof_indices (reference_elem, dof_indices);
    uint n_dofs = uint(dof_indices.size());
    
    // Set element's initial solution and solution sensitivity to zero
    RealVectorX elem_solution = RealVectorX::Zero(n_dofs);
    elem->set_solution(elem_solution);
    elem->set_solution(elem_solution, true);
    
    // Set element's initial acceleration and acceleration sensitivity to zero
    RealVectorX elem_accel = RealVectorX::Zero(n_dofs);
    elem->set_acceleration(elem_accel);
    elem->set_acceleration(elem_accel, true);
    
    const Real V0 = reference_elem->volume();
    
    // Calculate inertial residual and jacobians
    RealVectorX residual = RealVectorX::Zero(n_dofs);
    RealMatrixX jac0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xddot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xdot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    elem->inertial_residual(true, residual, jac_xddot0, jac_xdot0, jac0);
            
    double val_margin = (jac_xddot0.array().abs()).mean() * 1.490116119384766e-08;
    
    SECTION("inertial_jacobian_finite_difference_check")                   
    {
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_inertial_jacobian_with_finite_difference(*elem, elem_accel, jacobian_fd);
        
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        
        std::vector<double> test =  eigen_matrix_to_std_vector(jac_xddot0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian_fd);
        
        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    
    SECTION("inertial_jacobian_symmetry_check")
    {
        // Element inertial jacobian should be symmetric
        std::vector<double> test =  eigen_matrix_to_std_vector(jac_xddot0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jac_xddot0.transpose());
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("inertial_jacobian_determinant_check")
    {
        // Determinant of inertial jacobian should be positive
        REQUIRE( jac_xddot0.determinant() > 0.0 );
    }
    
    
    SECTION("inertial_jacobian_eigenvalue_check")
    {
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(jac_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}


TEST_CASE("quad4_linear_structural_inertial_lumped", 
          "[quad],[quad4],[dynamic],[2D],[element]")
{
    const int n_elems = 1;
    const int n_nodes = 4;
    
    // Point Coordinates
    RealMatrixX temp = RealMatrixX::Zero(3,4);
    temp << -1.0,  1.0, 1.0, -1.0, 
            -1.0, -1.0, 1.0,  1.0, 
             0.0,  0.0, 0.0,  0.0;
    const RealMatrixX X = temp;
    
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
        mesh.add_point(libMesh::Point(X(0,i), X(1,i), X(2,i)));
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
    Real true_volume = get_shoelace_area(X);
    
    // Ensure the libMesh element has the expected volume
    REQUIRE( elem_volume == true_volume );
            
    /**
     *  Setup the material and section properties to be used in the element
     */
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    material.add(rho_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    section.add(kappa_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    // Set the strain type to linear for the section
    section.set_strain(MAST::LINEAR_STRAIN);
    
    // Set the bending operator
    section.set_bending_model(MAST::MINDLIN);
    
    // Set the mass matrix to be lumped as opposed to the default, consistent
    section.set_diagonal_mass_matrix(true);
    REQUIRE( section.if_diagonal_mass_matrix() );
    
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
    
    // Get element DOFs
    const libMesh::DofMap& dof_map = assembly.system().get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    dof_map.dof_indices (reference_elem, dof_indices);
    uint n_dofs = uint(dof_indices.size());
    
    // Set element's initial solution and solution sensitivity to zero
    RealVectorX elem_solution = RealVectorX::Zero(n_dofs);
    elem->set_solution(elem_solution);
    elem->set_solution(elem_solution, true);
    
    // Set element's initial acceleration and acceleration sensitivity to zero
    RealVectorX elem_accel = RealVectorX::Zero(n_dofs);
    elem->set_acceleration(elem_accel);
    elem->set_acceleration(elem_accel, true);
    
    const Real V0 = reference_elem->volume();
    
    // Calculate inertial residual and jacobians
    RealVectorX residual = RealVectorX::Zero(n_dofs);
    RealMatrixX jac0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xddot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xdot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    elem->inertial_residual(true, residual, jac_xddot0, jac_xdot0, jac0);
            
    double val_margin = (jac_xddot0.array().abs()).mean() * 1.490116119384766e-08;
    
    SECTION("inertial_jacobian_finite_difference_check")                   
    {
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_inertial_jacobian_with_finite_difference(*elem, elem_accel, jacobian_fd);
        
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        
        std::vector<double> test =  eigen_matrix_to_std_vector(jac_xddot0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian_fd);
        
        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    
    SECTION("inertial_jacobian_symmetry_check")
    {
        // Element inertial jacobian should be symmetric
        std::vector<double> test =  eigen_matrix_to_std_vector(jac_xddot0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jac_xddot0.transpose());
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("inertial_jacobian_determinant_check")
    {
        // Determinant of inertial jacobian should be positive
        REQUIRE( jac_xddot0.determinant() > 0.0 );
    }
    
    
    SECTION("inertial_jacobian_eigenvalue_check")
    {
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(jac_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}
