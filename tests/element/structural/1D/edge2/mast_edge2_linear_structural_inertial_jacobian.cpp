// C++ Stanard Includes
#include <math.h>

// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_section_element_property_card.h"
#include "elasticity/structural_element_1d.h"
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


TEST_CASE("edge2_linear_structural_inertial_consistent",
          "[1D],[dynamic],[edge],[edge2],[element]")
{
    const int n_elems = 1;
    const int n_nodes = 2;
    
    RealMatrixX temp = RealMatrixX::Zero(3,n_nodes);
    temp << -1.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    const RealMatrixX X0 = temp;
    
    /**
     * Create the mesh with the one element we are testing
     */
    libMesh::ReplicatedMesh mesh(p_global_init->comm());
    mesh.set_mesh_dimension(1);
    mesh.set_spatial_dimension(3);
    mesh.reserve_elem(n_elems);
    mesh.reserve_nodes(n_nodes);
    
    // Add nodes to the mesh
    for (uint i=0; i<n_nodes; i++)
    {
        mesh.add_point(libMesh::Point(X0(0,i), X0(1,i), X0(2,i)));
    }
    
    // Add the element to the mesh
    libMesh::Elem *reference_elem = new libMesh::Edge2;
    reference_elem->set_id(0);    
    reference_elem->subdomain_id() = 0;
    reference_elem = mesh.add_elem(reference_elem);
    for (int i=0; i<n_nodes; i++)
    {
        reference_elem->set_node(i) = mesh.node_ptr(i);
    }
    
    // Prepare the mesh for use
    mesh.prepare_for_use();
    
    // Ensure the libMesh element has the expected volume
    const Real elem_volume = reference_elem->volume();
    Real true_volume = 2.0;
    REQUIRE( elem_volume == true_volume );
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness_y("thy_param", 0.8);   // Section thickness in y-direction
    MAST::Parameter thickness_z("thz_param", 0.7);   // Section thickness in z-direction
    MAST::Parameter offset_y("offy_param", 0.5);    // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", 0.4);    // Section offset in z-direction
    MAST::Parameter kappa_zz("kappa_zz", 5.0/6.0);    // Shear coefficient
    MAST::Parameter kappa_yy("kappa_yy", 2.0/6.0);    // Shear coefficient
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    MAST::ConstantFieldFunction thicknessy_f("hy", thickness_y);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction thicknessz_f("hz", thickness_z);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    MAST::ConstantFieldFunction kappa_zz_f("Kappazz", kappa_zz);
    MAST::ConstantFieldFunction kappa_yy_f("Kappayy", kappa_yy);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);
    material.add(nu_f);
    material.add(rho_f);
    material.add(alpha_f);
    material.add(k_f);
    material.add(cp_f);
    
    // Initialize the section
    MAST::Solid1DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thicknessy_f);
    section.add(offsety_f);
    section.add(thicknessz_f);
    section.add(offsetz_f);
    section.add(kappa_zz_f);
    section.add(kappa_yy_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(2) = 1.0;
    section.y_vector() = orientation;
    
    // Set the strain type to linear for the section
    section.set_strain(MAST::LINEAR_STRAIN);
    
    // Set the bending operator to Timoshenko
    section.set_bending_model(MAST::TIMOSHENKO);
    
    // Now initialize the section
    section.init();
    
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
    MAST::StructuralElement1D* elem = (dynamic_cast<MAST::StructuralElement1D*>(elem_base.get()));
    
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
    
    // Calculate residual and jacobian
    RealVectorX residual = RealVectorX::Zero(n_dofs);
    RealMatrixX jac0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xddot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xdot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    elem->inertial_residual(true, residual, jac_xddot0, jac_xdot0, jac0);
            
    double val_margin = (jac_xddot0.array().abs()).mean() * 1.490116119384766e-08;
    
    //libMesh::out << "Jac_xddot0:\n" << jac_xddot0 << std::endl;
    
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
        /**
         * A lumped mass matrix should have all positive eigenvalues since it
         * is a diagonal matrix and masses should not be zero or negative.
         */
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(jac_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        //libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}


TEST_CASE("edge2_linear_structural_inertial_lumped",
          "[1D],[dynamic],[edge],[edge2]")
{
    const int n_elems = 1;
    const int n_nodes = 2;
    
    RealMatrixX temp = RealMatrixX::Zero(3,n_nodes);
    temp << -1.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    const RealMatrixX X0 = temp;
    
    /**
     * Create the mesh with the one element we are testing
     */
    libMesh::ReplicatedMesh mesh(p_global_init->comm());
    mesh.set_mesh_dimension(1);
    mesh.set_spatial_dimension(3);
    mesh.reserve_elem(n_elems);
    mesh.reserve_nodes(n_nodes);
    
    // Add nodes to the mesh
    for (uint i=0; i<n_nodes; i++)
    {
        mesh.add_point(libMesh::Point(X0(0,i), X0(1,i), X0(2,i)));
    }
    
    // Add the element to the mesh
    libMesh::Elem *reference_elem = new libMesh::Edge2;
    reference_elem->set_id(0);    
    reference_elem->subdomain_id() = 0;
    reference_elem = mesh.add_elem(reference_elem);
    for (int i=0; i<n_nodes; i++)
    {
        reference_elem->set_node(i) = mesh.node_ptr(i);
    }
    
    // Prepare the mesh for use
    mesh.prepare_for_use();
    
    // Ensure the libMesh element has the expected volume
    const Real elem_volume = reference_elem->volume();
    Real true_volume = 2.0;
    REQUIRE( elem_volume == true_volume );
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness_y("thy_param", 0.8);   // Section thickness in y-direction
    MAST::Parameter thickness_z("thz_param", 0.7);   // Section thickness in z-direction
    MAST::Parameter offset_y("offy_param", 0.5);    // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", 0.4);    // Section offset in z-direction
    MAST::Parameter kappa_zz("kappa_zz", 5.0/6.0);    // Shear coefficient
    MAST::Parameter kappa_yy("kappa_yy", 2.0/6.0);    // Shear coefficient
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    MAST::ConstantFieldFunction thicknessy_f("hy", thickness_y);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction thicknessz_f("hz", thickness_z);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    MAST::ConstantFieldFunction kappa_zz_f("Kappazz", kappa_zz);
    MAST::ConstantFieldFunction kappa_yy_f("Kappayy", kappa_yy);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);
    material.add(nu_f);
    material.add(rho_f);
    material.add(alpha_f);
    material.add(k_f);
    material.add(cp_f);
    
    // Initialize the section
    MAST::Solid1DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thicknessy_f);
    section.add(offsety_f);
    section.add(thicknessz_f);
    section.add(offsetz_f);
    section.add(kappa_zz_f);
    section.add(kappa_yy_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(2) = 1.0;
    section.y_vector() = orientation;
    
    // Set the strain type to linear for the section
    section.set_strain(MAST::LINEAR_STRAIN);
    
    // Set the bending operator to Euler-Bernoulli
    section.set_bending_model(MAST::TIMOSHENKO);
    
    // Set the mass matrix to be lumped as opposed to the default, consistent
    section.set_diagonal_mass_matrix(true);
    REQUIRE( section.if_diagonal_mass_matrix() );
    
    // Now initialize the section
    section.init();
    
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
    MAST::StructuralElement1D* elem = (dynamic_cast<MAST::StructuralElement1D*>(elem_base.get()));
    
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
    
    // Calculate residual and jacobian
    RealVectorX residual = RealVectorX::Zero(n_dofs);
    RealMatrixX jac0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xddot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    RealMatrixX jac_xdot0 = RealMatrixX::Zero(n_dofs, n_dofs);
    elem->inertial_residual(true, residual, jac_xddot0, jac_xdot0, jac0);
            
    double val_margin = (jac_xddot0.array().abs()).mean() * 1.490116119384766e-08;
    
    //libMesh::out << "Jac_xddot0:\n" << jac_xddot0 << std::endl;
    
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
        // Element interial jacobian should be symmetric
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
        /**
         * A lumped mass matrix should have all positive eigenvalues since it
         * is a diagonal matrix and masses should not be zero or negative.
         */
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(jac_xddot0, false);
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        //libMesh::out << "Eigenvalues are:\n" << eigenvalues << std::endl;
        REQUIRE(eigenvalues.minCoeff()>0.0);
    }
}
