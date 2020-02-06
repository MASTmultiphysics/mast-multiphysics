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

TEST_CASE("quad4_linear_extension_bending_structural", 
          "[quad],[quad4],[linear],[structural],[2D],[element]")
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
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    
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
    
    const Real V0 = reference_elem->volume();
    
    /**
     *  Below, we start building the Jacobian up for a very basic element that 
     *  is already in an isoparametric format.  The following four steps:
     *      Extension Only Stiffness
     *      Extension & Bending Stiffness
     *      Extension, Bending, & Transverse Shear Stiffness
     *      Extension, Bending, & Extension-Bending Coupling Stiffness
     *      Extension, Bending, Transverse Shear & Extension-Bending Coupling Stiffness
     * 
     *  Testing the Jacobian incrementally this way allows errors in the
     *  Jacobian to be located more precisely.
     * 
     *  It would probabyl be even better to test each stiffness contribution 
     *  separately, but currently, we don't have a way to disable extension
     *  stiffness and transverse shear stiffness doesn't exist without bending 
     *  and extension-bending coupling stiffness don't make sense without both
     *  extension and/or bending. 
     */
    
    // Set shear coefficient to zero to disable transverse shear stiffness
    // FIXME: In MAST, when kappa=0.0, this results in zero rows/columns in element stiffness matrix
    kappa = 0.0;
    
    // Set the offset to zero to disable extension-bending coupling
    offset = 0.0;
    
    // Calculate residual and jacobian
    RealVectorX residual = RealVectorX::Zero(n_dofs);
    RealMatrixX jacobian0 = RealMatrixX::Zero(n_dofs, n_dofs);
    elem->internal_residual(true, residual, jacobian0);
            
    double val_margin = (jacobian0.array().abs()).mean() * 1.490116119384766e-08;
    
    
    SECTION("internal_jacobian_finite_difference_check")                   
    {
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        
        std::vector<double> test =  eigen_matrix_to_std_vector(jacobian0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian_fd);
        
        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    SECTION("internal_jacobian_symmetry_check")
    {
        // Element stiffness matrix should be symmetric
        std::vector<double> test =  eigen_matrix_to_std_vector(jacobian0);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian0.transpose());
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("internal_jacobian_determinant_check")
    {
        // Determinant of undeformed element stiffness matrix should be zero
        REQUIRE( jacobian0.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    
    SECTION("internal_jacobian_displacement_invariant")
    {
        // Calculate residual and jacobian at arbitrary displacement
        RealVectorX elem_sol = RealVectorX::Zero(n_dofs);
        elem_sol << -0.04384355,  0.03969142, -0.09470648, -0.05011107, 
                    -0.02989082, -0.01205296,  0.08846868,  0.04522207,  
                     0.06435953, -0.07282706,  0.09307561, -0.06250143,  
                     0.03332844, -0.00040089, -0.00423108, -0.07258241,  
                     0.06636534, -0.08421098, -0.0705489 , -0.06004976,
                     0.03873095, -0.09194373,  0.00055061,  0.046831;
        elem->set_solution(elem_sol);
        
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
                
        std::vector<double> test =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian0);

        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    SECTION("internal_jacobian_shifted_x_invariant")
    {
        // Shifted in x-direction
        transform_element(mesh, X, 5.2, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
                
        std::vector<double> test =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian0);

        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    SECTION("internal_jacobian_shifted_y_invariant")
    {
        // Shifted in y-direction
        transform_element(mesh, X, 0.0, -11.5, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        std::vector<double> test =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> truth = eigen_matrix_to_std_vector(jacobian0);
        
        REQUIRE_THAT( test, Catch::Approx<double>(truth).margin(val_margin) );
    }
    
    SECTION("internal_jacobian_rotated_about_z")
    {
        // Rotated 63.4 about z-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 63.4);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_rotated_about_y")
    {
        // Rotated 35.8 about y-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 35.8, 0.0);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_rotated_about_x")
    {
        // Rotated 15.8 about x-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 15.8, 0.0, 0.0);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_sheared_in_x")
    {
        // Rotated 63.4 about z-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 6.7, 0.0);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_sheared_in_y")
    {
        // Rotated 63.4 about z-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -11.2);
        REQUIRE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_scaled_x")
    {
        // Rotated 63.4 about z-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 3.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        REQUIRE_FALSE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_scaled_y")
    {
        // Rotated 63.4 about z-axis at element's centroid
        transform_element(mesh, X, 0.0, 0.0, 0.0, 1.0, 0.64, 0.0, 0.0, 0.0, 0.0, 0.0);
        REQUIRE_FALSE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_arbitrary_transformation")
    {
        // Arbitrary transformations applied to the element
        transform_element(mesh, X, -5.0, 7.8, -13.1, 2.7, 6.4, 20.0, 47.8, 
                          -70.1, 5.7, -6.3);
        REQUIRE_FALSE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
    
    SECTION("internal_jacobian_arbitrary_with_displacements")
    {
        RealMatrixX X = RealMatrixX::Zero(3,4);
        X << -3.2,  2.4, 1.5, -2.1, 
             -1.0, -0.2, 1.4,  0.8, 
              0.0,  0.0, 0.0,  0.0;
        // Update the element with the new node Coordinates
        for (int i=0; i<X.cols(); i++)
        {
            (*mesh.node_ptr(i)) = libMesh::Point(X(0,i), X(1,i), X(2,i));
        }
        
        // Arbitrary transformations applied to the element
        transform_element(mesh, X, 4.1, -6.3, 7.5, 4.2, 1.5, -18.0, -24.8, 
                          30.1, -3.2, 5.4);
        
        // Calculate residual and jacobian at arbitrary displacement
        RealVectorX elem_sol = RealVectorX::Zero(n_dofs);
        elem_sol << -0.04384355,  0.03969142, -0.09470648, -0.05011107, 
                    -0.02989082, -0.01205296,  0.08846868,  0.04522207,  
                     0.06435953, -0.07282706,  0.09307561, -0.06250143,  
                     0.03332844, -0.00040089, -0.00423108, -0.07258241,  
                     0.06636534, -0.08421098, -0.0705489 , -0.06004976,
                     0.03873095, -0.09194373,  0.00055061,  0.046831;
        elem->set_solution(elem_sol);
        
        REQUIRE_FALSE( reference_elem->volume() == Approx(V0) );
        
        // Calculate residual and jacobian
        RealVectorX residual = RealVectorX::Zero(n_dofs);
        RealMatrixX jacobian = RealMatrixX::Zero(n_dofs, n_dofs);
        elem->internal_residual(true, residual, jacobian);
        
        // Approximate Jacobian with Finite Difference
        RealMatrixX jacobian_fd = RealMatrixX::Zero(n_dofs, n_dofs);
        approximate_internal_jacobian_with_finite_difference(*elem, elem_solution, jacobian_fd);
        
        // This is necessary because MAST manually (hard-coded) adds a small 
        // value to the diagonal to prevent singularities at inactive DOFs
        //double val_margin = (jacobian_fd.array().abs()).maxCoeff() * 1.490116119384766e-08;
        val_margin = (jacobian_fd.array().abs()).mean() * 1.490116119384766e-08;
        //std::cout << "val_margin = " << val_margin << std::endl;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> J =  eigen_matrix_to_std_vector(jacobian);
        std::vector<double> J_fd = eigen_matrix_to_std_vector(jacobian_fd);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( J, Catch::Approx<double>(J_fd).margin(val_margin) );
        
        // Symmetry check
        std::vector<double> Jt = eigen_matrix_to_std_vector(jacobian.transpose());
        REQUIRE_THAT( Jt, Catch::Approx<double>(J) );
        
        // Determinant check
        REQUIRE( jacobian.determinant() == Approx(0.0).margin(1e-06) );
    }
}
