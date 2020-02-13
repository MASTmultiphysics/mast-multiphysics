// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/orthotropic_material_property_card.h"

// Custom includes
#include "test_helpers.h"

// TODO: Need to test temperature dependent material property
// TODO: Need to test stress dependent material property (plasticity)
// TODO: Need to test material property where field functions are not constant
// TODO: Need to implement tests for 3D elements


TEST_CASE("constant_orthotropic_thermoelastic_material_1d",
          "[orthotropic],[material],[constant],[1D],[thermoelastic]")
{
    MAST::Parameter alpha11("alpha11_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::ConstantFieldFunction alpha11_f("alpha11_expansion", alpha11);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;      
    
    material.add(alpha11_f);
    
    REQUIRE( material.contains("alpha11_expansion") );
    
    SECTION("1D material thermal expansion matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_texp =
            material.thermal_expansion_matrix(1);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_texp;
        mat_texp(point, time, D_texp);
        
        // Hard-coded in the true value of material thermal expansion matrix
        RealMatrixX D_texp_true = RealMatrixX::Zero(2,1);
        D_texp_true(0,0) = 5.43e-05;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_texp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_texp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_structural_material_1d", 
          "[orthotropic],[material],[constant],[1D],[structural]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E11("E11_param", 72.0e9);   // Modulus of Elasticity
    MAST::Parameter G12("G12_param", 68.3e9);   // Shear Modulus
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E11_f("E11", E11);
    MAST::ConstantFieldFunction G12_f("G12", G12);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E11_f);    
    material.add(G12_f);    

    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E11") );
    REQUIRE( material.contains("G12") );
    
    SECTION("1D material constitutive matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(1);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D;
        mat_stiffness(point, time, D);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX D_true = RealMatrixX::Zero(2,2);
        D_true(0,0) = 72.0e9;
        D_true(1,1) = 68.3e9;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_heat_transfer_material_1d", 
          "[orthotropic],[material],[1D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k11("k11_param",     237.0);     // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k11_f("k11_th", k11);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k11_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k11_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k11) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("1D thermal conductivity matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(1);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_k;
        mat_conduct(point, time, D_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_k_true = RealMatrixX::Identity(1,1);
        D_k_true *= 237.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_k);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_k_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_transient_heat_transfer_material_1d", 
          "[orthotropic],[material],[1D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    material.add(cp_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    REQUIRE( material.contains("cp") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
        CHECK( material.depends_on(cp) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("1D thermal capacitance matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(1);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_cp;
        mat_capacit(point, time, D_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_cp_true = RealMatrixX::Identity(1,1);
        D_cp_true *= (908.0 * 1234.5);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_thermoelastic_material_2d", 
          "[orthotropic],[material],[2D],[thermoelastic][constant]")
{
    MAST::Parameter alpha11("alpha11_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::Parameter alpha22("alpha22_param", 8.79e-06);   // Coefficient of thermal expansion
    
    MAST::ConstantFieldFunction alpha11_f("alpha11_expansion", alpha11);
    MAST::ConstantFieldFunction alpha22_f("alpha22_expansion", alpha22);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;      
    
    material.add(alpha11_f);
    material.add(alpha22_f);
    
    REQUIRE( material.contains("alpha11_expansion") );
    REQUIRE( material.contains("alpha22_expansion") );
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D material thermal expansion matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_texp =
            material.thermal_expansion_matrix(2);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_texp;
        mat_texp(point, time, D_texp);
        
        // Hard-coded in the true value of material thermal expansion matrix
        RealMatrixX D_texp_true = RealMatrixX::Zero(3,1);
        D_texp_true(0,0) = 5.43e-05;
        D_texp_true(1,0) = 8.79e-06;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_texp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_texp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_structural_material_2d", 
          "[orthotropic],[material],[2D],[structural],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E11("E11", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter E22("E22", 55.0e9);           
    MAST::Parameter E33("E33", 33.4e9);
    MAST::Parameter nu12("nu12", 0.33);             // Poisson's ratio
    MAST::Parameter nu23("nu23", 0.25);
    MAST::Parameter nu31("nu31", 0.45);
    MAST::Parameter G12("G12", 68.3e9);
        
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E11_f("E11", E11);
    MAST::ConstantFieldFunction E22_f("E22", E22);
    MAST::ConstantFieldFunction E33_f("E33", E33);
    MAST::ConstantFieldFunction nu12_f("nu12", nu12);
    MAST::ConstantFieldFunction nu23_f("nu23", nu23);
    MAST::ConstantFieldFunction nu31_f("nu31", nu31);
    MAST::ConstantFieldFunction G12_f("G12", G12);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E11_f);         
    material.add(E22_f);
    material.add(E33_f);
    material.add(nu12_f);
    material.add(nu23_f);
    material.add(nu31_f);
    material.add(G12_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E11") );
    REQUIRE( material.contains("E22") );
    REQUIRE( material.contains("E33") );
    REQUIRE( material.contains("nu12") );
    REQUIRE( material.contains("nu23") );
    REQUIRE( material.contains("nu31") );
    REQUIRE( material.contains("G12") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(E11) );
        CHECK( material.depends_on(E22) );
        CHECK( material.depends_on(E33) );
        CHECK( material.depends_on(nu12) );
        CHECK( material.depends_on(nu23) );
        CHECK( material.depends_on(nu31) );
        CHECK( material.depends_on(G12) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D plane stress material constitutive matrix")
    {
        const bool is_plane_stress = true;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX K_mat;
        mat_stiffness(point, time, K_mat);
        
        // Check for symmetry
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(K_mat.transpose());
        std::vector<double> truth = eigen_matrix_to_std_vector(K_mat);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
        
        // Check for positive definiteness
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(K_mat);
        if (eigensolver.info() != Success)
        {
            std::cout << "eigensolver failed to converge!" << std::endl;
            REQUIRE(false); // Force test failure.
        }
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        for (uint i=0; i<K_mat.cols(); i++)
        {
            REQUIRE( eigenvalues(i) > 0.0 );
        }
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX K_mat_true = RealMatrixX::Zero(3,3);
        K_mat_true(0,0) = 7.853296066534869e+10;
        K_mat_true(1,1) = 5.999045606380803e+10;
        K_mat_true(2,2) = 6.830000000000000e+10;
        K_mat_true(0,1) = K_mat_true(1,0) =1.979685050105665e+10;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        test =  eigen_matrix_to_std_vector(K_mat);
        truth = eigen_matrix_to_std_vector(K_mat_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D plane strain material constitutive matrix")
    {
        const bool is_plane_stress = false;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX K_mat;
        mat_stiffness(point, time, K_mat);
        
        // Check for symmetry
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(K_mat.transpose());
        std::vector<double> truth = eigen_matrix_to_std_vector(K_mat);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
        
        // Check for positive definiteness
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(K_mat);
        if (eigensolver.info() != Success)
        {
            std::cout << "eigensolver failed to converge!" << std::endl;
            REQUIRE(false); // Force test failure.
        }
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        for (uint i=0; i<K_mat.cols(); i++)
        {
            REQUIRE( eigenvalues(i) > 0.0 );
        }
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX K_mat_true = RealMatrixX::Zero(3,3);
        K_mat_true(0,0) = 1.881848591463047e11;
        K_mat_true(1,1) = 0.841961884847417e11;
        K_mat_true(2,2) = 0.683000000000000e11;
        K_mat_true(0,1) = K_mat_true(1,0) = 0.713158228712175e11;        
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        test =  eigen_matrix_to_std_vector(K_mat);
        truth = eigen_matrix_to_std_vector(K_mat_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("transverse shear material stiffness matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_trans_shear =
            material.transverse_shear_stiffness_matrix();
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_trans_shear;
        mat_trans_shear(point, time, D_trans_shear);
        
        // Hard-coded in the true value of the material transverse shear stiffness
        RealMatrixX D_trans_shear_true = RealMatrixX::Zero(2,2);
        D_trans_shear_true(0,0) = 6.830000000000000e+10;
        D_trans_shear_true(1,1) = 6.830000000000000e+10;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(D_trans_shear.data(), D_trans_shear.data()+D_trans_shear.rows()*D_trans_shear.cols());
        std::vector<double> truth(D_trans_shear_true.data(), D_trans_shear_true.data()+D_trans_shear_true.rows()*D_trans_shear_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_heat_transfer_material_2d", 
          "[orthotropic],[material],[2D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k11("k11_param",     237.0);   // Thermal Conductivity
    MAST::Parameter k22("k22_param",     542.0);   // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k11_f("k11_th", k11);
    MAST::ConstantFieldFunction k22_f("k22_th", k22);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k11_f);
    material.add(k22_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k11_th") );
    REQUIRE( material.contains("k22_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k11) );
        CHECK( material.depends_on(k22) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal conductivity matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(2);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_k;
        mat_conduct(point, time, D_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_k_true = RealMatrixX::Identity(2,2);
        D_k_true(0,0) = 237.0;
        D_k_true(1,1) = 542.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_k);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_k_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_transient_heat_transfer_material_2d", 
          "[orthotropic],[material],[2D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    material.add(cp_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    REQUIRE( material.contains("cp") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
        CHECK( material.depends_on(cp) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal capacitance matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(2);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_cp;
        mat_capacit(point, time, D_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_cp_true = RealMatrixX::Identity(1,1);
        D_cp_true *= (908.0 * 1234.5);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_thermoelastic_material_3d", 
          "[orthotropic],[material],[3D],[thermoelastic][constant]")
{
    MAST::Parameter alpha11("alpha11_param", 5.43e-05);
    MAST::Parameter alpha22("alpha22_param", 8.79e-06);
    MAST::Parameter alpha33("alpha33_param", 2.14e-05);
    
    MAST::ConstantFieldFunction alpha11_f("alpha11_expansion", alpha11);
    MAST::ConstantFieldFunction alpha22_f("alpha22_expansion", alpha22);
    MAST::ConstantFieldFunction alpha33_f("alpha33_expansion", alpha33);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;      
    
    material.add(alpha11_f);
    material.add(alpha22_f);
    material.add(alpha33_f);
    
    REQUIRE( material.contains("alpha11_expansion") );
    REQUIRE( material.contains("alpha22_expansion") );
    REQUIRE( material.contains("alpha33_expansion") );
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("3D material thermal expansion matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_texp =
            material.thermal_expansion_matrix(3);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_texp;
        mat_texp(point, time, D_texp);
        
        // Hard-coded in the true value of material thermal expansion matrix
        RealMatrixX D_texp_true = RealMatrixX::Zero(6,1);
        D_texp_true(0,0) = 5.43e-05;
        D_texp_true(1,0) = 8.79e-06;
        D_texp_true(2,0) = 2.14e-05;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_texp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_texp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_structural_material_3d", 
          "[orthotropic],[material],[3D],[structural],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E11("E11", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter E22("E22", 55.0e9);           
    MAST::Parameter E33("E33", 33.4e9);
    MAST::Parameter nu12("nu12", 0.33);             // Poisson's ratio
    MAST::Parameter nu23("nu23", 0.25);
    MAST::Parameter nu31("nu31", 0.45);
    MAST::Parameter G12("G12", 68.3e9);             // Shear Modulus
    MAST::Parameter G23("G23", 74.5e9);
    MAST::Parameter G13("G13", 55.1e9);
        
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E11_f("E11", E11);
    MAST::ConstantFieldFunction E22_f("E22", E22);
    MAST::ConstantFieldFunction E33_f("E33", E33);
    MAST::ConstantFieldFunction nu12_f("nu12", nu12);
    MAST::ConstantFieldFunction nu23_f("nu23", nu23);
    MAST::ConstantFieldFunction nu31_f("nu31", nu31);
    MAST::ConstantFieldFunction G12_f("G12", G12);
    MAST::ConstantFieldFunction G23_f("G23", G23);
    MAST::ConstantFieldFunction G13_f("G13", G13);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E11_f);         
    material.add(E22_f);
    material.add(E33_f);
    material.add(nu12_f);
    material.add(nu23_f);
    material.add(nu31_f);
    material.add(G12_f);
    material.add(G23_f);
    material.add(G13_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E11") );
    REQUIRE( material.contains("E22") );
    REQUIRE( material.contains("E33") );
    REQUIRE( material.contains("nu12") );
    REQUIRE( material.contains("nu23") );
    REQUIRE( material.contains("nu31") );
    REQUIRE( material.contains("G12") );
    REQUIRE( material.contains("G23") );
    REQUIRE( material.contains("G13") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(E11) );
        CHECK( material.depends_on(E22) );
        CHECK( material.depends_on(E33) );
        CHECK( material.depends_on(nu12) );
        CHECK( material.depends_on(nu23) );
        CHECK( material.depends_on(nu31) );
        CHECK( material.depends_on(G12) );
        CHECK( material.depends_on(G23) );
        CHECK( material.depends_on(G13) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("3D material constitutive matrix")
    {
        const bool is_plane_stress = true;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(3);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX K_mat;
        mat_stiffness(point, time, K_mat);
        
        // Check for symmetry
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(K_mat.transpose());
        std::vector<double> truth = eigen_matrix_to_std_vector(K_mat);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        REQUIRE_THAT( test, Catch::Approx<double>(truth) );
        
        // Check for positive definiteness
        SelfAdjointEigenSolver<RealMatrixX> eigensolver(K_mat);
        if (eigensolver.info() != Success)
        {
            std::cout << "eigensolver failed to converge!" << std::endl;
            REQUIRE(false); // Force test failure.
        }
        RealVectorX eigenvalues = eigensolver.eigenvalues();
        for (uint i=0; i<K_mat.cols(); i++)
        {
            REQUIRE( eigenvalues(i) > 0.0 );
        }
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX K_mat_true = RealMatrixX::Zero(6,6);
        K_mat_true(0,0) = 1.881848591463047e+11;
        K_mat_true(1,1) = 0.841961884847417e+11;
        K_mat_true(2,2) = 0.831923864531179e+11;
        K_mat_true(3,3) = 0.683000000000000e+11;
        K_mat_true(4,4) = 0.745000000000000e+11;
        K_mat_true(5,5) = 0.551000000000000e+11;
    
        K_mat_true(0,1) = K_mat_true(1,0) = 0.713158228712175e+11;
        K_mat_true(0,2) = K_mat_true(2,0) = 0.955102251790129e+11;
        K_mat_true(1,2) = K_mat_true(2,1) = 0.448746325438223e+11;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        test =  eigen_matrix_to_std_vector(K_mat);
        truth = eigen_matrix_to_std_vector(K_mat_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_heat_transfer_material_3d", 
          "[orthotropic],[material],[3D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k11("k11_param",     237.0);   // Thermal Conductivity
    MAST::Parameter k22("k22_param",     542.0);   // Thermal Conductivity
    MAST::Parameter k33("k33_param",     444.4);   // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k11_f("k11_th", k11);
    MAST::ConstantFieldFunction k22_f("k22_th", k22);
    MAST::ConstantFieldFunction k33_f("k33_th", k33);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k11_f);
    material.add(k22_f);
    material.add(k33_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k11_th") );
    REQUIRE( material.contains("k22_th") );
    REQUIRE( material.contains("k33_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k11) );
        CHECK( material.depends_on(k22) );
        CHECK( material.depends_on(k33) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("3D thermal conductivity matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(3);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_k;
        mat_conduct(point, time, D_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_k_true = RealMatrixX::Identity(3,3);
        D_k_true(0,0) = 237.0;
        D_k_true(1,1) = 542.0;
        D_k_true(2,2) = 444.4;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_k);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_k_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_transient_heat_transfer_material_3d", 
          "[orthotropic],[material],[3D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    material.add(cp_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    REQUIRE( material.contains("cp") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
        CHECK( material.depends_on(cp) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal capacitance matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(3);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_cp;
        mat_capacit(point, time, D_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_cp_true = RealMatrixX::Identity(1,1);
        D_cp_true *= (908.0 * 1234.5);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_orthotropic_dynamic_material_3d", 
          "[orthotropic],[material],[3D],[dynamic],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    
    // Initialize the material
    MAST::OrthotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("3D inertia matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.inertia_matrix(3);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_inertia;
        mat_capacit(point, time, D_inertia);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_inertia_true = RealMatrixX::Identity(3,3);
        D_inertia_true *= 1234.5;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(D_inertia);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_inertia_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}
