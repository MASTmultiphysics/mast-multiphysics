// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"

// Custom includes
#include "test_helpers.h"

// TODO: Need to test temperature dependent material property
// TODO: Need to test stress dependent material property (plasticity)

TEST_CASE("constant_isotropic_thermoelastic_material_1d",
          "[isotropic],[material],[constant],[1D],[thermoelastic]")
{
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;      
    
    material.add(alpha_f);
    
    REQUIRE( material.contains("alpha_expansion") );
    
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


TEST_CASE("constant_isotropic_structural_material_1d", 
          "[isotropic],[material],[constant],[1D],[structural]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E") );
    REQUIRE( material.contains("nu") );
    
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
        D_true(1,1) = 2.706766917293233e+10;
        
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


TEST_CASE("constant_isotropic_heat_transfer_material_1d", 
          "[isotropic],[material],[1D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k("k_param",     237.0);             // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k) );
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


TEST_CASE("constant_isotropic_transient_heat_transfer_material_1d", 
          "[isotropic],[material],[1D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
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


TEST_CASE("constant_isotropic_thermoelastic_material_2d", 
          "[isotropic],[material],[2D],[thermoelastic][constant]")
{
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;      
    
    material.add(alpha_f);
    
    REQUIRE( material.contains("alpha_expansion") );
    
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
        D_texp_true(1,0) = 5.43e-05;
        
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


TEST_CASE("constant_isotropic_structural_material_2d", 
          "[isotropic],[material],[2D],[structural],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    
    const Real G = E() / (2.0 * (1.0+nu()));
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E") );
    REQUIRE( material.contains("nu") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(E) );
        CHECK( material.depends_on(nu) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D plane stress material constitutive matrix is correct")
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
        K_mat_true(0,0) = 8.079901245651442e+10;
        K_mat_true(1,1) = 8.079901245651442e+10;
        K_mat_true(0,1) = 2.666367411064976e+10;
        K_mat_true(1,0) = 2.666367411064976e+10;
        K_mat_true(2,2) = 2.706766917293233e+10;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        test =  eigen_matrix_to_std_vector(K_mat);
        truth = eigen_matrix_to_std_vector(K_mat_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D plane strain material constitutive matrix is correct")
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
        K_mat_true(0,0) = 1.066784608580274e+11;
        K_mat_true(1,1) = 1.066784608580274e+11;
        K_mat_true(0,1) = 0.525431225121628e+11;
        K_mat_true(1,0) = 0.525431225121628e+11;
        K_mat_true(2,2) = 0.270676691729323e+11;
        
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
        D_trans_shear_true(0,0) = G;//2.255639097744361e+10;
        D_trans_shear_true(1,1) = G;//2.255639097744361e+10;
        
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


TEST_CASE("constant_isotropic_heat_transfer_material_2d", 
          "[isotropic],[material],[2D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k("k_param",     237.0);             // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k) );
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


TEST_CASE("constant_isotropic_transient_heat_transfer_material_2d", 
          "[isotropic],[material],[2D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
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


TEST_CASE("constant_isotropic_thermoelastic_material_3d", 
          "[isotropic],[material],[3D],[thermoelastic][constant]")
{
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;      
    
    material.add(alpha_f);
    
    REQUIRE( material.contains("alpha_expansion") );
    
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
        D_texp_true(1,0) = 5.43e-05;
        D_texp_true(2,0) = 5.43e-05;
        
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


TEST_CASE("constant_isotropic_structural_material_3d", 
          "[isotropic],[material],[3D],[structural],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    
    const Real G = E() / (2.0 * (1.0+nu()));
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E") );
    REQUIRE( material.contains("nu") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(E) );
        CHECK( material.depends_on(nu) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("3D material constitutive matrix is correct")
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
        K_mat_true(0,0) = 1.066784608580274e+11;
        K_mat_true(1,1) = 1.066784608580274e+11;
        K_mat_true(2,2) = 1.066784608580274e+11;
        K_mat_true(3,3) = 0.541353383458647e+11/2.;
        K_mat_true(4,4) = 0.541353383458647e+11/2.;
        K_mat_true(5,5) = 0.541353383458647e+11/2.;
        
        K_mat_true(0,1) = K_mat_true(1,0) = 0.525431225121628e+11;
        K_mat_true(0,2) = K_mat_true(2,0) = 0.525431225121628e+11;
        K_mat_true(1,2) = K_mat_true(2,1) = 0.525431225121628e+11;
        
        
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


TEST_CASE("constant_isotropic_heat_transfer_material_3d", 
          "[isotropic],[material],[3D],[heat_transfer][constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter k("k_param",     237.0);             // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(k_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("k_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(k) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal conductivity matrix")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(3);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_k;
        mat_conduct(point, time, D_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX D_k_true = RealMatrixX::Identity(3,3);
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


TEST_CASE("constant_isotropic_transient_heat_transfer_material_3d", 
          "[isotropic],[material],[3D],[heat_transfer],[constant],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
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


TEST_CASE("constant_isotropic_dynamic_material_3d", 
          "[isotropic],[material],[3D],[dynamic],[constant]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
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
