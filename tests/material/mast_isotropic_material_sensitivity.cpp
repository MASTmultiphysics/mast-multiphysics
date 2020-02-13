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
// TODO: Need to test material property where field functions are not constant
// TODO: Need to implement tests for 3D elements

// TODO: Check this. May need to remove the effect of Kappa from the calculations.

TEST_CASE("constant_isotropic_thermoelastic_material_1d_sensitivity",
          "[isotropic],[material],[constant],[1D],[thermoelastic],[sensitivity]")
{
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;      
    
    material.add(alpha_f);
    
    REQUIRE( material.contains("alpha_expansion") );
    
    SECTION("1D material thermal expansion matrix derivative w.r.t. coefficient of thermal expansion")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_texp =
            material.thermal_expansion_matrix(1);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_texp;
        mat_texp.derivative(alpha, point, time, dD_texp);
        
        // Hard-coded in the true value of material thermal expansion matrix
        RealMatrixX dD_texp_true = RealMatrixX::Zero(2,1);
        dD_texp_true(0,0) = 1.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(dD_texp);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_texp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_structural_material_1d_sensitivity", 
          "[isotropic],[material],[constant],[1D],[structural],[sensitivity]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    material.add(kappa_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E") );
    REQUIRE( material.contains("nu") );
    REQUIRE( material.contains("kappa") );
    
    SECTION("1D material constitutive matrix derivative w.r.t. elastic modulus")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(1);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD;
        mat_stiffness.derivative(E, point, time, dD);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dD_true = RealMatrixX::Zero(2,2);
        dD_true(0,0) = 1.0;
        dD_true(1,1) = 0.375939849624060;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(dD);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("1D material constitutive matrix derivative w.r.t. poisson's ratio")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(1);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD;
        mat_stiffness.derivative(nu, point, time, dD);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dD_true = RealMatrixX::Zero(2,2);
        dD_true(0,0) = 0.0;
        dD_true(1,1) = -0.203516309570920e+11;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(dD);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_heat_transfer_material_1d_sensitivity", 
          "[isotropic],[material],[1D],[heat_transfer],[constant],[sensitivity]")
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
    
    SECTION("1D thermal conductivity matrix derivative w.r.t. thermal conductivity")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(1);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_k;
        mat_conduct.derivative(k, point, time, dD_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX dD_k_true = RealMatrixX::Identity(1,1);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(dD_k);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_k_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_transient_heat_transfer_material_1d_sensitivity", 
          "[isotropic],[material],[1D],[heat_transfer],[constant],[transient],[sensitivity]")
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
    
    SECTION("1D thermal capacitance matrix derivative w.r.t. specifc heat capacity")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(1);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_cp;
        mat_capacit.derivative(cp, point, time, dD_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX dD_cp_true = RealMatrixX::Identity(1,1);
        dD_cp_true *= 1234.5;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(dD_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_thermoelastic_material_2d_sensitivity", 
          "[isotropic],[material],[2D],[thermoelastic],[constant],[sensitivity]")
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
    
    SECTION("2D material thermal expansion matrix derivative w.r.t. coefficient of thermal expansion")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_texp =
            material.thermal_expansion_matrix(2);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_texp;
        mat_texp.derivative(alpha, point, time, dD_texp);
        
        // Hard-coded in the true value of material thermal expansion matrix
        RealMatrixX dD_texp_true = RealMatrixX::Zero(3,1);
        dD_texp_true(0,0) = 1.0;
        dD_texp_true(1,0) = 1.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(dD_texp);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_texp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_structural_material_2d_sensitivity", 
          "[isotropic],[material],[2D],[structural],[constant],[sensitivity]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);                                             
    material.add(nu_f);
    material.add(kappa_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("E") );
    REQUIRE( material.contains("nu") );
    REQUIRE( material.contains("kappa") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(E) );
        CHECK( material.depends_on(nu) );
        CHECK( material.depends_on(kappa) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D plane stress material constitutive matrix derivative w.r.t. elastic modulus is correct")
    {
        const bool is_plane_stress = true;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dK_mat;
        mat_stiffness.derivative(E, point, time, dK_mat);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dK_mat_true = RealMatrixX::Zero(3,3);
        dK_mat_true(0,0) = 1.122208506340478;
        dK_mat_true(1,1) = 1.122208506340478;
        dK_mat_true(0,1) = 0.370328807092358;
        dK_mat_true(1,0) = 0.370328807092358;
        dK_mat_true(2,2) = 0.375939849624060;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dK_mat.data(), dK_mat.data()+dK_mat.rows()*dK_mat.cols());
        std::vector<double> truth(dK_mat_true.data(), dK_mat_true.data()+dK_mat_true.rows()*dK_mat_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane stress material constitutive matrix derivative w.r.t. poisson's ratio is correct")
    {
        const bool is_plane_stress = true;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dK_mat;
        mat_stiffness.derivative(nu, point, time, dK_mat);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dK_mat_true = RealMatrixX::Zero(3,3);
        dK_mat_true(0,0) = 0.598444037945231e+11;
        dK_mat_true(1,1) = 0.598444037945231e+11;
        dK_mat_true(0,1) = 1.005476657087070e+11;
        dK_mat_true(1,0) = 1.005476657087070e+11;
        dK_mat_true(2,2) = -0.203516309570920e+11;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dK_mat.data(), dK_mat.data()+dK_mat.rows()*dK_mat.cols());
        std::vector<double> truth(dK_mat_true.data(), dK_mat_true.data()+dK_mat_true.rows()*dK_mat_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D plane strain material constitutive matrix derivative w.r.t. elastic modulus is correct")
    {
        const bool is_plane_stress = false;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dK_mat;
        mat_stiffness.derivative(E, point, time, dK_mat);
        
//         // 4th order central finite difference approximation
//         RealMatrixX K_mat_h, K_mat_h2, K_mat_n, K_mat_n2;
//         Real delta = 1e9;
//         const Real E0 = E();
//         
//         E = E0 + delta;
//         mat_stiffness.derivative(E, point, time, K_mat_h);
//         
//         E = E0 + 2.0*delta;
//         mat_stiffness.derivative(E, point, time, K_mat_h2);
//         
//         E = E0 - delta;
//         mat_stiffness.derivative(E, point, time, K_mat_n);
//         
//         E = E0 - 2.0*delta;
//         mat_stiffness.derivative(E, point, time, K_mat_n2);
//         
//         E = E0;
//         
//         RealMatrixX dK_mat_true = (K_mat_n2 - 8.0*K_mat_n + 8.0*K_mat_h - K_mat_h2)/(12.0*delta);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dK_mat_true = RealMatrixX::Zero(3,3);
        dK_mat_true(0,0) = 1.481645289694825;
        dK_mat_true(1,1) = 1.481645289694825;
        dK_mat_true(0,1) = 0.729765590446705;
        dK_mat_true(1,0) = 0.729765590446705;
        dK_mat_true(2,2) = 0.375939849624060;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dK_mat.data(), dK_mat.data()+dK_mat.rows()*dK_mat.cols());
        std::vector<double> truth(dK_mat_true.data(), dK_mat_true.data()+dK_mat_true.rows()*dK_mat_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D plane strain material constitutive matrix derivative w.r.t. poisson's ratio is correct")
    {
        const bool is_plane_stress = false;
        const MAST::FieldFunction<RealMatrixX>& mat_stiffness = 
            material.stiffness_matrix(2, is_plane_stress);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dK_mat;
        mat_stiffness.derivative(nu, point, time, dK_mat);
        
        // Hard-code in the true value of the material stiffness
        RealMatrixX dK_mat_true = RealMatrixX::Zero(3,3);
        dK_mat_true(0,0) = 3.880894055520205e+11;
        dK_mat_true(1,1) = 3.880894055520205e+11;
        dK_mat_true(0,1) = 4.287926674662045e+11;
        dK_mat_true(1,0) = 4.287926674662045e+11;
        dK_mat_true(2,2) = -0.203516309570920e+11;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dK_mat.data(), dK_mat.data()+dK_mat.rows()*dK_mat.cols());
        std::vector<double> truth(dK_mat_true.data(), dK_mat_true.data()+dK_mat_true.rows()*dK_mat_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("transverse shear material stiffness matrix derivative w.r.t. elastic modulus")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_trans_shear =
            material.transverse_shear_stiffness_matrix();
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_trans_shear;
        mat_trans_shear.derivative(E, point, time, dD_trans_shear);
        
        // Hard-coded in the true value of the material transverse shear stiffness
        RealMatrixX dD_trans_shear_true = RealMatrixX::Zero(2,2);
        dD_trans_shear_true(0,0) = 0.313283208020050;
        dD_trans_shear_true(1,1) = 0.313283208020050;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dD_trans_shear.data(), dD_trans_shear.data()+dD_trans_shear.rows()*dD_trans_shear.cols());
        std::vector<double> truth(dD_trans_shear_true.data(), dD_trans_shear_true.data()+dD_trans_shear_true.rows()*dD_trans_shear_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("transverse shear material stiffness matrix derivative w.r.t. poisson's ratio")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_trans_shear =
            material.transverse_shear_stiffness_matrix();
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_trans_shear;
        mat_trans_shear.derivative(nu, point, time, dD_trans_shear);
        
        // Hard-coded in the true value of the material transverse shear stiffness
        RealMatrixX dD_trans_shear_true = RealMatrixX::Zero(2,2);
        dD_trans_shear_true(0,0) = -1.695969246424331e+10;
        dD_trans_shear_true(1,1) = -1.695969246424331e+10;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test(dD_trans_shear.data(), dD_trans_shear.data()+dD_trans_shear.rows()*dD_trans_shear.cols());
        std::vector<double> truth(dD_trans_shear_true.data(), dD_trans_shear_true.data()+dD_trans_shear_true.rows()*dD_trans_shear_true.cols());
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("transverse shear material stiffness matrix derivative w.r.t. shear coefficient")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_trans_shear =
            material.transverse_shear_stiffness_matrix();
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_trans_shear;
        mat_trans_shear.derivative(kappa, point, time, dD_trans_shear);
        
        // Hard-coded in the true value of the material transverse shear stiffness
        RealMatrixX dD_trans_shear_true = RealMatrixX::Zero(2,2);
        dD_trans_shear_true(0,0) = 2.706766917293233e+10;
        dD_trans_shear_true(1,1) = 2.706766917293233e+10;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(dD_trans_shear);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_trans_shear_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_transient_heat_transfer_material_2d_sensitivity",
          "[isotropic],[material],[2D],[heat_transfer],[constant],[sensitivity],[transient]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);             // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    material.add(cp_f);
    material.add(k_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    REQUIRE( material.contains("cp") );
    REQUIRE( material.contains("k_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
        CHECK( material.depends_on(cp) );
        CHECK( material.depends_on(k) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal capacitance matrix derivative w.r.t. specifc heat capacity")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(2);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_cp;
        mat_capacit.derivative(cp, point, time, dD_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX dD_cp_true = RealMatrixX::Identity(1,1);
        dD_cp_true *= 1234.5;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(dD_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D thermal capacitance matrix derivative w.r.t. density")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_capacit =
            material.capacitance_matrix(2);
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_cp;
        mat_capacit.derivative(rho, point, time, dD_cp);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX dD_cp_true = RealMatrixX::Identity(1,1);
        dD_cp_true *= 908.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(dD_cp);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_cp_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("constant_isotropic_heat_transfer_material_2d_sensitivity", 
          "[isotropic],[material],[2D],[heat_transfer][constant],[sensitivity]")
{
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1234.5);            // Density
    MAST::Parameter cp("cp_param",   908.0);             // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);             // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);                                             
    material.add(cp_f);
    material.add(k_f);
    
    // Ensure that the material properties were added before doing other checks
    REQUIRE( material.contains("rho") );
    REQUIRE( material.contains("cp") );
    REQUIRE( material.contains("k_th") );
    
    SECTION("material depends on the parameters that it should")
    {
        CHECK( material.depends_on(rho) );
        CHECK( material.depends_on(cp) );
        CHECK( material.depends_on(k) );
    }
    
    SECTION("material does not depend on other parameters")
    {
        MAST::Parameter dummy("dummy", 1.0);
        CHECK_FALSE( material.depends_on(dummy) );
    }
    
    SECTION("2D thermal conductivity matrix derivative w.r.t. thermal conductivity")
    {
        const MAST::FieldFunction<RealMatrixX>& mat_conduct = 
            material.conductance_matrix(2);
            
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX dD_k;
        mat_conduct.derivative(k, point, time, dD_k);
        
        // Hard-coded values for thermal conductivity matrix
        RealMatrixX dD_k_true = RealMatrixX::Identity(2,2);
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test = eigen_matrix_to_std_vector(dD_k);
        std::vector<double> truth = eigen_matrix_to_std_vector(dD_k_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}
