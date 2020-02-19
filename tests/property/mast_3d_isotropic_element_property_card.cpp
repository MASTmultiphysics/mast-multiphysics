// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/isotropic_element_property_card_3D.h"

// Custom includes
#include "test_helpers.h"

// TODO: Need to test with other types of materials.
// TODO: Need to test function that gets material.

extern libMesh::LibMeshInit* p_global_init;

TEST_CASE("element_property_card_constant_heat_transfer_isotropic_3d",
          "[heat_transfer],[3D],[isotropic],[constant],[property]")
{
    const uint dim = 3;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k_f("k_th", k);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(k_f);
    
    // Initialize the section
    MAST::IsotropicElementPropertyCard3D section;
        
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim);
    REQUIRE( section.depends_on(k) );
    
    SECTION("3D section thermal conductance matrix")
    {
        /** 
         * As of Dec. 17, 2019, inertia_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * inertia_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> conduct_mat = section.thermal_conductance_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_conduc;
        conduct_mat->operator()(point, time, D_sec_conduc);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_conduc_true = RealMatrixX::Zero(3,3);
        D_sec_conduc_true(0,0) = 237.0;
        D_sec_conduc_true(1,1) = 237.0;
        D_sec_conduc_true(2,2) = 237.0;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_conduc);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_conduc_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("element_property_card_constant_transient_heat_transfer_isotropic_3d",
          "[heat_transfer],[3D],[isotropic],[constant],[property],[transient]")
{
    const uint dim = 3;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
        
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(cp_f);
    
    // Initialize the section
    MAST::IsotropicElementPropertyCard3D section;
    
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim);
    REQUIRE( section.depends_on(cp) );
    REQUIRE( section.depends_on(rho) );
    
    SECTION("3D section thermal capacitance matrix")
    {
        /** 
         * As of Dec. 17, 2019, inertia_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * inertia_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> capaci_mat = section.thermal_capacitance_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_capac;
        capaci_mat->operator()(point, time, D_sec_capac);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_capac_true = RealMatrixX::Zero(1,1);
        D_sec_capac_true(0,0) = 908.0*1420.5;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_capac);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_capac_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("element_property_card_constant_thermoelastic_isotropic_3d",
          "[thermoelastic],[3D],[isotropic],[constant],[property]")
{
    const uint dim = 3;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
        
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);
    
    // Initialize the section
    MAST::IsotropicElementPropertyCard3D section;
        
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(alpha) );
    
    SECTION("3D plane stress thermal expansion A matrix")
    {
        /*!
         *  thermal expansion A matrix is defined as D * alpha_vec
         */
                        
        /** 
         * As of Dec. 17, 2019, thermal_expansion_A_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * thermal_expansion_A_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> texp_A_mat = section.thermal_expansion_A_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_texpA;
        texp_A_mat->operator()(point, time, D_sec_texpA);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_texpA_true = RealMatrixX::Zero(6,1);
        D_sec_texpA_true(0,0) = 1.149882352941177e+07;
        D_sec_texpA_true(1,0) = 1.149882352941177e+07;
        D_sec_texpA_true(2,0) = 1.149882352941177e+07;
        
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpA);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpA_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("3D plane stress thermal expansion B matrix")
    {
        /*!
         *  thermal expansion B matrix is defined as D * alpha_vec
         */
                
        /** 
         * As of Dec. 17, 2019, thermal_expansion_B_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * thermal_expansion_B_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> texp_B_mat = section.thermal_expansion_B_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_texpB;
        texp_B_mat->operator()(point, time, D_sec_texpB);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_texpB_true = RealMatrixX::Zero(6,1);
        D_sec_texpB_true(0,0) = 1.149882352941177e+07;
        D_sec_texpB_true(1,0) = 1.149882352941177e+07;
        D_sec_texpB_true(2,0) = 1.149882352941177e+07;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpB);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpB_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}


TEST_CASE("element_property_card_constant_dynamic_isotropic_3d",
          "[dynamic],[3D],[isotropic],[constant],[property]")
{
    const uint dim = 3;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
        
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(E_f);                                             
    material.add(nu_f);
    
    // Initialize the section
    MAST::IsotropicElementPropertyCard3D section;
    
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(E) );
    REQUIRE( section.depends_on(nu) );
    
    REQUIRE_FALSE( section.if_diagonal_mass_matrix() );
    
    section.set_diagonal_mass_matrix(true);
    REQUIRE( section.if_diagonal_mass_matrix() );
    
    SECTION("3D section inertia matrix")
    {
        /** 
         * As of Dec. 17, 2019, inertia_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * inertia_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> inertia_mat = section.inertia_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_inertia;
        inertia_mat->operator()(point, time, D_sec_inertia);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_inertia_true = RealMatrixX::Identity(3,3);        
        
        D_sec_inertia_true *= 1420.5;

        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_inertia);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_inertia_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}



TEST_CASE("element_property_card_constant_structural_isotropic_3d",
          "[structural],[3D],[isotropic],[constant],[property]")
{
    const uint dim = 3;
    
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
    
    // Initialize the section
    MAST::IsotropicElementPropertyCard3D section;
        
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(E) );
    REQUIRE( section.depends_on(nu) );
    
    SECTION("solid_2d_section is isotropic")
    {
        CHECK( section.if_isotropic() );
    }
    
    SECTION("set_get_strain_type")
    {
        // Check that the default is linear strain
        REQUIRE( section.strain_type() == MAST::LINEAR_STRAIN );
        
        section.set_strain(MAST::LINEAR_STRAIN);
        REQUIRE( section.strain_type() == MAST::LINEAR_STRAIN );
        REQUIRE_FALSE( section.strain_type() == MAST::NONLINEAR_STRAIN );
        
        section.set_strain(MAST::NONLINEAR_STRAIN);
        REQUIRE( section.strain_type() == MAST::NONLINEAR_STRAIN );
        REQUIRE_FALSE( section.strain_type() == MAST::LINEAR_STRAIN );
    }
    
//     SECTION("quadrature_order")
//     {
//         section.set_bending_model(MAST::MINDLIN);
//         REQUIRE( section.extra_quadrature_order(elem) == 0 );
//         
//         section.set_bending_model(MAST::DKT);
//         REQUIRE( section.extra_quadrature_order(elem) == 2 );
//     }
    
    SECTION("3D stiffness matrix")
    {
        /** 
         * As of Dec. 17, 2019, stiffness_A_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * stiffness_A_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> stiffness_mat = section.stiffness_A_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_stiff;
        stiffness_mat->operator()(point, time, D_stiff);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_stiff_true = RealMatrixX::Zero(6,6);
        D_stiff_true(0,0) = 1.066784608580274e+11;
        D_stiff_true(1,1) = 1.066784608580274e+11;
        D_stiff_true(2,2) = 1.066784608580274e+11;
        D_stiff_true(3,3) = 0.541353383458647e+11/2.0;
        D_stiff_true(4,4) = 0.541353383458647e+11/2.0;
        D_stiff_true(5,5) = 0.541353383458647e+11/2.0;
        
        D_stiff_true(0,1) = D_stiff_true(1,0) = 0.525431225121628e+11;
        D_stiff_true(0,2) = D_stiff_true(2,0) = 0.525431225121628e+11;
        D_stiff_true(1,2) = D_stiff_true(2,1) = 0.525431225121628e+11;
        
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_stiff);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_stiff_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
   
}
