// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_2d_section_element_property_card.h"

// Custom includes
#include "test_helpers.h"

// TODO: Need to test with other types of materials.
// TODO: Need to test function that gets material.

extern libMesh::LibMeshInit* p_global_init;

TEST_CASE("element_property_card_constant_heat_transfer_2d",
          "[heat_transfer],[2D],[isotropic],[constant],[property]")
{
    const uint dim = 2;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction k_f("k_th", k);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(k_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(thickness) );
    REQUIRE( section.depends_on(offset) );
    REQUIRE( section.depends_on(k) );
    
    SECTION("2D section thermal conductance matrix")
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
        RealMatrixX D_sec_conduc_true = RealMatrixX::Zero(2,2);
        D_sec_conduc_true(0,0) = 237.0*0.06;
        D_sec_conduc_true(1,1) = 237.0*0.06;
        
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


TEST_CASE("element_property_card_constant_transient_heat_transfer_2d",
          "[heat_transfer],[2D],[isotropic],[constant],[property],[transient]")
{
    const uint dim = 2;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(cp_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(thickness) );
    REQUIRE( section.depends_on(offset) );
    REQUIRE( section.depends_on(cp) );
    REQUIRE( section.depends_on(rho) );
    
    SECTION("2D section thermal capacitance matrix")
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
        D_sec_capac_true(0,0) = 908.0*1420.5*0.06;
        
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


TEST_CASE("element_property_card_constant_thermoelastic_2d",
          "[thermoelastic],[2D],[isotropic],[constant],[property]")
{
    const uint dim = 2;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter alpha("alpha_param", 5.43e-05);   // Coefficient of thermal expansion
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);
    
    // Initialize the section
    MAST::Solid2DSectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(thickness_f);
    section.add(offset_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(thickness) );
    REQUIRE( section.depends_on(offset) );
    REQUIRE( section.depends_on(alpha) );
    
    SECTION("2D plane stress thermal expansion A matrix")
    {
        /*!
         *  thermal expansion A matrix is defined as D * alpha_vec * thickness
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
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
        RealMatrixX D_sec_texpA_true = RealMatrixX::Zero(3,1);
        D_sec_texpA_true(0,0) = 3.501134328358209e+05;
        D_sec_texpA_true(1,0) = 3.501134328358209e+05;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpA);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpA_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane strain thermal expansion A matrix")
    {
        /*!
         *  thermal expansion A matrix is defined as D * alpha_vec * thickness
         */
        
        section.set_plane_stress(false);
        REQUIRE_FALSE( section.plane_stress() );
                
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
        RealMatrixX D_sec_texpA_true = RealMatrixX::Zero(3,1);
        D_sec_texpA_true(0,0) = 5.187439186200796e+05;
        D_sec_texpA_true(1,0) = 5.187439186200796e+05;
                
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpA);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpA_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane stress thermal expansion B matrix")
    {
        /*!
         *  thermal expansion B matrix is defined as D * alpha_vec * thickness * offset
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
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
        RealMatrixX D_sec_texpB_true = RealMatrixX::Zero(3,1);
        D_sec_texpB_true(0,0) = 1.050340298507463e+04;
        D_sec_texpB_true(1,0) = 1.050340298507463e+04;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpB);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpB_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("2D plane strain thermal expansion B matrix")
    {
        /*!
         *  thermal expansion B matrix is defined as D * alpha_vec * thickness * offset
         */
        
        section.set_plane_stress(false);
        REQUIRE_FALSE( section.plane_stress() );
                
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
        RealMatrixX D_sec_texpB_true = RealMatrixX::Zero(3,1);
        D_sec_texpB_true(0,0) = 1.556231755860239e+04;
        D_sec_texpB_true(1,0) = 1.556231755860239e+04;
        
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


TEST_CASE("element_property_card_constant_dynamic_2d",
          "[dynamic],[2D],[isotropic],[constant],[property]")
{
    const uint dim = 2;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter thickness("th_param", 0.06);      // Section thickness
    MAST::Parameter offset("off_param", 0.03);        // Section offset
    
    // Create field functions to dsitribute these constant parameters throughout the model
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
    MAST::ConstantFieldFunction thickness_f("h", thickness);
    MAST::ConstantFieldFunction offset_f("off", offset);
    
    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;                   
    
    // Add the material property constant field functions to the material card
    material.add(rho_f);
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
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(thickness) );
    REQUIRE( section.depends_on(offset) );
    REQUIRE( section.depends_on(E) );
    REQUIRE( section.depends_on(nu) );
    REQUIRE( section.depends_on(kappa) );
    
    REQUIRE_FALSE( section.if_diagonal_mass_matrix() );
    
    section.set_diagonal_mass_matrix(true);
    REQUIRE( section.if_diagonal_mass_matrix() );
    
    SECTION("2D section inertia matrix")
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
        RealMatrixX D_sec_iner;
        inertia_mat->operator()(point, time, D_sec_iner);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_iner_true = RealMatrixX::Zero(6,6);
        D_sec_iner_true(0,0) = D_sec_iner_true(1,1) = D_sec_iner_true(2,2) = 0.06;
        D_sec_iner_true(0,4) = D_sec_iner_true(4,0) = 0.03*0.06;
        D_sec_iner_true(1,3) = D_sec_iner_true(3,1) = -0.03*0.06;
        D_sec_iner_true(3,3) = D_sec_iner_true(4,4) = pow(0.06,3.0)/12.0 + 0.06*0.03*0.03;
        D_sec_iner_true(5,5) = pow(0.06,3.0)/12.0 * 1.0e-6; // Ignore this term
        
        D_sec_iner_true *= 1420.5;

        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_iner);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_iner_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}



TEST_CASE("element_property_card_constant_structural_2d",
          "[structural],[2D],[isotropic],[constant],[property]")
{
    const uint dim = 2;
    
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
    
    REQUIRE( section.dim() == dim); // Ensure section is 2 dimensional
    REQUIRE( section.depends_on(thickness) );
    REQUIRE( section.depends_on(offset) );
    REQUIRE( section.depends_on(E) );
    REQUIRE( section.depends_on(nu) );
    REQUIRE( section.depends_on(kappa) );
    
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
    
    SECTION("set_get_bending_model")
    {
        // NOTE: MAST::BERNOULLI AND MAST::TIMOSHENKO are not valid options for 2D sections, even though their input is accepted.
        section.set_bending_model(MAST::DKT);
        section.set_bending_model(MAST::DEFAULT_BENDING);
        section.set_bending_model(MAST::NO_BENDING);
        section.set_bending_model(MAST::MINDLIN);
        
        // TODO: Implement element creation for testing of getting bending_model and checking default
        //REQUIRE( section.bending_model()
    }
    
    SECTION("set_get_plane_stress_strain")
    {
        // Check that default is plane stress
        REQUIRE( section.plane_stress() == true );
        
        section.set_plane_stress(false);
        REQUIRE( section.plane_stress() == false );
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() == true);
    }
    
    
//     SECTION("quadrature_order")
//     {
//         section.set_bending_model(MAST::MINDLIN);
//         REQUIRE( section.extra_quadrature_order(elem) == 0 );
//         
//         section.set_bending_model(MAST::DKT);
//         REQUIRE( section.extra_quadrature_order(elem) == 2 );
//     }
    
    SECTION("2D plane stress extension stiffness matrix")
    {
        /*!
         *  extension sitffness matrix is defined as D * thickness
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
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
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> extension_stiffness_mat = section.stiffness_A_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_ext;
        extension_stiffness_mat->operator()(point, time, D_sec_ext);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_ext_true = RealMatrixX::Zero(3,3);
        D_sec_ext_true(0,0) = 4.847940747390865e+09;
        D_sec_ext_true(1,1) = 4.847940747390865e+09;
        D_sec_ext_true(0,1) = 1.599820446638986e+09;
        D_sec_ext_true(1,0) = 1.599820446638986e+09;
        D_sec_ext_true(2,2) = 1.624060150375940e+09;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_ext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_ext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane strain extension stiffness matrix")
    {
        /*!
         *  extension sitffness matrix is defined as D * thickness
         */
        
        section.set_plane_stress(false);
        REQUIRE_FALSE( section.plane_stress() );
                
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
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> extension_stiffness_mat = section.stiffness_A_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_ext;
        extension_stiffness_mat->operator()(point, time, D_sec_ext);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_ext_true = RealMatrixX::Zero(3,3);
        D_sec_ext_true(0,0) = 6.400707651481644e+09;
        D_sec_ext_true(1,1) = 6.400707651481644e+09;
        D_sec_ext_true(0,1) = 3.152587350729766e+09;
        D_sec_ext_true(1,0) = 3.152587350729766e+09;
        D_sec_ext_true(2,2) = 1.624060150375940e+09;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_ext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_ext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane stress bending section stiffness matrix")
    {
        /*!
         *  bending section stiffness matrix is defined as D * (thickness^3 / 12 +  thickness*offset^2)
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
        /** 
         * As of Dec. 17, 2019, stiffness_D_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * stiffness_D_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bending_stiffness_mat = section.stiffness_D_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_bnd;
        bending_stiffness_mat->operator()(point, time, D_sec_bnd);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bnd_true = RealMatrixX::Zero(3,3);
        D_sec_bnd_true(0,0) = 5.817528896869037e+06;
        D_sec_bnd_true(1,1) = 5.817528896869037e+06;
        D_sec_bnd_true(0,1) = 1.919784535966782e+06;
        D_sec_bnd_true(1,0) = 1.919784535966782e+06;
        D_sec_bnd_true(2,2) = 1.948872180451127e+06;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bnd);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bnd_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane strain bending section stiffness matrix")
    {
        /*!
         *  bending section stiffness matrix is defined as D * (thickness^3 / 12 +  thickness*offset^2)
         */
        
        section.set_plane_stress(false);
        REQUIRE_FALSE( section.plane_stress() );
                
        /** 
         * As of Dec. 17, 2019, stiffness_D_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * stiffness_D_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bending_stiffness_mat = section.stiffness_D_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_bnd;
        bending_stiffness_mat->operator()(point, time, D_sec_bnd);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bnd_true = RealMatrixX::Zero(3,3);
        D_sec_bnd_true(0,0) = 7.680849181777972e+06;
        D_sec_bnd_true(1,1) = 7.680849181777972e+06;
        D_sec_bnd_true(0,1) = 3.783104820875719e+06;
        D_sec_bnd_true(1,0) = 3.783104820875719e+06;
        D_sec_bnd_true(2,2) = 1.948872180451128e+06;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bnd);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bnd_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane stress extension-bending section stiffness matrix")
    {
        /*!
         *  extension-bending section stiffness matrix is defined as 
         *  D * thickness * offset
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
        /** 
         * As of Dec. 17, 2019, stiffness_B_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * stiffness_B_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bndext_stiffness_mat = section.stiffness_B_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_bndext;
        bndext_stiffness_mat->operator()(point, time, D_sec_bndext);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bndext_true = RealMatrixX::Zero(3,3);
        D_sec_bndext_true(0,0) = 1.454382224217259e+08;
        D_sec_bndext_true(1,1) = 1.454382224217259e+08;
        D_sec_bndext_true(0,1) = 0.479946133991696e+08;
        D_sec_bndext_true(1,0) = 0.479946133991696e+08;
        D_sec_bndext_true(2,2) = 0.487218045112782e+08;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bndext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bndext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D plane strain extension-bending section stiffness matrix")
    {
        /*!
         *  extension-bending section stiffness matrix is defined as 
         *  D * thickness * offset
         */
        
        section.set_plane_stress(false);
        REQUIRE_FALSE( section.plane_stress() );
                
        /** 
         * As of Dec. 17, 2019, stiffness_B_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * stiffness_B_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bndext_stiffness_mat = section.stiffness_B_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_bndext;
        bndext_stiffness_mat->operator()(point, time, D_sec_bndext);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bndext_true = RealMatrixX::Zero(3,3);
        D_sec_bndext_true(0,0) = 1.920212295444493e+08;
        D_sec_bndext_true(1,1) = 1.920212295444493e+08;
        D_sec_bndext_true(0,1) = 0.945776205218930e+08;
        D_sec_bndext_true(1,0) = 0.945776205218930e+08;
        D_sec_bndext_true(2,2) = 0.487218045112782e+08;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bndext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bndext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    
    SECTION("2D transverse shear section stiffness matrix")
    {
        /*!
         * transverse shear section stiffness matrix is defined as
         * D * thickness * kappa
         */
        
        section.set_plane_stress(true);
        REQUIRE( section.plane_stress() );
                
        /** 
         * As of Dec. 17, 2019, transverse_shear_stiffness_matrix requires the input of an
         * MAST::ElementBase object, but does not actually use it. To get 
         * around requiring the creation of such an object (and therefore 
         * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
         * and MAST::GeomElem objects as well).
         * 
         * To remedy this, an additional method was added to MAST which allows
         * transverse_shear_stiffness_matrix to be obtained without any input arguments.
         */
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> trans_shear_stiffness_mat = section.transverse_shear_stiffness_matrix();
        
        const libMesh::Point point(2.3, 3.1, 5.2);
        const Real time = 2.34;
        RealMatrixX D_sec_shr;
        trans_shear_stiffness_mat->operator()(point, time, D_sec_shr);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_shr_true = RealMatrixX::Zero(2,2);
        D_sec_shr_true(0,0) = 1.353383458646616e+09;
        D_sec_shr_true(1,1) = 1.353383458646616e+09;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_shr);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_shr_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
}
