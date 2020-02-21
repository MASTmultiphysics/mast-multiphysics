// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_tube_section_element_property_card.h"

// Custom includes
#include "test_helpers.h"

extern libMesh::LibMeshInit* p_global_init;

#define PI 3.1415926535897932


TEST_CASE("tube_element_property_card_constant_base_1d",
          "[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real cE = 72.0e9;
    const Real cnu = 0.33;
    const Real cG = cE / (2.0 * (1.0+cnu) );
    const Real calpha = 5.43e-05;
    const Real crho = 1420.5;
    const Real ccp = 908.0;
    const Real ck = 237.0;
    const Real ckappa = 5.3284279639775167e-01;
    
    const Real coff_y = 0.35;
    const Real coff_z = 0.26;

    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param",            cE);     // Modulus of Elasticity
    MAST::Parameter nu("nu_param",          cnu);    // Poisson's ratio
    MAST::Parameter alpha("alpha_param",    calpha); // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param",        crho);   // Density
    MAST::Parameter cp("cp_param",          ccp);    // Specific Heat Capacity
    MAST::Parameter k("k_param",            ck);     // Thermal Conductivity

    // Define Section Properties as MAST Parameters
    MAST::Parameter r_o("DIM1", 1.125);   // Outer radius
    MAST::Parameter r_i("DIM2", 0.750);   // Inner radius
    MAST::Parameter offset_y("offy_param", coff_y);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", coff_z);     // Section offset in z-direction

    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", r_o);
    MAST::ConstantFieldFunction DIM2_f("DIM2", r_i);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);

    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;    

    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);

    // Initialize the section
    MAST::Solid1DTubeSectionElementPropertyCard section;


    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);

    // Add the material card to the section card
    section.set_material(material);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Now initialize the section
    section.init(*p_global_init);
    
    // True values
    const Real ro = r_o();
    const Real ri = r_i();
    const Real area_true = PI*(pow(ro,2.0) - pow(ri,2.0));
    const Real first_area_moment_z_true = area_true * coff_y;
    const Real first_area_moment_y_true = area_true * coff_z;
    const Real Izzc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Iyyc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Izyc_true = 0.0;
    const Real Ipc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/2.0;
    const Real second_area_moment_zz_true = Izzc_true + area_true * coff_y * coff_y;
    const Real second_area_moment_yy_true = Iyyc_true + area_true * coff_z * coff_z;
    const Real second_area_moment_zy_true = Izyc_true + area_true * coff_y * coff_z;
    const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
    
    const Real torsion_constant_true = 2.0184430643017226e+00;
    const Real warping_constant_true = 2.7224768298299495e-11;
    const Real kappa_z_true = 5.3284267813657360e-01;
    const Real kappa_y_true = 5.3284267085782810e-01;
    const Real xs_true = coff_z;
    const Real ys_true = coff_y;
    const Real xc_true = coff_z;
    const Real yc_true = coff_y;
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    REQUIRE( section.dim() == dim); // Ensure section is 1 dimensional
    REQUIRE( section.depends_on(r_o) );
    REQUIRE( section.depends_on(offset_y) );
    REQUIRE( section.depends_on(r_i) );
    REQUIRE( section.depends_on(offset_z) );
    REQUIRE( section.depends_on(k) );
    REQUIRE( section.depends_on(cp) );
    REQUIRE( section.depends_on(rho) );
    REQUIRE( section.if_isotropic() );
    
    REQUIRE_FALSE( section.if_diagonal_mass_matrix() );
    
    section.set_diagonal_mass_matrix(true);
    REQUIRE( section.if_diagonal_mass_matrix() );
}


TEST_CASE("tube_element_property_card_constant_heat_transfer_1d",
          "[heat_transfer],[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real cE = 72.0e9;
    const Real cnu = 0.33;
    const Real cG = cE / (2.0 * (1.0+cnu) );
    const Real calpha = 5.43e-05;
    const Real crho = 1420.5;
    const Real ccp = 908.0;
    const Real ck = 237.0;
    const Real ckappa = 5.3284279639775167e-01;
    
    const Real coff_y = 0.35;
    const Real coff_z = 0.26;

    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param",            cE);     // Modulus of Elasticity
    MAST::Parameter nu("nu_param",          cnu);    // Poisson's ratio
    MAST::Parameter alpha("alpha_param",    calpha); // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param",        crho);   // Density
    MAST::Parameter cp("cp_param",          ccp);    // Specific Heat Capacity
    MAST::Parameter k("k_param",            ck);     // Thermal Conductivity

    // Define Section Properties as MAST Parameters
    MAST::Parameter r_o("DIM1", 1.125);   // Outer radius
    MAST::Parameter r_i("DIM2", 0.750);   // Inner radius
    MAST::Parameter offset_y("offy_param", coff_y);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", coff_z);     // Section offset in z-direction

    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", r_o);
    MAST::ConstantFieldFunction DIM2_f("DIM2", r_i);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);

    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;    

    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);

    // Initialize the section
    MAST::Solid1DTubeSectionElementPropertyCard section;


    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);

    // Add the material card to the section card
    section.set_material(material);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Now initialize the section
    section.init(*p_global_init);
    
    // True values
    const Real ro = r_o();
    const Real ri = r_i();
    const Real area_true = PI*(pow(ro,2.0) - pow(ri,2.0));
    const Real first_area_moment_z_true = area_true * coff_y;
    const Real first_area_moment_y_true = area_true * coff_z;
    const Real Izzc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Iyyc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Izyc_true = 0.0;
    const Real Ipc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/2.0;
    const Real second_area_moment_zz_true = Izzc_true + area_true * coff_y * coff_y;
    const Real second_area_moment_yy_true = Iyyc_true + area_true * coff_z * coff_z;
    const Real second_area_moment_zy_true = Izyc_true + area_true * coff_y * coff_z;
    const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
    
    const Real torsion_constant_true = 2.0184430643017226e+00;
    const Real warping_constant_true = 2.7224768298299495e-11;
    const Real kappa_z_true = 5.3284267813657360e-01;
    const Real kappa_y_true = 5.3284267085782810e-01;
    const Real xs_true = coff_z;
    const Real ys_true = coff_y;
    const Real xc_true = coff_z;
    const Real yc_true = coff_y;
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    Real area;
    const MAST::FieldFunction<Real>& Area = section.A();
    Area(point, time, area);
    REQUIRE( area == Approx(area_true) );
    

    SECTION("1D section thermal conductance matrix")
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
        
        RealMatrixX D_sec_conduc;
        conduct_mat->operator()(point, time, D_sec_conduc);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_conduc_true = RealMatrixX::Zero(1,1);
        D_sec_conduc_true(0,0) = ck * area_true;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_conduc);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_conduc_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("1D section thermal capacitance matrix")
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
        
        RealMatrixX D_sec_capac;
        capaci_mat->operator()(point, time, D_sec_capac);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_capac_true = RealMatrixX::Zero(1,1);
        D_sec_capac_true(0,0) = ccp * crho * area_true;
        
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


TEST_CASE("tube_element_property_card_constant_thermoelastic_1d",
          "[thermoelastic],[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real cE = 72.0e9;
    const Real cnu = 0.33;
    const Real cG = cE / (2.0 * (1.0+cnu) );
    const Real calpha = 5.43e-05;
    const Real crho = 1420.5;
    const Real ccp = 908.0;
    const Real ck = 237.0;
    const Real ckappa = 5.3284279639775167e-01;
    
    const Real coff_y = 0.35;
    const Real coff_z = 0.26;

    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param",            cE);     // Modulus of Elasticity
    MAST::Parameter nu("nu_param",          cnu);    // Poisson's ratio
    MAST::Parameter alpha("alpha_param",    calpha); // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param",        crho);   // Density
    MAST::Parameter cp("cp_param",          ccp);    // Specific Heat Capacity
    MAST::Parameter k("k_param",            ck);     // Thermal Conductivity

    // Define Section Properties as MAST Parameters
    MAST::Parameter r_o("DIM1", 1.125);   // Outer radius
    MAST::Parameter r_i("DIM2", 0.750);   // Inner radius
    MAST::Parameter offset_y("offy_param", coff_y);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", coff_z);     // Section offset in z-direction

    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", r_o);
    MAST::ConstantFieldFunction DIM2_f("DIM2", r_i);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);

    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;    

    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);

    // Initialize the section
    MAST::Solid1DTubeSectionElementPropertyCard section;


    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);

    // Add the material card to the section card
    section.set_material(material);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Now initialize the section
    section.init(*p_global_init);
    
    // True values
    const Real ro = r_o();
    const Real ri = r_i();
    const Real area_true = PI*(pow(ro,2.0) - pow(ri,2.0));
    const Real first_area_moment_z_true = area_true * coff_y;
    const Real first_area_moment_y_true = area_true * coff_z;
    const Real Izzc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Iyyc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Izyc_true = 0.0;
    const Real Ipc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/2.0;
    const Real second_area_moment_zz_true = Izzc_true + area_true * coff_y * coff_y;
    const Real second_area_moment_yy_true = Iyyc_true + area_true * coff_z * coff_z;
    const Real second_area_moment_zy_true = Izyc_true + area_true * coff_y * coff_z;
    const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
    
    const Real torsion_constant_true = 2.0184430643017226e+00;
    const Real warping_constant_true = 2.7224768298299495e-11;
    const Real kappa_z_true = 5.3284267813657360e-01;
    const Real kappa_y_true = 5.3284267085782810e-01;
    const Real xs_true = coff_z;
    const Real ys_true = coff_y;
    const Real xc_true = coff_z;
    const Real yc_true = coff_y;
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    SECTION("1D thermal expansion A matrix")
    {
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
        
        Real area;
        const MAST::FieldFunction<Real>& Area = section.A();
        Area(point, time, area);
        CHECK( area == Approx(area_true) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> texp_A_mat = section.thermal_expansion_A_matrix();
        
        RealMatrixX D_sec_texpA;
        texp_A_mat->operator()(point, time, D_sec_texpA);
                
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_texpA_true = RealMatrixX::Zero(2,1);
        D_sec_texpA_true(0,0) = cE * calpha * area_true;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_texpA);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_texpA_true);
        
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("1D thermal expansion B matrix")
    {
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
        
        Real first_area_moment_y;
        const MAST::FieldFunction<Real>& Ay = section.Ay();
        Ay(point, time, first_area_moment_y);
        CHECK( first_area_moment_y == Approx(first_area_moment_y_true) );
        
        Real first_area_moment_z;
        const MAST::FieldFunction<Real>& Az = section.Az();
        Az(point, time, first_area_moment_z);
        CHECK( first_area_moment_z == Approx(first_area_moment_z_true) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> texp_B_mat = section.thermal_expansion_B_matrix();
        
        RealMatrixX D_sec_texpB;
        texp_B_mat->operator()(point, time, D_sec_texpB);
        
        libMesh::out << "texp_B_mat =\n" << D_sec_texpB << std::endl;
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_texpB_true = RealMatrixX::Zero(2,1);
        D_sec_texpB_true(0,0) = cE * calpha * first_area_moment_z;
        D_sec_texpB_true(1,0) = cE * calpha * first_area_moment_y;
        
        
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


TEST_CASE("tube_element_property_card_constant_dynamic_1d",
          "[dynamic],[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real cE = 72.0e9;
    const Real cnu = 0.33;
    const Real cG = cE / (2.0 * (1.0+cnu) );
    const Real calpha = 5.43e-05;
    const Real crho = 1420.5;
    const Real ccp = 908.0;
    const Real ck = 237.0;
    const Real ckappa = 5.3284279639775167e-01;
    
    const Real coff_y = 0.35;
    const Real coff_z = 0.26;

    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param",            cE);     // Modulus of Elasticity
    MAST::Parameter nu("nu_param",          cnu);    // Poisson's ratio
    MAST::Parameter alpha("alpha_param",    calpha); // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param",        crho);   // Density
    MAST::Parameter cp("cp_param",          ccp);    // Specific Heat Capacity
    MAST::Parameter k("k_param",            ck);     // Thermal Conductivity

    // Define Section Properties as MAST Parameters
    MAST::Parameter r_o("DIM1", 1.125);   // Outer radius
    MAST::Parameter r_i("DIM2", 0.750);   // Inner radius
    MAST::Parameter offset_y("offy_param", coff_y);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", coff_z);     // Section offset in z-direction

    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", r_o);
    MAST::ConstantFieldFunction DIM2_f("DIM2", r_i);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);

    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;    

    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);

    // Initialize the section
    MAST::Solid1DTubeSectionElementPropertyCard section;


    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);

    // Add the material card to the section card
    section.set_material(material);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Now initialize the section
    section.init(*p_global_init);
    
    // True values
    const Real ro = r_o();
    const Real ri = r_i();
    const Real area_true = PI*(pow(ro,2.0) - pow(ri,2.0));
    const Real first_area_moment_z_true = area_true * coff_y;
    const Real first_area_moment_y_true = area_true * coff_z;
    const Real Izzc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Iyyc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Izyc_true = 0.0;
    const Real Ipc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/2.0;
    const Real second_area_moment_zz_true = Izzc_true + area_true * coff_y * coff_y;
    const Real second_area_moment_yy_true = Iyyc_true + area_true * coff_z * coff_z;
    const Real second_area_moment_zy_true = Izyc_true + area_true * coff_y * coff_z;
    const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
    
    const Real torsion_constant_true = 2.0184430643017226e+00;
    const Real warping_constant_true = 2.7224768298299495e-11;
    const Real kappa_z_true = 5.3284267813657360e-01;
    const Real kappa_y_true = 5.3284267085782810e-01;
    const Real xs_true = coff_z;
    const Real ys_true = coff_y;
    const Real xc_true = coff_z;
    const Real yc_true = coff_y;
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    SECTION("1D section inertia matrix")
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
        
        RealMatrixX D_sec_iner;
        inertia_mat->operator()(point, time, D_sec_iner);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_iner_true = RealMatrixX::Zero(6,6);
        
        Real area;
        const MAST::FieldFunction<Real>& Area = section.A();
        Area(point, time, area);
        REQUIRE( area == Approx(area_true) );
        D_sec_iner_true(0,0) = D_sec_iner_true(1,1) = D_sec_iner_true(2,2) = area_true;
        
        Real Ip;
        const MAST::FieldFunction<Real>& PolarInertia = section.Ip();
        PolarInertia(point, time, Ip);
        CHECK( Ip == Approx(second_area_moment_polar_true) );
        D_sec_iner_true(3,3) = second_area_moment_polar_true;
        
        Real first_area_moment_y;
        const MAST::FieldFunction<Real>& Ay = section.Ay();
        Ay(point, time, first_area_moment_y);
        CHECK( first_area_moment_y == Approx(first_area_moment_y_true) );
        D_sec_iner_true(0,4) = D_sec_iner_true(4,0) = first_area_moment_y_true;
        
        Real first_area_moment_z;
        const MAST::FieldFunction<Real>& Az = section.Az();
        Az(point, time, first_area_moment_z);
        CHECK( first_area_moment_z == Approx(first_area_moment_z_true) );
        D_sec_iner_true(0,5) = D_sec_iner_true(5,0) = first_area_moment_z_true;
        
        RealMatrixX I;
        const MAST::FieldFunction<RealMatrixX>& Inertias = section.I();
        Inertias(point, time, I);
        REQUIRE( I(0,1) == I(1,0) );
        Real Iyy = I(1,1);
        Real Izz = I(0,0);
        Real Izy = I(0,1);
        REQUIRE( Izz == Approx(second_area_moment_zz_true) );
        REQUIRE( Iyy == Approx(second_area_moment_yy_true) );
        REQUIRE( Izy == Approx(second_area_moment_zy_true) );

        D_sec_iner_true(4,4) = Iyy; // TODO: Should this be Izz?
        D_sec_iner_true(4,5) = Izy;
        D_sec_iner_true(5,4) = Izy;
        D_sec_iner_true(5,5) = Izz; // TODO: Should this be Iyy??
        
        D_sec_iner_true *= rho();

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


TEST_CASE("tube_element_property_card_constant_structural_1d",
          "[structural],[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real cE = 72.0e9;
    const Real cnu = 0.33;
    const Real cG = cE / (2.0 * (1.0+cnu) );
    const Real calpha = 5.43e-05;
    const Real crho = 1420.5;
    const Real ccp = 908.0;
    const Real ck = 237.0;
    const Real ckappa = 5.3284279639775167e-01;
    
    const Real coff_y = 0.35;
    const Real coff_z = 0.26;

    // Define Material Properties as MAST Parameters
    MAST::Parameter E("E_param",            cE);     // Modulus of Elasticity
    MAST::Parameter nu("nu_param",          cnu);    // Poisson's ratio
    MAST::Parameter alpha("alpha_param",    calpha); // Coefficient of thermal expansion
    MAST::Parameter rho("rho_param",        crho);   // Density
    MAST::Parameter cp("cp_param",          ccp);    // Specific Heat Capacity
    MAST::Parameter k("k_param",            ck);     // Thermal Conductivity

    // Define Section Properties as MAST Parameters
    MAST::Parameter r_o("DIM1", 1.125);   // Outer radius
    MAST::Parameter r_i("DIM2", 0.750);   // Inner radius
    MAST::Parameter offset_y("offy_param", coff_y);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", coff_z);     // Section offset in z-direction

    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction DIM1_f("DIM1", r_o);
    MAST::ConstantFieldFunction DIM2_f("DIM2", r_i);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction alpha_f("alpha_expansion", alpha);
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction cp_f("cp", cp);
    MAST::ConstantFieldFunction k_f("k_th", k);

    // Initialize the material
    MAST::IsotropicMaterialPropertyCard material;    

    // Add the material property constant field functions to the material card
    material.add(rho_f);
    material.add(k_f);
    material.add(cp_f);
    material.add(E_f);
    material.add(nu_f);
    material.add(alpha_f);

    // Initialize the section
    MAST::Solid1DTubeSectionElementPropertyCard section;


    // Add the section property constant field functions to the section card
    section.add(DIM1_f);
    section.add(DIM2_f);
    section.add(offsety_f);
    section.add(offsetz_f);

    // Add the material card to the section card
    section.set_material(material);

    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;

    // Now initialize the section
    section.init(*p_global_init);
    
    // True values
    const Real ro = r_o();
    const Real ri = r_i();
    const Real area_true = PI*(pow(ro,2.0) - pow(ri,2.0));
    const Real first_area_moment_z_true = area_true * coff_y;
    const Real first_area_moment_y_true = area_true * coff_z;
    const Real Izzc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Iyyc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/4.0;
    const Real Izyc_true = 0.0;
    const Real Ipc_true = PI*(pow(ro,4.0) - pow(ri,4.0))/2.0;
    const Real second_area_moment_zz_true = Izzc_true + area_true * coff_y * coff_y;
    const Real second_area_moment_yy_true = Iyyc_true + area_true * coff_z * coff_z;
    const Real second_area_moment_zy_true = Izyc_true + area_true * coff_y * coff_z;
    const Real second_area_moment_polar_true = second_area_moment_zz_true + second_area_moment_yy_true;
    
    const Real torsion_constant_true = 2.0184430643017226e+00;
    const Real warping_constant_true = 2.7224768298299495e-11;
    const Real kappa_z_true = 5.3284267813657360e-01;
    const Real kappa_y_true = 5.3284267085782810e-01;
    const Real xs_true = coff_z;
    const Real ys_true = coff_y;
    const Real xc_true = coff_z;
    const Real yc_true = coff_y;
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    
    SECTION("set_get_bending_model")
    {
        // NOTE: MAST::DKT and MAST::MINDLIN are not valid options for 1D sections, even though their input is accepted.
        section.set_bending_model(MAST::BERNOULLI);
        section.set_bending_model(MAST::DEFAULT_BENDING);
        section.set_bending_model(MAST::NO_BENDING);
        section.set_bending_model(MAST::TIMOSHENKO);
        
        // TODO: Implement element creation for testing of getting bending_model and checking default
        //REQUIRE( section.bending_model()
    }
    
//     SECTION("quadrature_order")
//     {
//         section.set_bending_model(MAST::MINDLIN);
//         REQUIRE( section.extra_quadrature_order(elem) == 0 );
//         
//         section.set_bending_model(MAST::DKT);
//         REQUIRE( section.extra_quadrature_order(elem) == 2 );
//     }
    
    SECTION("1D extension stiffness matrix")
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
        Real area;
        const MAST::FieldFunction<Real>& Area = section.A();
        Area(point, time, area);
        REQUIRE( area == Approx(area_true) );
        
        Real torsion_constant;
        const MAST::FieldFunction<Real>& TorsionConstant = section.J();
        TorsionConstant(point, time, torsion_constant);
        REQUIRE( torsion_constant == Approx(torsion_constant_true).epsilon(0.05) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> extension_stiffness_mat = section.stiffness_A_matrix();
        
        RealMatrixX D_sec_ext;
        extension_stiffness_mat->operator()(point, time, D_sec_ext);
        
        libMesh::out << "D_sec_ext\n" << D_sec_ext << std::endl;
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_ext_true = RealMatrixX::Zero(2,2);
        D_sec_ext_true(0,0) = cE * area_true;
        D_sec_ext_true(1,1) = cG * torsion_constant_true;
        
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_ext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_ext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth).epsilon(0.05) );
    }
    
    
    SECTION("1D bending section stiffness matrix")
    {
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
        RealMatrixX I;
        const MAST::FieldFunction<RealMatrixX>& Inertias = section.I();
        Inertias(point, time, I);
        REQUIRE( I(0,1) == I(1,0) );
        Real Iyy = I(1,1);
        Real Izz = I(0,0);
        Real Izy = I(0,1);
        REQUIRE( Izz == Approx(second_area_moment_zz_true) );
        REQUIRE( Iyy == Approx(second_area_moment_yy_true) );
        REQUIRE( Izy == Approx(second_area_moment_zy_true) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bending_stiffness_mat = section.stiffness_D_matrix();
        
        RealMatrixX D_sec_bnd;
        bending_stiffness_mat->operator()(point, time, D_sec_bnd);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bnd_true = RealMatrixX::Zero(2,2);
        D_sec_bnd_true(0,0) = cE * second_area_moment_zz_true;
        D_sec_bnd_true(1,1) = cE * second_area_moment_yy_true;
        D_sec_bnd_true(0,1) = cE * second_area_moment_zy_true;
        D_sec_bnd_true(1,0) = cE * second_area_moment_zy_true;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bnd);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bnd_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("1D extension-bending section stiffness matrix")
    {
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
        Real first_area_moment_y;
        const MAST::FieldFunction<Real>& Ay = section.Ay();
        Ay(point, time, first_area_moment_y);
        CHECK( first_area_moment_y == Approx(first_area_moment_y_true) );
        
        Real first_area_moment_z;
        const MAST::FieldFunction<Real>& Az = section.Az();
        Az(point, time, first_area_moment_z);
        CHECK( first_area_moment_z == Approx(first_area_moment_z_true) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> bndext_stiffness_mat = section.stiffness_B_matrix();
        
        RealMatrixX D_sec_bndext;
        bndext_stiffness_mat->operator()(point, time, D_sec_bndext);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_bndext_true = RealMatrixX::Zero(2,2);
        D_sec_bndext_true(0,0) = cE * first_area_moment_z_true;
        D_sec_bndext_true(0,1) = cE * first_area_moment_y_true;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_bndext);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_bndext_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth) );
    }
    
    SECTION("1D transverse shear section stiffness matrix")
    {
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
        Real area;
        const MAST::FieldFunction<Real>& Area = section.A();
        Area(point, time, area);
        REQUIRE( area == Approx(area_true) );
        
        std::unique_ptr<MAST::FieldFunction<RealMatrixX>> trans_shear_stiffness_mat = section.transverse_shear_stiffness_matrix();
        
        RealMatrixX D_sec_shr;
        trans_shear_stiffness_mat->operator()(point, time, D_sec_shr);
        
        // Hard-coded value of the section's extension stiffness
        RealMatrixX D_sec_shr_true = RealMatrixX::Zero(2,2);
        D_sec_shr_true(0,0) = cG * kappa_z_true * area_true;
        D_sec_shr_true(1,1) = cG * kappa_y_true * area_true;
        
        // Convert the test and truth Eigen::Matrix objects to std::vector
        // since Catch2 has built in methods to compare vectors
        std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_shr);
        std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_shr_true);
        
        // Floating point approximations are diffcult to compare since the
        // values typically aren't exactly equal due to numerical error.
        // Therefore, we use the Approx comparison instead of Equals
        CHECK_THAT( test, Catch::Approx<double>(truth).epsilon(0.01) );
    }
    
    
//     SECTION("1D spring section stiffness matrix")
//     {
//         /** 
//          * As of Dec. 17, 2019, spring_stiffness_matrix requires the input of an
//          * MAST::ElementBase object, but does not actually use it. To get 
//          * around requiring the creation of such an object (and therefore 
//          * the creation of a MAST::SystemInitialization, MAST::AssemblyBase,
//          * and MAST::GeomElem objects as well).
//          * 
//          * To remedy this, an additional method was added to MAST which allows
//          * transverse_shear_stiffness_matrix to be obtained without any input arguments.
//          */
//         std::unique_ptr<MAST::FieldFunction<RealMatrixX>> spring_stiffness_mat = section.stiffness_S_matrix();
//         
//         RealMatrixX D_sec_spring;
//         spring_stiffness_mat->operator()(point, time, D_sec_spring);
//         
//         // Hard-coded value of the section's extension stiffness
//         // NOTE: Should be all zero's for non-bushing sections
//         RealMatrixX D_sec_spring_true = RealMatrixX::Zero(4,6);
//         
//         // Convert the test and truth Eigen::Matrix objects to std::vector
//         // since Catch2 has built in methods to compare vectors
//         std::vector<double> test =  eigen_matrix_to_std_vector(D_sec_spring);
//         std::vector<double> truth = eigen_matrix_to_std_vector(D_sec_spring_true);
//         
//         // Floating point approximations are diffcult to compare since the
//         // values typically aren't exactly equal due to numerical error.
//         // Therefore, we use the Approx comparison instead of Equals
//         CHECK_THAT( test, Catch::Approx<double>(truth) );
//     }
}
