// Catch2 includes
#include "catch.hpp"

// libMesh includes
#include "libmesh/point.h"

// MAST includes
#include "base/parameter.h"
#include "base/constant_field_function.h"
#include "property_cards/isotropic_material_property_card.h"
#include "property_cards/solid_1d_arbitrary_section_element_property_card.h"

// Custom includes
#include "test_helpers.h"

extern libMesh::LibMeshInit* p_global_init;


TEST_CASE("arbitrary_element_property_card_constant_base_1d",
          "[1D],[isotropic],[constant],[property]")
{
    const uint dim = 1;
    
    const Real A_test = 11.3;
    const Real Iyyc_test = 3.1;
    const Real Izzc_test = 4.2;
    const Real Izyc_test = 2.6;
    const Real Ipc_test = Iyyc_test + Izzc_test;
    const Real J_test = 5.2;
    const Real W_test = 3.2;
    const Real kappa_zz_test = 0.66;
    const Real kappa_yy_test = 0.23;
    
    //const Real offset_y_test = 0.547;
    //const Real offset_z_test = -0.258;
    
    const Real offset_y_test = 0.0;
    const Real offset_z_test = 0.0;
    
    const Real Qz_true = A_test * offset_y_test;
    const Real Qy_true = A_test * offset_z_test;
    
    const Real Iyy_true = Iyyc_test + A_test * offset_z_test * offset_z_test;
    const Real Izz_true = Izzc_test + A_test * offset_y_test * offset_y_test;
    const Real Izy_true = Izyc_test + A_test * offset_y_test * offset_z_test;
    const Real Ip_true  = Iyy_true + Izz_true;
    
    // Define Material Properties as MAST Parameters
    MAST::Parameter rho("rho_param", 1420.5);         // Density
    MAST::Parameter E("E_param", 72.0e9);             // Modulus of Elasticity
    MAST::Parameter nu("nu_param", 0.33);             // Poisson's ratio
    MAST::Parameter kappa("kappa_param", 5.0/6.0);    // Shear coefficient
    MAST::Parameter cp("cp_param",   908.0);          // Specific Heat Capacity
    MAST::Parameter k("k_param",     237.0);          // Thermal Conductivity
    
    // Define Section Properties as MAST Parameters
    MAST::Parameter A("A", A_test);
    MAST::Parameter Iyyc("Iyyc", Iyyc_test);   
    MAST::Parameter Izzc("Izzc", Izzc_test);   
    MAST::Parameter Izyc("Izyc", Izyc_test);
    MAST::Parameter Ipc("Ipc", Ipc_test);
    MAST::Parameter J("J", J_test);
    MAST::Parameter W("W", W_test);
    MAST::Parameter Kzz("Kzz", kappa_zz_test);
    MAST::Parameter Kyy("Kyy", kappa_yy_test);
    
    MAST::Parameter offset_y("offy_param", offset_y_test);     // Section offset in y-direction
    MAST::Parameter offset_z("offz_param", offset_z_test);     // Section offset in z-direction
    
    // Create field functions to dsitribute these constant parameters throughout the model
    // Section Property Field Functions
    MAST::ConstantFieldFunction A_f("A", A);
    MAST::ConstantFieldFunction Iyy_f("Iyy", Iyyc);
    MAST::ConstantFieldFunction Izz_f("Izz", Izzc);
    MAST::ConstantFieldFunction Izy_f("Izy", Izyc);
    MAST::ConstantFieldFunction Ip_f("Ip", Ipc);
    MAST::ConstantFieldFunction J_f("J", J);
    MAST::ConstantFieldFunction W_f("W", W);
    MAST::ConstantFieldFunction Kzz_f("Kappazz", Kzz);
    MAST::ConstantFieldFunction Kyy_f("Kappayy", Kyy);
    MAST::ConstantFieldFunction offsety_f("hy_off", offset_y);
    MAST::ConstantFieldFunction offsetz_f("hz_off", offset_z);
    // Material Property Field Functions
    MAST::ConstantFieldFunction rho_f("rho", rho);
    MAST::ConstantFieldFunction E_f("E", E);
    MAST::ConstantFieldFunction nu_f("nu", nu);
    MAST::ConstantFieldFunction kappa_f("kappa", kappa);
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
    material.add(kappa_f);
    
    // Initialize the section
    MAST::Solid1DArbitrarySectionElementPropertyCard section;
    
    // Add the section property constant field functions to the section card
    section.add(A_f);
    section.add(Iyy_f);
    section.add(Izz_f);
    section.add(Izy_f);
    section.add(Ip_f);
    section.add(J_f);
    section.add(W_f);
    section.add(Kzz_f);
    section.add(Kyy_f);
    section.add(offsety_f);
    section.add(offsetz_f);
    
    // Add the material card to the section card
    section.set_material(material);
    
    // Specify a section orientation point and add it to the section.
    RealVectorX orientation = RealVectorX::Zero(3);
    orientation(1) = 1.0;
    section.y_vector() = orientation;
    
    // Now initialize the section
    section.init();
    
    const libMesh::Point point(4.3, -3.5, -6.7);
    const Real time = 8.22;
    const libMesh::Point Zp(0.0, 0.0, 0.0);
    
    SECTION("cross_sectional_properties")
    {
        REQUIRE( section.dim() == dim); // Ensure section is 1 dimensional
        REQUIRE( section.depends_on(A) );
        REQUIRE( section.depends_on(Iyyc) );
        REQUIRE( section.depends_on(Izzc) );
        REQUIRE( section.depends_on(Izyc) );
        REQUIRE( section.depends_on(Ipc) );
        REQUIRE( section.depends_on(J) );
        REQUIRE( section.depends_on(W) );
        REQUIRE( section.depends_on(Kzz) );
        REQUIRE( section.depends_on(Kyy) );
        REQUIRE( section.depends_on(offset_y) );
        REQUIRE( section.depends_on(offset_z) );
        REQUIRE( section.depends_on(k) );
        REQUIRE( section.depends_on(cp) );
        REQUIRE( section.depends_on(rho) );
        REQUIRE( section.if_isotropic() );
        
        Real area;
        const MAST::FieldFunction<Real>& Area = section.A();
        Area(point, time, area);
        REQUIRE( area == Approx(A_test) );
        
        Real first_area_moment_y;
        const MAST::FieldFunction<Real>& Ay = section.Ay();
        Ay(point, time, first_area_moment_y);
        CHECK( first_area_moment_y == Approx(Qy_true) );
        
        Real first_area_moment_z;
        const MAST::FieldFunction<Real>& Az = section.Az();
        Az(point, time, first_area_moment_z);
        CHECK( first_area_moment_z == Approx(Qz_true) );
        
        RealMatrixX I;
        const MAST::FieldFunction<RealMatrixX>& Inertias = section.I();
        Inertias(point, time, I);
        REQUIRE( I(0,1) == I(1,0) );
        Real Iyyv = I(1,1);
        Real Izzv = I(0,0);
        Real Izyv = I(0,1);
        REQUIRE( Izzv == Approx(Izz_true) );
        REQUIRE( Iyyv == Approx(Iyy_true) );
        REQUIRE( Izyv == Approx(Izy_true) );
        
        Real Ipv;
        const MAST::FieldFunction<Real>& PolarInertia = section.Ip();
        PolarInertia(point, time, Ipv);
        REQUIRE( Ipv == Approx(Ip_true) );
        
        Real torsion_constant;
        const MAST::FieldFunction<Real>& TorsionConstant = section.J();
        TorsionConstant(point, time, torsion_constant);
        REQUIRE( torsion_constant == Approx(J_test));
        
        Real warping_constant;
        const MAST::FieldFunction<Real>& WarpingConstant = section.Gam();
        WarpingConstant(point, time, warping_constant);
        REQUIRE( warping_constant == Approx(W_test) );
        
        RealMatrixX shear_coefficients;
        const MAST::FieldFunction<RealMatrixX>& ShearCoefficientMatrix = section.Kap();
        ShearCoefficientMatrix(point, time, shear_coefficients);
        REQUIRE( shear_coefficients(0,0) == Approx(kappa_zz_test) );
        REQUIRE( shear_coefficients(1,1) == Approx(kappa_yy_test) );
    }
    
    SECTION("diagonal_mass_matrix_flag")
    {
        REQUIRE_FALSE( section.if_diagonal_mass_matrix() );
    
        section.set_diagonal_mass_matrix(true);
        REQUIRE( section.if_diagonal_mass_matrix() );
    }
    
    SECTION("stress_points")
    {
        const libMesh::Point Cp(0.1, 2.0, 0.0);
        const libMesh::Point Dp(0.5, -3.1, 0.0);
        const libMesh::Point Ep(-0.6, -0.7, 0.0);
        const libMesh::Point Fp(-0.7, 1.1, 0.0);
        
        section.add_stress_point(Cp);
        section.add_stress_point(Dp);
        section.add_stress_point(Ep);
        section.add_stress_point(Fp);
        
        std::vector<libMesh::Point> stress_points = section.get_stress_points(point, time, Zp);
    
        REQUIRE( stress_points.size() == 4 );
        REQUIRE( stress_points[0] == Cp );
        REQUIRE( stress_points[1] == Dp );
        REQUIRE( stress_points[2] == Ep );
        REQUIRE( stress_points[3] == Fp );
        
        std::vector<libMesh::Point> dstress_points = section.get_stress_points_derivative(rho, point, time, Zp);
        for (uint i=0; i<dstress_points.size(); i++)
        {
            REQUIRE(dstress_points[i] == Zp );
        }
    }
    
    
    SECTION("centroid location")
    {
        libMesh::Point centroid = section.get_centroid(point, time);
        REQUIRE( centroid == Zp );
        
        libMesh::Point dcentroid = section.get_centroid_derivative(rho, point, time);
        REQUIRE( dcentroid == Zp );
    }
    
    
    SECTION("shear center location")
    {
        libMesh::Point shear_center = section.get_shear_center(point, time);
        REQUIRE( shear_center == Zp);
        
        libMesh::Point dshear_center = section.get_shear_center_derivative(rho, point, time);
        REQUIRE( dshear_center == Zp );
    }
}
