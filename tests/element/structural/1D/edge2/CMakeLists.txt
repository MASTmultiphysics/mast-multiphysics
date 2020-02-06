target_sources( mast_catch_tests
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_internal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_internal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_coupling_internal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_extension_bending_shear_internal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_internal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_inertial_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_linear_structural_thermal_jacobian.cpp
                ${CMAKE_CURRENT_LIST_DIR}/mast_edge2_nonlinear_structural_thermal_jacobian.cpp)

#1D Structural Element Test Basic Tests
add_test(NAME Element_Edge2_Linear_Structural_Extension
         COMMAND mast_catch_tests "edge2_linear_extension_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_1D_Structural_Basic_Tests
                     FIXTURES_SETUP  Element_Edge2_Linear_Structural_Extension)
                     
add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending
         COMMAND mast_catch_tests "edge2_linear_extension_bending_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension
                     FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending)

add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Coupling
COMMAND mast_catch_tests "edge2_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Coupling
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending
                     FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Coupling)
                     
add_test(NAME Element_Edge2_Linear_Structural_Extension_Bending_Shear
COMMAND mast_catch_tests "edge2_linear_extension_bending_shear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_Extension_Bending_Shear
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Edge2_Linear_Structural_Extension_Bending
                     FIXTURES_SETUP     Element_Edge2_Linear_Structural_Extension_Bending_Shear)
                     
add_test(NAME Element_Edge2_Linear_Structural_All
         COMMAND mast_catch_tests "edge2_linear_structural")
set_tests_properties(Element_Edge2_Linear_Structural_All
                     PROPERTIES 
                     FIXTURES_REQUIRED "Element_Edge2_Linear_Structural_Extension_Bending_Coupling;Element_Edge2_Linear_Structural_Extension_Bending_Shear"
                     )
                     
add_test(NAME Element_Edge2_Linear_Structural_Inertial_Consistent
         COMMAND mast_catch_tests "edge2_linear_structural_inertial_consistent")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Consistent
                     PROPERTIES 
                     FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests)
                     
add_test(NAME Element_Edge2_Linear_Structural_Inertial_Lumped
         COMMAND mast_catch_tests "edge2_linear_structural_inertial_lumped")
set_tests_properties(Element_Edge2_Linear_Structural_Inertial_Lumped
                     PROPERTIES 
                     FIXTURES_REQUIRED Element_1D_Structural_Basic_Tests)
                     
add_test(NAME Element_Edge2_Linear_Structural_Thermal
         COMMAND mast_catch_tests "edge2_linear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Linear_Structural_Thermal
                     PROPERTIES 
                     FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load"
                     )
                     
add_test(NAME Element_Edge2_Nonlinear_Structural_Thermal
         COMMAND mast_catch_tests "edge2_nonlinear_structural_thermal_jacobian")
set_tests_properties(Element_Edge2_Nonlinear_Structural_Thermal
                     PROPERTIES 
                     FIXTURES_REQUIRED "Element_1D_Structural_Basic_Tests;Thermoelastic_Load"
                     )
