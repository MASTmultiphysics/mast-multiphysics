target_sources(mast_catch_tests
               PRIVATE
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_structural_shape_functions.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_strain_displacement_matrix.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_internal_jacobian.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_bending_shear_internal_jacobian.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_bending_coupling_internal_jacobian.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_internal_jacobian.cpp
               ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_inertial_jacobian.cpp)
               
add_test(NAME Element_Quad4_Structural_Shape_Functions
         COMMAND mast_catch_tests "quad4_structural_shape_functions")
set_tests_properties(Element_Quad4_Structural_Shape_Functions
                     PROPERTIES 
                     FIXTURES_SETUP     Element_Quad4_Structural_Shape_Functions)
                     
add_test(NAME Element_Quad4_Linear_Structural_Strain_Displacement_Matrix
         COMMAND mast_catch_tests "quad4_linear_structural_strain_displacement_matrix")
set_tests_properties(Element_Quad4_Linear_Structural_Strain_Displacement_Matrix
                     PROPERTIES
                     FIXTURES_REQUIRED  Element_Quad4_Structural_Shape_Functions
                     FIXTURES_SETUP     Element_Quad4_Linear_Structural_Strain_Displacement_Matrix)

add_test(NAME Element_Quad4_Linear_Structural_Extension
         COMMAND mast_catch_tests "quad4_linear_extension_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_2D_Structural_Basic_Tests
                     FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension)
                     
add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending
         COMMAND mast_catch_tests "quad4_linear_extension_bending_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension
                     FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending)
                     
add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Shear
         COMMAND mast_catch_tests "quad4_linear_extension_bending_shear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Shear
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending
                     FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Shear)
                     
add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Coupling
         COMMAND mast_catch_tests "quad4_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Coupling
                     PROPERTIES 
                     FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending
                     FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Coupling)
                     
add_test(NAME Element_Quad4_Linear_Structural_All
         COMMAND mast_catch_tests "quad4_linear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_All
                     PROPERTIES 
                     FIXTURES_REQUIRED "Element_Quad4_Linear_Structural_Extension_Bending_Shear;Element_Quad4_Linear_Structural_Extension_Bending_Coupling"
                     )
                     
add_test(NAME Element_Quad4_Linear_Structural_Inertial_Consistent
         COMMAND mast_catch_tests "quad4_linear_structural_inertial_consistent")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Consistent
                     PROPERTIES
                     FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests)
                     
add_test(NAME Element_Quad4_Linear_Structural_Inertial_Lumped
         COMMAND mast_catch_tests "quad4_linear_structural_inertial_lumped")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Lumped
                     PROPERTIES
                     FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests)
