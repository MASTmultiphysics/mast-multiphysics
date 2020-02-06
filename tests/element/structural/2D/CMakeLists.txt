target_sources(mast_catch_tests
               PRIVATE
               ${CMAKE_CURRENT_LIST_DIR}/mast_structural_element_2d.cpp)

#add_subdirectory(quad4)

#2D Structural Element Test Basic Tests
add_test(NAME Element_2D_Structural_Basic_Tests
         COMMAND mast_catch_tests "structural_element_2d_base_tests")
set_tests_properties(Element_2D_Structural_Basic_Tests
                     PROPERTIES 
                     FIXTURES_REQUIRED  "Element_Property_Card_2D_Structural;libMesh_Mesh_Generation_2d"
                     FIXTURES_SETUP     Element_2D_Structural_Basic_Tests)
