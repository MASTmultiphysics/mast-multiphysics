target_sources(mast_catch_tests
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_structural_shape_functions.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_strain_displacement_matrix.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_bending_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_bending_shear_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_extension_bending_coupling_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_internal_jacobian.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mast_quad4_linear_structural_inertial_jacobian.cpp)
               
add_test(NAME Element_Quad4_Structural_Shape_Functions
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_structural_shape_functions")
set_tests_properties(Element_Quad4_Structural_Shape_Functions
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_SETUP     Element_Quad4_Structural_Shape_Functions)

add_test(NAME Element_Quad4_Structural_Shape_Functions_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_structural_shape_functions")
set_tests_properties(Element_Quad4_Structural_Shape_Functions_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_SETUP     Element_Quad4_Structural_Shape_Functions_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Strain_Displacement_Matrix
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_strain_displacement_matrix")
set_tests_properties(Element_Quad4_Linear_Structural_Strain_Displacement_Matrix
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Quad4_Structural_Shape_Functions
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Strain_Displacement_Matrix)

add_test(NAME Element_Quad4_Linear_Structural_Strain_Displacement_Matrix_mpi
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_strain_displacement_matrix")
set_tests_properties(Element_Quad4_Linear_Structural_Strain_Displacement_Matrix_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Quad4_Structural_Shape_Functions_mpi
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Strain_Displacement_Matrix_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Extension
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_2D_Structural_Basic_Tests
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension)

add_test(NAME Element_Quad4_Linear_Structural_Extension_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_2D_Structural_Basic_Tests_mpi
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending)

add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_mpi
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Shear
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_shear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Shear
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Shear)

add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Shear_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_shear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Shear_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending_mpi
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Shear_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Coupling
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Coupling
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Coupling)

add_test(NAME Element_Quad4_Linear_Structural_Extension_Bending_Coupling_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_extension_bending_coupling_structural")
set_tests_properties(Element_Quad4_Linear_Structural_Extension_Bending_Coupling_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED  Element_Quad4_Linear_Structural_Extension_Bending_mpi
        FIXTURES_SETUP     Element_Quad4_Linear_Structural_Extension_Bending_Coupling_mpi)



add_test(NAME Element_Quad4_Linear_Structural_All
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_All
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED "Element_Quad4_Linear_Structural_Extension_Bending_Shear;Element_Quad4_Linear_Structural_Extension_Bending_Coupling")

add_test(NAME Element_Quad4_Linear_Structural_All_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural")
set_tests_properties(Element_Quad4_Linear_Structural_All_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED "Element_Quad4_Linear_Structural_Extension_Bending_Shear_mpi;Element_Quad4_Linear_Structural_Extension_Bending_Coupling_mpi")



add_test(NAME Element_Quad4_Linear_Structural_Inertial_Consistent
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_inertial_consistent")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Consistent
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests)

add_test(NAME Element_Quad4_Linear_Structural_Inertial_Consistent_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_inertial_consistent")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Consistent_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests_mpi)



add_test(NAME Element_Quad4_Linear_Structural_Inertial_Lumped
    COMMAND $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_inertial_lumped")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Lumped
    PROPERTIES
        LABELS "SEQ"
        FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests)

add_test(NAME Element_Quad4_Linear_Structural_Inertial_Lumped_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} -np 2 $<TARGET_FILE:mast_catch_tests> -w NoTests "quad4_linear_structural_inertial_lumped")
set_tests_properties(Element_Quad4_Linear_Structural_Inertial_Lumped_mpi
    PROPERTIES
        LABELS "MPI"
        FIXTURES_REQUIRED Element_2D_Structural_Basic_Tests_mpi)
