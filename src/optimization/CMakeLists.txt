target_sources(mast
        PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/function_evaluation.cpp
        ${CMAKE_CURRENT_LIST_DIR}/function_evaluation.h
        ${CMAKE_CURRENT_LIST_DIR}/optimization_interface.h)

if(${ENABLE_DOT})
    target_sources(mast
            PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/dot_optimization_interface.cpp
            ${CMAKE_CURRENT_LIST_DIR}/dot_optimization_interface.h)
    # TODO: Create FindDOT CMake module and add target_link_libraries(mast PRIVATE ${DOT_LIBRARIES})
endif()

if(${ENABLE_GCMMA})
    target_sources(mast
            PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/gcmma_optimization_interface.cpp
            ${CMAKE_CURRENT_LIST_DIR}/gcmma_optimization_interface.h)
    # TODO: Create FindGCMMA CMake module and add target_link_libraries(mast PRIVATE ${GCMMA_LIBRARIES})
endif()

if(${ENABLE_NPSOL})
    target_sources(mast
            PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/npsol_optimization_interface.cpp
            ${CMAKE_CURRENT_LIST_DIR}/npsol_optimization_interface.h)
    # TODO: Create FindNPSOL CMake module and add target_link_libraries(mast PRIVATE ${NPSOL_LIBRARIES})
endif()

# Note: In order to add the target_link_libraries(mast ...) inside in this directory, I think the add_subdirectory()
#       calls one directory up need to be switched to:
#
# include(aeroelasticity/CMakeLists.txt)
# include(base/CMakeLists.txt)
# include(boundary_condition/CMakeLists.txt)
# include(coordinates/CMakeLists.txt)
# include(elasticity/CMakeLists.txt)
# include(fluid/CMakeLists.txt)
# include(heat_conduction/CMakeLists.txt)
# include(level_set/CMakeLists.txt)
# include(mesh/CMakeLists.txt)
# include(numerics/CMakeLists.txt)
# include(optimization/CMakeLists.txt)
# include(property_cards/CMakeLists.txt)
# include(solver/CMakeLists.txt)
#
# See: https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/optimization
        FILES_MATCHING PATTERN "*.h")