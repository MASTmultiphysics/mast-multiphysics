# This CMakeLists.txt defines the MAST library target (mast).
# - Source files are added to the `mast` target in subdirectories using the
#   target_sources(mast PRIVATE <source_files...>) command.
#
# - This has the advantage of keeping definitions localized.

# Define the target.
add_library(mast)

# Add include directories.
target_include_directories(mast
        PUBLIC
            ${CMAKE_CURRENT_LIST_DIR}
            ${MPI_CXX_INCLUDE_PATH}
            ${PETSc_INCLUDE_DIRS}
            ${SLEPc_INCLUDE_DIRS}
            ${libMesh_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS}) # Could possibly move this into src/optimization since it is closer to where Boost is actually used.

# Add libraries.
target_link_libraries(mast
        PUBLIC
            ${MPI_CXX_LIBRARIES}
            ${PETSc_LIBRARIES}
            ${SLEPc_LIBRARIES}
            ${libMesh_LIBRARIES})

# NOTE: Use of PUBLIC keyword above means other CMake target (like an example)
#       that target_link_library(XXXX mast) will inheret these properties.

# Add subdirectories containing sources.
add_subdirectory(aeroelasticity)
add_subdirectory(base)
add_subdirectory(boundary_condition)
add_subdirectory(coordinates)
add_subdirectory(elasticity)
add_subdirectory(fluid)
add_subdirectory(heat_conduction)
add_subdirectory(level_set)
add_subdirectory(mesh)
add_subdirectory(numerics)
add_subdirectory(optimization)
add_subdirectory(property_cards)
add_subdirectory(solver)

# NOTE: There is an argument to be made for why we might want to switch the
#       add_subdirectory(<dir>) above to include(<dir>/CMakeLists.txt).
#       - https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/
#
#       This would enable us to call target_link_libraries() closer to to the
#       code that is using it. (There is a documented limitation with CMake
#       that prevents this with using add_subdirectory).

# MAST library installation.
install(TARGETS mast
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})