# This CMake defines the MAST library target.
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
            ${PETSc_INCLUDE_DIRS}
            ${SLEPc_INCLUDE_DIRS}
            ${libMesh_INCLUDE_DIRS})

# Add libraries.
target_link_libraries(mast
        PUBLIC
            ${MPI_CXX_LIBRARIES}
            ${PETSc_LIBRARIES}
            ${SLEPc_LIBRARIES}
            ${libMesh_LIBRARIES})

# Add subdirectories containing sources.
add_subdirectory(base)
add_subdirectory(solver)

# NOTE: Use of PUBLIC keyword above means other CMake target (like an example)
#       that target_link_library(XXXX mast) will inheret these properties.