# This CMakeLists.txt defines the MAST library target (mast).
# - Source files are added to the `mast` target in subdirectories using the
#   target_sources(mast PRIVATE <source_files...>) command.
#
# - This has the advantage of keeping definitions localized.


# Define the libraries that MAST will link against
# (Also used for setting RPATH)
set(MAST_LINK_LIBS 
    ${MPI_CXX_LIBRARIES}
    ${GCMMA_LIBRARY}
    ${DOT_LIBRARY}
    ${NLOPT_LIBRARY}
    ${SNOPT_LIBRARY}
    ${LAPACK_LIBRARIES}
    ${PETSc_LIBRARIES}
    ${SLEPc_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},RelWithDebInfo>:${libMesh_dev_LIBRARIES}>
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},RelWithDebInfo>:${timpi_dev_LIBRARIES}>
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:${libMesh_dbg_LIBRARIES}>
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:${timpi_dbg_LIBRARIES}>
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},Release>:${libMesh_opt_LIBRARIES}>
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},Release>:${timpi_opt_LIBRARIES}>
    ${LIBGFORTRAN_LIBRARIES})


# Define the RPATH which will be set on the installed MAST library
if (USE_RPATH)
    set(CMAKE_INSTALL_RPATH "")
    message("Setting RPATH/RUNPATH...")
    foreach(library_full_path ${MAST_LINK_LIBS})
        get_filename_component(lib_dir "${library_full_path}" PATH)
        string(REPLACE "$<$<STREQUAL:${CMAKE_BUILD_TYPE},RelWithDebInfo>:" "" lib_dir "${lib_dir}")
        string(REPLACE "$<$<STREQUAL:${CMAKE_BUILD_TYPE},Debug>:" "" lib_dir "${lib_dir}")
        string(REPLACE "$<$<STREQUAL:${CMAKE_BUILD_TYPE},Release>:" "" lib_dir "${lib_dir}") 
        if(NOT "${lib_dir}" IN_LIST CMAKE_INSTALL_RPATH)
            message(STATUS "Adding ${lib_dir} to CMAKE_INSTALL_RPATH")
            list( APPEND CMAKE_INSTALL_RPATH "${lib_dir}")
        endif()
    endforeach()
endif()

# Define the target.
add_library(mast SHARED)

# Add include directories.
target_include_directories(mast
        PUBLIC
            ${CMAKE_SOURCE_DIR}/contrib
            ${CMAKE_CURRENT_LIST_DIR}
            ${MPI_CXX_INCLUDE_PATH}
            ${PETSc_INCLUDE_DIRS}
            ${SLEPc_INCLUDE_DIRS}
            ${HDF5_INCLUDE_DIRS}
            ${libMesh_INCLUDE_DIRS}
            ${EIGEN3_INCLUDE_DIR}
            ${NLOPT_INCLUDE_DIR}
            ${Boost_INCLUDE_DIRS}) # Could possibly move this into src/optimization since it is closer to where Boost is actually used.

# Add libraries.
target_link_libraries(mast
            PUBLIC
                ${MAST_LINK_LIBS})

# # Replaced with the short command above since the target link libraries are
# # defined int he MAST_LINK_LIBS list and the build type is handled using
# # generator expressions within that list.
# if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" OR "${USE_LIBMESH_DEVEL}")
#     message(STATUS "Attempting to use libmesh_devel library.")

#     target_link_libraries(mast
#             PUBLIC
#                 ${MPI_CXX_LIBRARIES}
#                 ${GCMMA_LIBRARY}
#                 ${DOT_LIBRARY}
#                 ${NLOPT_LIBRARY}
#                 ${SNOPT_LIBRARY}
#                 ${LAPACK_LIBRARIES}
#                 ${PETSc_LIBRARIES}
#                 ${SLEPc_LIBRARIES}
#                 ${HDF5_LIBRARIES}
#                 ${Boost_IOSTREAMS_LIBRARY}
#                 ${Boost_FILESYSTEM_LIBRARY}
#                 ${Boost_SYSTEM_LIBRARY}
#                 ${libMesh_dev_LIBRARIES}
#                 ${timpi_dev_LIBRARIES}
#                 ${LIBGFORTRAN_LIBRARIES})
# else()
#     if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
#         message(STATUS "Attempting to use libmesh_debug library.")
#     else()
#         message(STATUS "Attempting to use libmesh_opt library.")
#     endif()
#     target_link_libraries(mast
#             PUBLIC
#                 ${MPI_CXX_LIBRARIES}
#                 ${GCMMA_LIBRARY}
#                 ${DOT_LIBRARY}
#                 ${NLOPT_LIBRARY}
#                 ${SNOPT_LIBRARY}
#                 ${LAPACK_LIBRARIES}
#                 ${PETSc_LIBRARIES}
#                 ${SLEPc_LIBRARIES}
#                 ${HDF5_LIBRARIES}
#                 ${Boost_IOSTREAMS_LIBRARY}
#                 ${Boost_FILESYSTEM_LIBRARY}
#                 ${Boost_SYSTEM_LIBRARY}
#                 debug ${libMesh_dbg_LIBRARIES}
#                 debug ${timpi_dbg_LIBRARIES}
#                 optimized ${libMesh_opt_LIBRARIES}
#                 optimized ${timpi_opt_LIBRARIES}
#                 ${LIBGFORTRAN_LIBRARIES})
# endif()

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
add_subdirectory(utility)

# NOTE: There is an argument to be made for why we might want to switch the
#       add_subdirectory(<dir>) above to include(<dir>/CMakeLists.txt).
#       - https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/
#
#       This would enable us to call target_link_libraries() closer to to the
#       code that is using it. (There is a documented limitation with CMake
#       that prevents this with using add_subdirectory).

# Add third party contrib folder sources.
target_sources(mast
    PUBLIC
        ${CMAKE_SOURCE_DIR}/contrib/libfort/fort.c
        ${CMAKE_SOURCE_DIR}/contrib/libfort/fort.h
        ${CMAKE_SOURCE_DIR}/contrib/libfort/fort.hpp)

# MAST library installation.
install(TARGETS mast
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
