# This module relies on libMesh_DIR being set.
#
# libMesh_FOUND - system has libMesh
# libMesh_INCLUDE_DIRS - libMesh include directories.
# libMesh_dbg_LIBRARIES - libMesh libraries (libmesh_dbg)
# libMesh_dev_LIBRARIES - libMesh libraries (libmesh_devel)
# libMesh_opt_LIBRARIES - libMesh libraries (libmesh_opt)

# Find the headers.
find_path(libMesh_INCLUDE_DIR libmesh/libmesh_config.h
    HINTS ${libMesh_DIR}/include)

# Find the optimized libraries.
find_library(libMesh_opt_LIBRARY
     NAMES mesh_opt
     HINTS ${libMesh_DIR}/lib)

find_library(timpi_opt_LIBRARY
     NAMES timpi_opt
     HINTS ${libMesh_DIR}/lib)

# Find the development libraries.
find_library(libMesh_dev_LIBRARY
     NAMES mesh_devel
     HINTS ${libMesh_DIR}/lib)

find_library(timpi_dev_LIBRARY
     NAMES timpi_devel
     HINTS ${libMesh_DIR}/lib)

# Find the debug libraries.
find_library(libMesh_dbg_LIBRARY
     NAMES mesh_dbg
     HINTS ${libMesh_DIR}/lib)

find_library(timpi_dbg_LIBRARY
     NAMES timpi_dbg
     HINTS ${libMesh_DIR}/lib)
            
# If debug library is not available, then set it to the development library if available, else set it to the optimized library.
if(NOT libMesh_dbg_LIBRARY)
    if(NOT libMesh_dev_LIBRARY)
        if(libMesh_opt_LIBRARY)
            message(WARNING "Did not find libmesh_dbg or libmesh_devel, using libmesh_opt for debug version.")
            SET(libMesh_dbg_LIBRARY ${libMesh_opt_LIBRARY})
        else()
            message(FATAL_ERROR "Could not find any libMesh libraries.")
        endif()
    else()
        message(WARNING "Did not find libmesh_dbg using libmesh_devel for debug version.")
        SET(libMesh_dbg_LIBRARY ${libMesh_dev_LIBRARY})
    endif()
endif()

# If development library is not available, then set it to the debug library if available, else set it to the optimized library.
if(NOT libMesh_dev_LIBRARY)
    if(NOT libMesh_dbg_LIBRARY)
        message(WARNING "Did not find libMesh_devel libmesh_dbg, using libmesh_opt for development version.")
        SET(libMesh_dev_LIBRARY ${libMesh_opt_LIBRARY})
    else()
        message(WARNING "Did not find libMesh_devel, using libmesh_dbg for devel versions.")
        SET(libMesh_dev_LIBRARY ${libMesh_dbg_LIBRARY})
    endif()
endif()

# If optimized library is not available, then set it to the development library if available, else set it to the debug library.
if(NOT libMesh_opt_LIBRARY)
    if(NOT libMesh_devel_LIBRARY)
        message(WARNING "Did not find libmesh_opt or libmesh_devel, using libmesh_dbg for optimized version.")
        SET(libMesh_opt_LIBRARY ${libMesh_dbg_LIBRARY})
    else()
        message(WARNING "Did not find libmesh_opt, using libmesh_devel for optimized version.")
        SET(libMesh_opt_LIBRARY ${libMesh_dev_LIBRARY})
    endif()
endif()

if (NOT timpi_opt_LIBRARY)
    SET(timpi_opt_LIBRARY ${libMesh_opt_LIBRARY})
endif()

if (NOT timpi_dev_LIBRARY)
    SET(timpi_dev_LIBRARY ${libMesh_dev_LIBRARY})
endif()

if (NOT timpi_dbg_LIBRARY)
    SET(timpi_dbg_LIBRARY ${libMesh_dbg_LIBRARY})
endif()

# Find libMesh version.
if(libMesh_INCLUDE_DIR)
    set(HEADER "${libMesh_INCLUDE_DIR}/libmesh/libmesh_config.h")
    file(STRINGS "${HEADER}" major REGEX "define +LIBMESH_MAJOR_VERSION")
    file(STRINGS "${HEADER}" minor REGEX "define +LIBMESH_MINOR_VERSION")
    file(STRINGS "${HEADER}" patch REGEX "define +LIBMESH_MICRO_VERSION")
    string(REGEX MATCH "#define LIBMESH_MAJOR_VERSION *([0-9]*)" _ ${major})
    set(major ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define LIBMESH_MINOR_VERSION *([0-9]*)" _ ${minor})
    set(minor ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define LIBMESH_MICRO_VERSION *([0-9]*)" _ ${patch})
    set(patch ${CMAKE_MATCH_1})
    #string(REGEX REPLACE ".+([0-9]+)" "\\1" major ${major})
    #string(REGEX REPLACE ".+([0-9]+)" "\\1" minor ${minor})
    #string(REGEX REPLACE ".+([0-9]+)" "\\1" patch ${patch})
    string(STRIP "${major}" major)
    string(STRIP "${minor}" minor)
    string(STRIP "${patch}" patch)
    set(libMesh_VERSION "${major}.${minor}.${patch}")
endif()

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libMesh
    REQUIRED_VARS
        libMesh_dbg_LIBRARY
        libMesh_dev_LIBRARY
        libMesh_opt_LIBRARY
        timpi_dbg_LIBRARY
        timpi_dev_LIBRARY
        timpi_opt_LIBRARY
        libMesh_INCLUDE_DIR
        VERSION_VAR libMesh_VERSION)

mark_as_advanced(libMesh_INCLUDE_DIR
                 libMesh_dbg_LIBRARY
                 libMesh_dev_LIBRARY
                 libMesh_opt_LIBRARY
                 timpi_dbg_LIBRARY
                 timpi_dev_LIBRARY
                 timpi_opt_LIBRARY
                 libMesh_VERSION
                 libMesh_FOUND)

set(libMesh_dbg_LIBRARIES ${libMesh_dbg_LIBRARY})
set(libMesh_dev_LIBRARIES ${libMesh_dev_LIBRARY})
set(libMesh_opt_LIBRARIES ${libMesh_opt_LIBRARY})
set(timpi_dbg_LIBRARIES ${timpi_dbg_LIBRARY})
set(timpi_dev_LIBRARIES ${timpi_dev_LIBRARY})
set(timpi_opt_LIBRARIES ${timpi_opt_LIBRARY})
set(libMesh_INCLUDE_DIRS ${libMesh_INCLUDE_DIR})
