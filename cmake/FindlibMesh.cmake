# This module relies on libMesh_DIR being set.
#
# libMesh_FOUND - system has libMesh
# libMesh_INCLUDE_DIRS - libMesh include directories.
# libMesh_LIBRARIES - libMesh libraries (libmesh_opt)

# Find the headers.
find_path(libMesh_INCLUDE_DIR libmesh/libmesh_config.h
          HINTS ${libMesh_DIR}/include)

# Find the libraries.
find_library(libMesh_LIBRARY
             NAMES mesh_opt
             HINTS ${libMesh_DIR}/lib)

# Find libMesh version.
if(libMesh_INCLUDE_DIR)
    set(HEADER "${libMesh_INCLUDE_DIR}/libmesh/libmesh_config.h")
    file(STRINGS "${HEADER}" major REGEX "define +LIBMESH_MAJOR_VERSION")
    file(STRINGS "${HEADER}" minor REGEX "define +LIBMESH_MINOR_VERSION")
    file(STRINGS "${HEADER}" patch REGEX "define +LIBMESH_MICRO_VERSION")
    string(REGEX REPLACE ".+([0-9]+)" "\\1" major ${major})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" minor ${minor})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" patch ${patch})
    string(STRIP "${major}" major)
    string(STRIP "${minor}" minor)
    string(STRIP "${patch}" patch)
    set(libMesh_VERSION "${major}.${minor}.${patch}")
endif()

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libMesh
        REQUIRED_VARS libMesh_LIBRARY libMesh_INCLUDE_DIR
        VERSION_VAR libMesh_VERSION)

mark_as_advanced(libMesh_INCLUDE_DIR libMesh_LIBRARY libMesh_VERSION)

set(libMesh_LIBRARIES ${libMesh_LIBRARY})
set(libMesh_INCLUDE_DIRS ${libMesh_INCLUDE_DIR})
