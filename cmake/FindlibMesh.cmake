# This module relies on libMesh_DIR being set.
#
# libMesh_FOUND - system has libMesh
# libMesh_INCLUDE_DIRS - libMesh include directories.
# libMesh_dbg_LIBRARIES - libMesh libraries (libmesh_dbg)
# libMesh_opt_LIBRARIES - libMesh libraries (libmesh_opt)

# Find the headers.
find_path(libMesh_INCLUDE_DIR libmesh/libmesh_config.h
          HINTS ${libMesh_DIR}/include)

# Find the optimized libraries.
find_library(libMesh_opt_LIBRARY
             NAMES mesh_opt
             HINTS ${libMesh_DIR}/lib)

# Find the debug libraries.
find_library(libMesh_dbg_LIBRARY
             NAMES mesh_dbg
             HINTS ${libMesh_DIR}/lib)

# If debug library is not available then set it to the optimized library
if(NOT libMesh_dbg_LIBRARY)
   message("Did not fine libmesh_dbg using libmesh_opt for debug version.")
   set(libMesh_dbg_LIBRARY  ${libMesh_opt_LIBRARY})
endif()


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
        REQUIRED_VARS libMesh_dbg_LIBRARY libMesh_opt_LIBRARY libMesh_INCLUDE_DIR
        VERSION_VAR libMesh_VERSION)

mark_as_advanced(libMesh_INCLUDE_DIR libMesh_LIBRARY libMesh_VERSION)

set(libMesh_dbg_LIBRARIES ${libMesh_dbg_LIBRARY})
set(libMesh_opt_LIBRARIES ${libMesh_opt_LIBRARY})
set(libMesh_INCLUDE_DIRS ${libMesh_INCLUDE_DIR})
