# This module relies on NLOPT_DIR being set.
#
# NLOPT_FOUND - system has NLOPT.
# NLOPT_INCLUDE_DIRS - NLOPT include directories.
# NLOPT_LIBRARIES - NLOPT libraries.

# Find the headers.
find_path(NLOPT_INCLUDE_DIR nlopt.h
          HINTS ${NLOPT_DIR}/include)

# Find the libraries.
find_library(NLOPT_LIBRARY
             NAMES nlopt
             HINTS ${NLOPT_DIR}/lib)

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT
                                  REQUIRED_VARS NLOPT_LIBRARY NLOPT_INCLUDE_DIR)

if (NLOPT_FOUND)
    set (MAST_ENABLE_NLOPT 1)
else()
    set (MAST_ENABLE_NLOPT 0)
endif()

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY NLOPT_FOUND MAST_ENABLE_NLOPT)

set(NLOPT_LIBRARIES ${NLOPT_LIBRARY})
set(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR})
