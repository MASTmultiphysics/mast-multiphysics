# This module relies on NLOPT_DIR being set.
#
# NLOPT_FOUND - system has NLOPT.
# NLOPT_INCLUDE_DIRS - NLOPT include directories.
# NLOPT_LIBRARIES - NLOPT libraries.
# NLOPT_VERSION - NLOPT version.

# Find the headers.
find_path(NLOPT_INCLUDE_DIR nlopt.h
          HINTS ${NLOPT_DIR}/include)

# Find the libraries.
find_library(NLOPT_LIBRARY
             NAMES nlopt
             HINTS ${NLOPT_DIR}/lib)

# Find the version number.
if(EXISTS "${NLOPT_INCLUDE_DIR}/../lib64/pkgconfig/nlopt.pc")
    file(READ "${NLOPT_INCLUDE_DIR}/../lib64/pkgconfig/nlopt.pc" NLOPTVERINFO)
    string(REGEX MATCH "Version: ([A-Za-z0-9]+(\\.[A-Za-z0-9]+)+)" _ ${NLOPTVERINFO})
    set(NLOPT_VERSION ${CMAKE_MATCH_1})
elseif(EXISTS "${NLOPT_INCLUDE_DIR}/../lib64/cmake//nlopt/NLoptConfig.cmake")
    file(READ "${NLOPT_INCLUDE_DIR}/../lib64/cmake//nlopt/NLoptConfig.cmake" NLOPTVERINFO)
    string(REGEX MATCH "set \\(NLOPT_VERSION \"([A-Za-z0-9]+(\\.[A-Za-z0-9]+)+)\"\\)" _ ${NLOPTVERINFO})
    set(NLOPT_VERSION ${CMAKE_MATCH_1})
else()
    set(NLOPT_VERSION "")
endif()

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOPT
                                  REQUIRED_VARS NLOPT_LIBRARY NLOPT_INCLUDE_DIR
                                  VERSION_VAR NLOPT_VERSION)

if (NLOPT_FOUND)
    set (MAST_ENABLE_NLOPT 1)
else()
    set (MAST_ENABLE_NLOPT 0)
endif()

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY NLOPT_FOUND MAST_ENABLE_NLOPT)

set(NLOPT_LIBRARIES ${NLOPT_LIBRARY})
set(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR})
