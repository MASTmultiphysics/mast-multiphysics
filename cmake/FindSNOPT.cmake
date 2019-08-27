# This module relies on SNOPT_DIR being set.
#
# SNOPT_FOUND - system has SNPT.
# SNOPT_LIBRARIES - SNOPT libraries.

# Find the libraries.
find_library(SNOPT_LIBRARY
             NAMES snopt7
             HINTS ${SNOPT_DIR}/lib)

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SNOPT
                                  REQUIRED_VARS SNOPT_LIBRARY)

if (SNOPT_FOUND)
    set (MAST_ENABLE_SNOPT 1)
else()
    set (MAST_ENABLE_SNOPT 0)
endif()

mark_as_advanced(SNOPT_LIBRARY SNOPT_FOUND MAST_ENABLE_SNOPT)

