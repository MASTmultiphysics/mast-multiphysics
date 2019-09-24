# This module relies on DOT_DIR being set.
#
# DOT_FOUND - system has DOT.
# DOT_LIBRARIES - DOT libraries.

# Find the libraries.
find_library(DOT_LIBRARY
            NAMES dot
            HINTS ${DOT_DIR}/lib)

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DOT
REQUIRED_VARS DOT_LIBRARY)

if (DOT_FOUND)
    set (MAST_ENABLE_DOT 1)
else()
    set (MAST_ENABLE_DOT 0)
endif()

mark_as_advanced(DOT_LIBRARY DOT_FOUND MAST_ENABLE_DOT)

set(DOT_LIBRARIES ${DOT_LIBRARY})

