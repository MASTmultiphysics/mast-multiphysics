# This module relies on GCMMA_DIR being set.
#
# GCMMA_FOUND - system has GCMMA.
# GCMMA_LIBRARIES - GCMMA libraries.

# Find the libraries.
find_library(GCMMA_LIBRARY
             NAMES gcmma
             HINTS ${GCMMA_DIR}/lib)

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GCMMA
                                  REQUIRED_VARS GCMMA_LIBRARY)

if (GCMMA_FOUND)
  set (MAST_ENABLE_GCMMA 1)
else()
  set (MAST_ENABLE_GCMMA 0)
endif()

mark_as_advanced(GCMMA_LIBRARY GCMMA_FOUND MAST_ENABLE_GCMMA)

set(GCMMA_LIBRARIES ${GCMMA_LIBRARY})

