# This module relies on Eigen_DIR being set.
#
# Note that Eigen is a header only library.
#
# Eigen_FOUND - system has Eigen
# Eigen_INCLUDE_DIRS - Eigen include directories.


# Find the headers.
find_path(Eigen_INCLUDE_DIR Eigen/Core
        HINTS ${Eigen_DIR})

# Find Eigen version.
# Sometime would be good to add in a check to make sure we are finding Eigen 3.

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
        REQUIRED_VARS Eigen_INCLUDE_DIR)

mark_as_advanced(Eigen_INCLUDE_DIR)

set(Eigen_INCLUDE_DIRS ${Eigen_INCLUDE_DIR})
