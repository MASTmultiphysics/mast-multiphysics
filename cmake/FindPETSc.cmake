# This module relies on PETSc_DIR being set.
#
# PETSc_FOUND - system has PETSc.
# PETSc_INCLUDE_DIRS - PETSc include directories.
# PETSc_LIBRARIES - PETSc libraries.

# Find the headers.
find_path(PETSc_INCLUDE_DIR petsc.h
          HINTS ${PETSc_DIR}/include)

# Find the libraries.
find_library(PETSc_LIBRARY
             NAMES petsc
             HINTS ${PETSc_DIR}/lib)

# Find PETSc version.
if(PETSc_INCLUDE_DIR)
    set(HEADER "${PETSc_INCLUDE_DIR}/petscversion.h")
    file(STRINGS "${HEADER}" major REGEX "define +PETSC_VERSION_MAJOR")
    file(STRINGS "${HEADER}" minor REGEX "define +PETSC_VERSION_MINOR")
    file(STRINGS "${HEADER}" patch REGEX "define +PETSC_VERSION_SUBMINOR")
    string(REGEX REPLACE ".+([0-9]+)" "\\1" major ${major})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" minor ${minor})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" patch ${patch})
    string(STRIP "${major}" major)
    string(STRIP "${minor}" minor)
    string(STRIP "${patch}" patch)
    set(PETSc_VERSION "${major}.${minor}.${patch}")
endif()

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
        REQUIRED_VARS PETSc_LIBRARY PETSc_INCLUDE_DIR
        VERSION_VAR PETSc_VERSION)

mark_as_advanced(PETSc_INCLUDE_DIR PETSc_LIBRARY PETSc_VERSION PETSc_FOUND)

set(PETSc_LIBRARIES ${PETSc_LIBRARY})
set(PETSc_INCLUDE_DIRS ${PETSc_INCLUDE_DIR})
