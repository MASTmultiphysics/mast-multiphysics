# This module relies on PETSc_DIR and PETSc_ARCH being set in CMake or
# PETSC_DIR and PETSC_ARCH being set in the environment
#
# PETSc_FOUND - system has PETSc.
# PETSc_INCLUDE_DIRS - PETSc include directories.
# PETSc_LIBRARIES - PETSc libraries.

# Find the headers.
# Search CMake variable paths first (PETSc_DIR) and then environment variable paths next (PETSC_DIR)
find_path(PETSc_INCLUDE_DIR petsc.h
          HINTS "${PETSc_DIR}/include"
                "${PETSc_DIR}/${PETSc_ARCH}/include" 
                "$ENV{PETSC_DIR}/include"
                "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include")

# Find the libraries. 
# Search CMake variable paths first (PETSc_DIR) and then environment variable paths next (PETSC_DIR)
find_library(PETSc_LIBRARY
             NAMES petsc
             HINTS "${PETSc_DIR}/lib"
                   "${PETSc_DIR}/${PETSc_ARCH}/lib"
                   "$ENV{PETSC_DIR}/lib"
                   "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib")

# Find PETSc version.
if(PETSc_INCLUDE_DIR)
    set(HEADER "${PETSc_INCLUDE_DIR}/petscversion.h")
    file(STRINGS "${HEADER}" major REGEX "define +PETSC_VERSION_MAJOR")
    file(STRINGS "${HEADER}" minor REGEX "define +PETSC_VERSION_MINOR")
    file(STRINGS "${HEADER}" patch REGEX "define +PETSC_VERSION_SUBMINOR")
    string(REGEX MATCH "#define PETSC_VERSION_MAJOR *([0-9]*)" _ ${major})
    set(major ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define PETSC_VERSION_MINOR *([0-9]*)" _ ${minor})
    set(minor ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define PETSC_VERSION_SUBMINOR *([0-9]*)" _ ${patch})
    set(patch ${CMAKE_MATCH_1})
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
