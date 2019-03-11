# PREAMBLE
cmake_minimum_required(VERSION 3.2)
project(MAST
        VERSION 0.3
        LANGUAGES C CXX)

# PROJECT WIDE SETUP
# Get CMake modules and set language standards.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure installation sub-directories.
include(GNUInstallDirs)

# Build options with defaults.
option(ENABLE_GCMMA   "Build with GCMMA interface"  OFF)
option(ENABLE_DOT     "Build with DOT interface"    OFF)
option(ENABLE_NPSOL   "Build with NPSOL interface"  OFF)
option(ENABLE_NLOPT   "Build with NLOPT interface"  OFF)
option(ENABLE_CYTHON  "Build with CYTHON interface" OFF)
option(BUILD_DOC      "Build documentation"         OFF)

# Required dependency paths.
set(MAST_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(libMesh_DIR "libMesh_DIR" CACHE PATH "Directory containing libMesh include/ and lib/")
set(PETSc_DIR   "PETSc_DIR"   CACHE PATH "Directory containing PETSc include/ and lib/")
set(SLEPc_DIR   "SLEPc_DIR"   CACHE PATH "Directory containing SLEPc include/ and lib/")
set(EIGEN3_ROOT "Eigen_DIR"   CACHE PATH "Directory containing eigen header files")
set(NLOPT_DIR "NLOpt_DIR"   CACHE PATH "Directory containing NLOpt include/ and lib/")

# EXTERNALLY PROVIDED CONTENT
# None. Use this if we pull something in during the build in the future.

# FIND DEPENDENCIES
find_package(MPI REQUIRED)
find_package(LAPACK REQUIRED)
find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
find_package(libMesh REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# Find optional packages.
if (ENABLE_GCMMA)
    find_package(GCMMA REQUIRED)
endif()

if (ENABLE_DOT)
    find_package(DOT REQUIRED)
endif()

if (ENABLE_NLOPT)
   find_package(NLOpt REQUIRED)
endif()

# gfortran lib is required if GCMMA/DOT are used, since these are written
# in fortran
if (ENABLE_GCMMA OR ENABLE_DOT)
find_package(GFortranLibs REQUIRED)
endif()

# MAIN TARGETS
add_subdirectory(src)

# EXAMPLES
add_subdirectory(examples)

# DOCUMENTATION
if(BUILD_DOC)
    add_subdirectory(doc)
endif()

# TESTS
enable_testing()
add_subdirectory(tests)

