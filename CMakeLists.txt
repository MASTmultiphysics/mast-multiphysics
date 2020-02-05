# PREAMBLE
cmake_minimum_required(VERSION 3.13)
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
option(ENABLE_GCMMA     "Build with GCMMA interface"  OFF)
option(ENABLE_DOT       "Build with DOT interface"    OFF)
option(ENABLE_SNOPT     "Build with SNOPT interface"  OFF)
option(ENABLE_NLOPT     "Build with NLOPT interface"  OFF)
option(ENABLE_NASTRANIO "Build with support for reading Nastran meshes" OFF)
option(ENABLE_CYTHON    "Build with support for Cython development"     OFF)
option(BUILD_DOC        "Build documentation"         OFF)

# Dependency paths (optional) - these paths help CMake find dependencies.
set(MAST_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(libMesh_DIR "libMesh_DIR" CACHE PATH "Directory containing libMesh include/ and lib/")
set(PETSc_DIR   "PETSc_DIR"   CACHE PATH "Directory containing PETSc include/ and lib/")
set(PETSc_ARCH  "PETSc_ARCH"  CACHE STRING "Configuration/build of PETSc that should be used.")
set(SLEPc_DIR   "SLEPc_DIR"   CACHE PATH "Directory containing SLEPc include/ and lib/")
set(SLEPc_ARCH  "SLEPc_ARCH"  CACHE STRING "Configuration/build of SLEPc that should be used.")
set(EIGEN3_ROOT "Eigen_DIR"   CACHE PATH "Directory containing eigen header files")
set(DOT_DIR     "DOT_DIR"     CACHE PATH "Directory containing DOT lib/")
set(SNOPT_DIR   "SNOpt_DIR"   CACHE PATH "Directory containing SNOPT lib/")
set(NLOPT_DIR   "NLOpt_DIR"   CACHE PATH "Directory containing NLOpt include/ and lib/")

# EXTERNALLY PROVIDED CONTENT
# None. Use this if we pull something in during the build in the future.

# FIND DEPENDENCIES
find_package(MPI REQUIRED)
find_package(LAPACK REQUIRED)
message(STATUS "Found BLAS libs: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK libs: ${LAPACK_LIBRARIES}")
find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
find_package(HDF5 REQUIRED)
find_package(libMesh REQUIRED)
find_package(Eigen3 REQUIRED)
# Note: for some reason, letting CMake find a package configuration file for Boost from
#       a CMake/Spack-based installation of Boost causes trouble. Some of the include 
#       directories seem to be missing in this case. We disable the search for a 
#       Boost CMake configuration in the FindBoost module for now.
set(Boost_NO_BOOST_CMAKE ON) 
find_package(Boost COMPONENTS iostreams system filesystem unit_test_framework REQUIRED)
message(STATUS "Found Boost include: ${Boost_INCLUDE_DIRS}") # To ensure the found Boost is what user expected
message(STATUS "Found Boost libs:    ${Boost_LIBRARY_DIRS}") # To ensure the found Boost is what user expected

# Find optional packages.
if (ENABLE_GCMMA)
    find_package(GCMMA REQUIRED)
    find_package(GFortranLibs REQUIRED)
else()
    set (MAST_ENABLE_GCMMA 0)
endif()

if (ENABLE_DOT)
    find_package(DOT REQUIRED)
    find_package(GFortranLibs REQUIRED)
else()
    set (MAST_ENABLE_DOT 0)
endif()

if (ENABLE_SNOPT)
    find_package(SNOPT REQUIRED)
    find_package(GFortranLibs REQUIRED)
else()
    set (MAST_ENABLE_SNOPT 0)
endif()

if (ENABLE_NLOPT)
   find_package(NLOpt REQUIRED)
else()
    set (MAST_ENABLE_NLOPT 0)
endif()

if (ENABLE_NASTRANIO)
    find_package(Python3 REQUIRED)

    # Make sure Python has pyNastran.
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import pyNastran"
        RESULT_VARIABLE PYNASTRAN_ERROR)
    if(PYNASTRAN_ERROR)
        message(FATAL_ERROR "Error finding pyNastran package for Python3_EXECUTABLE!")
    else()
        message(STATUS "  pyNastran package found")
    endif()

    # Make sure Python has NumPy.
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE NUMPY_ERROR)
    if(NUMPY_ERROR)
        message(FATAL_ERROR "Error finding pyNastran package for Python3_EXECUTABLE!")
    else()
        message(STATUS "  Numpy package found")
    endif()
    
    set (MAST_ENABLE_NASTRANIO 1)
else()
    set (MAST_ENABLE_NASTRANIO 0)
endif()

# THIRD PARTY/CONTRIB
# - This directory contains files developed by third parties and are included in the
#   MAST source for convenience.
add_subdirectory(contrib)

# MAIN TARGETS
add_subdirectory(src)

# DOCUMENTATION
if(BUILD_DOC)
    add_subdirectory(doc)
endif()

# TESTS
enable_testing()
add_subdirectory(tests)

# EXAMPLES
add_subdirectory(examples)

