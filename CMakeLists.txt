# PREAMBLE
cmake_minimum_required(VERSION 3.2)
project(MAST
        VERSION 0.1
        LANGUAGES C CXX)

# PROJECT WIDE SETUP
# Get CMake modules and set language standards.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options with defaults.
option(ENABLE_GCMMA   "Build with GCMMA interface"  OFF)
option(ENABLE_DOT     "Build with DOT interface"    OFF)
option(ENABLE_NPSOL   "Build with NPSOL interface"  OFF)
option(ENABLE_CYTHON  "Build with CYTHON interface" OFF)
option(BUILD_DOC      "Build documentation"         OFF)

# Required dependency paths.
set(libMesh_DIR "libMesh_DIR" CACHE PATH "Directory containing libMesh include/ and lib/")
set(PETSc_DIR   "PETSc_DIR"   CACHE PATH "Directory containing PETSc include/ and lib/")
set(SLEPc_DIR   "SLEPc_DIR"   CACHE PATH "Directory containing SLEPc include/ and lib/")

# EXTERNALLY PROVIDED CONTENT
# None. Use this if we pull something in during the build in the future.

# FIND DEPENDENCIES
find_package(MPI REQUIRED)
find_package(PETSc)
find_package(SLEPc)
find_package(libMesh)

# MAIN TARGETS
add_subdirectory(src)

# EXAMPLES
add_subdirectory(ex)

# DOCUMENTATION
if(BUILD_DOC)
    add_subdirectory(doc)
endif()