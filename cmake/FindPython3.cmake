# This module relies on Python3_DIR being set.
#
# The reason that we include this module rather than the CMake provided
# FindPythonInterp, FindPythonLibs, or FindPython is that they seem to be unable
# consistently find matching Python interpreter and libraries when there
# are multiple Python distributions available (especially on macOS).
# We use a HINT in this module to try to hone in on a single Python distribution.
#
# Python3_FOUND - Python 3 required components are on the system.
# Python3_EXECUTABLE - Path to Python 3 interpreter.
# Python3_LIBRARIES - Paths to Python 3 libraries.
# Python3_INCLUDE_DIRS - Python 3 include directories (including site-packages).

# Acceptable Python versions.
set(_PY_VERSIONS 3.8 3.7 3.6)

# Here we ensure MacOS system provided Python distributions are found by CMake
# last because they typically don't have the functionality we want.
set(CMAKE_FIND_FRAMEWORK LAST)
set(CMAKE_FIND_APPBUNDLE LAST)

# Build list of acceptable Python executables that we will look for.
set(_PY_NAMES "")
foreach(_PY_VER ${_PY_VERSIONS})
  list(APPEND _PY_NAMES "python${_PY_VER}")
endforeach()

# Find path to Python interpreter. We call find_program() twice to make sure we
# search user provided Python3_DIR (CMake variable) before trying default
# system environment paths.
find_program(Python3_EXECUTABLE
        NAMES ${_PY_NAMES}
        HINTS
          ${Python3_DIR}
          PATH_SUFFIXES
          "bin"
        NO_DEFAULT_PATH)
find_program(Python3_EXECUTABLE
        NAMES ${_PY_NAMES})

# Find root directory of Python distribution. Stored in Python_ROOT_DIR.
get_filename_component(Python3_BIN_DIR ${Python3_EXECUTABLE} DIRECTORY)
get_filename_component(Python3_ROOT_DIR ${Python3_BIN_DIR} DIRECTORY)
# Also, get name of Python interpreter that we actually found above.
# Stored in Python3_NAME.
get_filename_component(Python3_NAME ${Python3_EXECUTABLE} NAME)

# Find the Python headers.
# -- Currently we search relative to Python_ROOT_DIR we found above.
find_path(Python3_INCLUDE_DIR Python.h
        HINTS
          ${Python3_ROOT_DIR}/include/${Python3_NAME}/
          ${Python3_ROOT_DIR}/include/${Python3_NAME}m/)

# Find the Python library. Store also library directory to help
# find modules/site-packages.
# -- Currently we search relative to Python_ROOT_DIR we found above.
find_library(Python3_LIBRARY
        NAMES ${Python3_NAME} ${Python3_NAME}m
        HINTS ${Python3_ROOT_DIR}/lib)
get_filename_component(Python3_LIBRARY_DIR ${Python3_LIBRARY} DIRECTORY)

# Find the mpi4py headers required for Cython wrapper.
# find_path(mpi4py_INCLUDE_DIR MPI.pxd
#         HINTS
#         ${Python_LIBRARY_DIR}/${Python_NAME}/site-packages/mpi4py)

# find_path(Python_MODULEDIR mpi4py/MPI.pxd
#         HINTS
#         ${Python_LIBRARY_DIR}/${Python_NAME}/site-packages)

# Output what we found.
message(STATUS "  Python3_EXECUTABLE:  ${Python3_EXECUTABLE}")
message(STATUS "  Python3_INCLUDE_DIR: ${Python3_INCLUDE_DIR}")
message(STATUS "  Python3_LIBRARY:     ${Python3_LIBRARY}")
message(STATUS "  Python3 Modules:     ${Python3_MODULEDIR}")

# Set standard CMake variables and output status.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Python3
    REQUIRED_VARS
        Python3_EXECUTABLE
        Python3_LIBRARY
        Python3_INCLUDE_DIR)
  #      mpi4py_INCLUDE_DIR)
#message("--       mpi4py Headers: ${mpi4py_INCLUDE_DIR}")

# Advance cache variables.
mark_as_advanced(Python3_INCLUDE_DIR Python3_LIBRARY Python3_EXECUTABLE)

# Roll-up libraries and headers into standard CMake variables.
set(Python3_LIBRARIES    ${Python3_LIBRARY})
set(Python3_INCLUDE_DIRS ${Python3_INCLUDE_DIR})# ${mpi4py_INCLUDE_DIR})