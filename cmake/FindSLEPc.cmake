# This module relies on SLEPc_DIR being set.
#
# SLEPc_FOUND - system has SLEPc.
# SLEPc_INCLUDE_DIRS - SLEPc include directories.
# SLEPc_LIBRARIES - SLEPc libraries.

# Find the headers.
# Search CMake variable paths first (SLEPc_DIR) and then environment variable paths next (SLEPC_DIR)
find_path(SLEPc_INCLUDE_DIR1 slepc.h
          HINTS "${SLEPc_DIR}/include"
                "${SLEPc_DIR}/${SLEPc_ARCH}/include" 
                "$ENV{SLEPC_DIR}/include"
                "$ENV{SLEPC_DIR}/$ENV{SLEPC_ARCH}/include")

# Search for slepcconf.h, which could be in a different location than slepc.h if an ARCH was specified
find_path(SLEPc_INCLUDE_DIR2 slepcconf.h
          HINTS "${SLEPc_DIR}/include"
                "${SLEPc_DIR}/${SLEPc_ARCH}/include" 
                "$ENV{SLEPC_DIR}/include"
                "$ENV{SLEPC_DIR}/$ENV{SLEPC_ARCH}/include")

# Use both SLEPc include directories found above
set(SLEPc_INCLUDE_DIR "${SLEPc_INCLUDE_DIR1};${SLEPc_INCLUDE_DIR2}")

# Find the libraries. 
# Search CMake variable paths first (SLEPc_DIR) and then environment variable paths next (SLEPC_DIR)
find_library(SLEPc_LIBRARY
             NAMES slepc
             HINTS "${SLEPc_DIR}/lib"
                   "${SLEPc_DIR}/${SLEPc_ARCH}/lib"
                   "$ENV{SLEPC_DIR}/lib"
                   "$ENV{SLEPC_DIR}/$ENV{SLEPC_ARCH}/lib")

# Find SLEPc version.
if(SLEPc_INCLUDE_DIR)
    set(HEADER "${SLEPc_INCLUDE_DIR1}/slepcversion.h")
    file(STRINGS "${HEADER}" major REGEX "define +SLEPC_VERSION_MAJOR")
    file(STRINGS "${HEADER}" minor REGEX "define +SLEPC_VERSION_MINOR")
    file(STRINGS "${HEADER}" patch REGEX "define +SLEPC_VERSION_SUBMINOR")
    string(REGEX MATCH "#define SLEPC_VERSION_MAJOR *([0-9]*)" _ ${major})
    set(major ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define SLEPC_VERSION_MINOR *([0-9]*)" _ ${minor})
    set(minor ${CMAKE_MATCH_1})
    string(REGEX MATCH "#define SLEPC_VERSION_SUBMINOR *([0-9]*)" _ ${patch})
    set(patch ${CMAKE_MATCH_1})
    string(STRIP "${major}" major)
    string(STRIP "${minor}" minor)
    string(STRIP "${patch}" patch)
    set(SLEPc_VERSION "${major}.${minor}.${patch}")
endif()

# Set variables.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc
        REQUIRED_VARS SLEPc_LIBRARY SLEPc_INCLUDE_DIR
        VERSION_VAR SLEPc_VERSION)

mark_as_advanced(SLEPc_INCLUDE_DIR SLEPc_LIBRARY SLEPc_VERSION SLEPc_FOUND)

set(SLEPc_LIBRARIES ${SLEPc_LIBRARY})
set(SLEPc_INCLUDE_DIRS ${SLEPc_INCLUDE_DIR})
