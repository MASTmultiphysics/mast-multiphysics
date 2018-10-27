# This module relies on SLEPc_DIR being set.
#
# SLEPc_FOUND - system has SLEPc.
# SLEPc_INCLUDE_DIRS - SLEPc include directories.
# SLEPc_LIBRARIES - SLEPc libraries.

# Find the headers.
find_path(SLEPc_INCLUDE_DIR slepc.h
          HINTS ${SLEPc_DIR}/include)

# Find the libraries.
find_library(SLEPc_LIBRARY
             NAMES slepc
             HINTS ${SLEPc_DIR}/lib)

# Find SLEPc version.
if(SLEPc_INCLUDE_DIR)
    set(HEADER "${SLEPc_INCLUDE_DIR}/slepcversion.h")
    file(STRINGS "${HEADER}" major REGEX "define +SLEPC_VERSION_MAJOR")
    file(STRINGS "${HEADER}" minor REGEX "define +SLEPC_VERSION_MINOR")
    file(STRINGS "${HEADER}" patch REGEX "define +SLEPC_VERSION_SUBMINOR")
    string(REGEX REPLACE ".+([0-9]+)" "\\1" major ${major})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" minor ${minor})
    string(REGEX REPLACE ".+([0-9]+)" "\\1" patch ${patch})
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

mark_as_advanced(SLEPc_INCLUDE_DIR SLEPc_LIBRARY SLEPc_VERSION)

set(SLEPc_LIBRARIES ${SLEPc_LIBRARY})
set(SLEPc_INCLUDE_DIRS ${SLEPc_INCLUDE_DIR})
