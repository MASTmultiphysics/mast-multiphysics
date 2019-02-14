# This file defines the build target for generating documentation with
# a combination of Doxygen/Python. It assumes that CMake will be able
# to find Doxygen and Python 3 on your system somewhere using
# CMake's (>= version 3.2) default modules for:
#    - Doxygen - https://cmake.org/cmake/help/v3.2/module/FindDoxygen.html
#    - Python  - https://cmake.org/cmake/help/v3.2/module/FindPythonInterp.html

# Doxygen is used to create API documentation for classes in src/. Python is to parse
# example_XX files in the examples/ directory and create HTML output for
# example problems/tutorials.

# Look for Doxygen executable on user's system. If successful, we get the
# following CMake cache variables.
#    - DOXYGEN_EXECUTABLE     = The path to the doxygen command.
#    - DOXYGEN_FOUND          = Was Doxygen found or not?
#    - DOXYGEN_VERSION        = The version reported by doxygen --version
#    - DOXYGEN_DOT_EXECUTABLE = The path to the dot program used by doxygen.
#    - DOXYGEN_DOT_FOUND      = Was Dot found or not?
find_package(Doxygen)

# Look for Python (at least version XX) on the user's system. If successful,
# we get the following CMake cache variables.
#    - PYTHONINTERP_FOUND    = Was the Python executable found
#    - PYTHON_EXECUTABLE     = path to the Python interpreter
#    - PYTHON_VERSION_STRING = Python version found e.g. 2.5.2
#    - PYTHON_VERSION_MAJOR  = Python major version found e.g. 2
#    - PYTHON_VERSION_MINOR  = Python minor version found e.g. 5
#    - PYTHON_VERSION_PATCH  = Python patch version found e.g. 2
find_package(PythonInterp 3 REQUIRED)

if(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doc_doxygen
            COMMAND ${PYTHON_EXECUTABLE} ${MAST_ROOT_DIR}/doc/tutorial_doc.py ${MAST_ROOT_DIR}/examples ${MAST_ROOT_DIR}/doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Python/Doxygen."
            VERBATIM)

else()
    message(WARNING "Doxygen and Python 3 are required to generate the documentation!")
endif()
