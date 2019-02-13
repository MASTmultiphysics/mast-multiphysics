# This file creates the build target for generating documentation with
# Doxygen. It assumes that CMake will be able to find Doxygen on your system
# somewhere.
#
# Currently we only generate API documentation for classes in src/, but we can
# eventually extend this to generate documented example problems/tutorials.

# Look for Doxygen executable on user's system. If successful, we get the
# following CMake cache variables.
#    - DOXYGEN_EXECUTABLE     = The path to the doxygen command.
#    - DOXYGEN_FOUND          = Was Doxygen found or not?
#    - DOXYGEN_VERSION        = The version reported by doxygen --version
#    - DOXYGEN_DOT_EXECUTABLE = The path to the dot program used by doxygen.
#    - DOXYGEN_DOT_FOUND      = Was Dot found or not?
find_package(Doxygen)
find_package(Python3)

if(DOXYGEN_FOUND AND Python3_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doc_doxygen
            COMMAND ${Python3_EXECUTABLE} ${MAST_ROOT_DIR}/doc/tutorial_doc.py ${MAST_ROOT_DIR}/examples ${MAST_ROOT_DIR}/doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Doxygen."
            VERBATIM)

else()
    message(WARNING "Doxygen is required to generate the documentation!")
endif()
