target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/fe_base.h
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.cpp
        ${CMAKE_CURRENT_LIST_DIR}/geom_elem.h
        ${CMAKE_CURRENT_LIST_DIR}/mesh_coupling_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/mesh_coupling_base.h)

# Only add nastran_io/Python data to `mast` target if we are building Nastran mesh support.
if(ENABLE_NASTRANIO)

    # Cython is currently not required unless you want to make changes in pynastran_io.pyx
    # and generate new pynastran_io.h/.cpp. If you want to modify pynastran_io.pyx, make sure
    # to set ENABLE_CYTHON=ON and CMake will ensure Cython is executed to generate updated
    # bindings prior to all `mast` target builds.
    if (ENABLE_CYTHON)

        # Make sure we have the Cython compiler.
        find_package(Cython REQUIRED)

        # Define the target to run Cython.
        add_custom_target(cython_pynastran
            DEPENDS
                ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.pyx
                COMMAND
                ${Cython_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.pyx --cplus -3 -Werror #-Wextra
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            BYPRODUCTS
                ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.h
                ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.cpp
            COMMENT
                "Running Cython to generate pyNastran bindings.")

        # Make the MAST library target `mast` dependency on running Cython.
        add_dependencies(mast cython_pynastran)
    endif()

    # Add sources to mast target.
    target_sources(mast
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.h
            ${CMAKE_CURRENT_LIST_DIR}/pynastran_io.cpp
            ${CMAKE_CURRENT_LIST_DIR}/nastran_io.cpp
            ${CMAKE_CURRENT_LIST_DIR}/nastran_io.h)

    # Find Numpy header path.
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
            RESULT_VARIABLE NUMPY_ERROR
            OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NUMPY_ERROR)
        message(FATAL_ERROR "Error finding Numpy header path for Python3_EXECUTABLE!")
    else()
        message(STATUS "  Numpy header path: ${NUMPY_INCLUDE_PATH}")
    endif()

    # Add header/library to mast target.
    target_include_directories(mast
        PUBLIC
            ${Python3_INCLUDE_DIRS}
            ${NUMPY_INCLUDE_PATH})
    target_link_libraries(mast PUBLIC ${Python3_LIBRARIES})
endif()

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/mesh
        FILES_MATCHING PATTERN "*.h")
