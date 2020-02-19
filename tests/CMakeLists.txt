# Create the executable unit test target that uses the Catch2 testing framework.
add_executable(mast_catch_tests
    test_helpers.cpp
    init_catch2.cpp)

# Add include directories to the unit test executable target.
target_include_directories(mast_catch_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/contrib/catch
        ${CMAKE_CURRENT_LIST_DIR})

# Link the unit test executable against the MAST library.
target_link_libraries(mast_catch_tests mast)

# TODO: May be better to use Catch2's built in CMake support rather than manually adding through ctest

add_subdirectory(base)
add_subdirectory(boundary_condition)
add_subdirectory(material)
add_subdirectory(property)
add_subdirectory(element)

message(NOTICE "It is recommended to run 'make check' instead of 'make test'. Alternatively, for 'ctest' or \
'make test' to output Catch2 error messages when a failure occurs, you must set the environment variable \
CTEST_OUTPUT_ON_FAILURE=1.")

# Add a customized "test" command.  We can't override "test", so we call it 
# "check" instead. This "check" command will generate a graph to visualize test
# dependencies, run ctest with output-on-failure enabled so that Catch2 output
# is available, and then generate a graph to visualize which tests passed, 
# failed, or were not ran.
#
# References 
# ----------
# 'CMake: setting an environment variable for ctest (or otherwise getting failed test output from ctest/make test automtically)', stackoverflow.com
# 'Replace ctest command with "ctest --output-on-failure" permanently for a specific project in CMakeLists.txt', stackoverflow.com
# 'CMake ignore return on custom target', stackoverflow.com
#
if (CMAKE_CONFIGURATION_TYPES)
    add_custom_target(check
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/GenerateTestDependencyGraph.py ${CMAKE_CURRENT_LIST_DIR}/CMakeLists.txt
        COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure --output-log ./test_results.log --build-config "$<CONFIGURATION>" || (exit 0)
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/GenerateTestResultsGraph.py ./test_dependency_graph.dot ./test_results.log)
else()
    add_custom_target(check
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/GenerateTestDependencyGraph.py ${CMAKE_CURRENT_LIST_DIR}/CMakeLists.txt
        COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure --output-log ./test_results.log || (exit 0)
        COMMAND python ${CMAKE_CURRENT_LIST_DIR}/GenerateTestResultsGraph.py ./test_dependency_graph.dot ./test_results.log)
endif()


# A custom target called "test_graph" that generates a graph of test 
# dependencies when run. This allows test dependencies to be visualized without
# having to run the entire "make check" command.
add_custom_target(test_graph
    COMMAND python ${CMAKE_CURRENT_LIST_DIR}/GenerateTestDependencyGraph.py ${CMAKE_CURRENT_LIST_DIR}/CMakeLists.txt)
