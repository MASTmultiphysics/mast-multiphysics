
# Define the target
add_executable(fluid_jacobian  check_fluid_jacobian.cpp)

target_include_directories(fluid_jacobian
                           PRIVATE
                           ${MAST_TEST_DIR})

target_link_libraries(fluid_jacobian
                      mast
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME fluid_jacobian COMMAND fluid_jacobian)

