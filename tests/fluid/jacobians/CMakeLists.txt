
# Define the target
add_executable(fluid_jacobian   check_fluid_jacobian.cpp)
add_executable(fluid_eigen_vec  check_fluid_eigenvectors.cpp)

target_include_directories(fluid_jacobian
                           PRIVATE
                           ${MAST_TEST_DIR})

target_include_directories(fluid_eigen_vec
                            PRIVATE
                            ${MAST_TEST_DIR})

target_link_libraries(fluid_jacobian
                      mast
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

target_link_libraries(fluid_eigen_vec
                        mast
                        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME fluid_jacobian COMMAND fluid_jacobian)
add_test(NAME fluid_eigen_vec COMMAND fluid_eigen_vec)



