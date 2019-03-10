
# Define the target
add_executable(fluid_jacobian)

target_sources(fluid_jacobian
                PRIVATE
                ${CMAKE_CURRENT_LIST_DIR}/check_fluid_jacobian.cpp)

target_include_directories(fluid_jacobian
                           PRIVATE
                           ${MAST_TEST_DIR})

target_link_libraries(fluid_jacobian
                      mast
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(fluid_jacobian_test fluid_jacobian)

