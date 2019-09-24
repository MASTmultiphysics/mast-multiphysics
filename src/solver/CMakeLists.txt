target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/arclength_continuation_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/arclength_continuation_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/complex_solver_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/complex_solver_base.h
        ${CMAKE_CURRENT_LIST_DIR}/continuation_solver_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/continuation_solver_base.h
        ${CMAKE_CURRENT_LIST_DIR}/first_order_newmark_transient_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/first_order_newmark_transient_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/generalized_alpha_transient_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/generalized_alpha_transient_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/multiphysics_nonlinear_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/multiphysics_nonlinear_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/pseudo_arclength_continuation_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/pseudo_arclength_continuation_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/second_order_newmark_transient_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/second_order_newmark_transient_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/slepc_eigen_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/slepc_eigen_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/stabilized_first_order_transient_sensitivity_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/stabilized_first_order_transient_sensitivity_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/transient_solver_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/transient_solver_base.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/solver
        FILES_MATCHING PATTERN "*.h")
