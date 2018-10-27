target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/flutter_root_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/flutter_root_base.h
        ${CMAKE_CURRENT_LIST_DIR}/flutter_root_crossover_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/flutter_root_crossover_base.h
        ${CMAKE_CURRENT_LIST_DIR}/flutter_solution_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/flutter_solution_base.h
        ${CMAKE_CURRENT_LIST_DIR}/flutter_solver_base.cpp
        ${CMAKE_CURRENT_LIST_DIR}/flutter_solver_base.h
        ${CMAKE_CURRENT_LIST_DIR}/frequency_function.cpp
        ${CMAKE_CURRENT_LIST_DIR}/frequency_function.h
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_root.cpp
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_root.h
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_root_crossover.cpp
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_root_crossover.h
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_solution.cpp
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_solution.h
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/pk_flutter_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_root.cpp
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_root.h
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_root_crossover.cpp
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_root_crossover.h
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_solution.cpp
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_solution.h
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/time_domain_flutter_solver.h
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_root.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_root.h
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_root_crossover.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_root_crossover.h
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_solution.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_solution.h
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_solver.cpp
        ${CMAKE_CURRENT_LIST_DIR}/ug_flutter_solver.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/aeroelasticity
        FILES_MATCHING PATTERN "*.h")