# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/tools                                          #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- Parameter generators ------------------------------------------------ #

if (MILPSolver_USE_CPLEX)
    # ----- CPX Solver tool ------------------------------------------------- #
    add_executable(cpx_pars cpx_pars.cpp)
    target_compile_features(cpx_pars PRIVATE cxx_std_17)
    target_link_libraries(cpx_pars PRIVATE SMS++::MILPSolver)

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS cpx_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_SCIP)
    # ----- SCIP Solver tool ------------------------------------------------ #
    add_executable(scip_pars scip_pars.cpp)
    target_compile_features(scip_pars PRIVATE cxx_std_17)
    target_link_libraries(scip_pars PUBLIC SMS++::MILPSolver)
    target_include_directories(scip_pars PRIVATE ${SCIP_INCLUDE_DIRS})

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS scip_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_GUROBI)
    # ----- GRB Solver tool ------------------------------------------------- #
    add_executable(grb_pars grb_pars.cpp)
    target_compile_features(grb_pars PRIVATE cxx_std_17)
    target_link_libraries(grb_pars PRIVATE SMS++::MILPSolver)

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS grb_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_HiGHS)
    # ----- HIGHS Solver tool ----------------------------------------------- #
    add_executable(highs_pars highs_pars.cpp)
    target_compile_features(highs_pars PRIVATE cxx_std_17)
    target_link_libraries(highs_pars PRIVATE SMS++::MILPSolver)

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS highs_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

# ----- MILP Solver tool ---------------------------------------------------- #
add_executable(milp_solver solver.cpp)
target_link_libraries(milp_solver PRIVATE SMS++::MILPSolver)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)
install(TARGETS milp_solver
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --------------------------------------------------------------------------- #
