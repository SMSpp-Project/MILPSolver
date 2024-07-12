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
    target_link_libraries(cpx_pars PRIVATE CPLEX::Cplex)

    # Custom command to run cpx_pars in the tools directory
    add_custom_command(
            TARGET cpx_pars
            POST_BUILD
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cpx_pars
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running cpx_pars to generate necessary headers"
            VERBATIM
    )

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS cpx_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_SCIP)
    # ----- SCIP Solver tool ------------------------------------------------ #
    add_executable(scip_pars scip_pars.cpp)
    target_compile_features(scip_pars PRIVATE cxx_std_17)
    target_link_libraries(scip_pars PRIVATE SCIP::SCIP)

    # Custom command to run scip_pars in the tools directory
    add_custom_command(
            TARGET scip_pars
            POST_BUILD
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/scip_pars
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running scip_pars to generate necessary headers"
            VERBATIM
    )

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS scip_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_GUROBI)
    # ----- GRB Solver tool ------------------------------------------------- #
    add_executable(grb_pars grb_pars.cpp)
    target_compile_features(grb_pars PRIVATE cxx_std_17)
    target_link_libraries(grb_pars PRIVATE GUROBI::Gurobi)

    # Custom command to run grb_pars in the tools directory
    add_custom_command(
            TARGET grb_pars
            POST_BUILD
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/grb_pars
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running grb_pars to generate necessary headers"
            VERBATIM
    )

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS grb_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (MILPSolver_USE_HiGHS)
    # ----- HIGHS Solver tool ----------------------------------------------- #
    add_executable(highs_pars highs_pars.cpp)
    target_compile_features(highs_pars PRIVATE cxx_std_17)
    target_link_libraries(highs_pars PRIVATE HiGHS::HiGHS)

    # Custom command to run highs_pars in the tools directory
    add_custom_command(
            TARGET highs_pars
            POST_BUILD
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/highs_pars
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running highs_pars to generate necessary headers"
            VERBATIM
    )

    # ----- Install instructions -------------------------------------------- #
    include(GNUInstallDirs)
    install(TARGETS highs_pars
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

# --------------------------------------------------------------------------- #