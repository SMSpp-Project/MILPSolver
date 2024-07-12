# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/tools                                          #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Function to remove dots from version strings
function(remove_dots VERSION OUTPUT_VAR)
    string(REPLACE "." "" _VERSION_NO_DOTS ${VERSION})
    set(${OUTPUT_VAR} ${_VERSION_NO_DOTS} PARENT_SCOPE)
endfunction()

set(MILPSolver_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

# ----- Parameter generators ------------------------------------------------ #

if (MILPSolver_USE_CPLEX)
    # Remove dots from CPLEX version
    remove_dots(${CPLEX_VERSION} CPLEX_VERSION_NO_DOTS)

    set(CPLEX_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/CPX${CPLEX_VERSION_NO_DOTS}_defs.h)
    set(CPLEX_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/CPX${CPLEX_VERSION_NO_DOTS}_maps.h)

    # Only add custom command if files do not exist
    if (NOT (EXISTS ${CPLEX_DEFS_FILE} AND EXISTS ${CPLEX_MAPS_FILE}))
        set(GENERATE_CPLEX_HEADERS ON CACHE BOOL "Generate CPLEX headers" FORCE)

        # ----- CPX Solver tool --------------------------------------------- #
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

        # ----- Install instructions ---------------------------------------- #
        include(GNUInstallDirs)
        install(TARGETS cpx_pars
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        set(GENERATE_CPLEX_HEADERS OFF CACHE BOOL "Generate CPLEX headers" FORCE)
    endif ()
endif ()

if (MILPSolver_USE_SCIP)
    # Remove dots from SCIP version
    remove_dots(${SCIP_VERSION} SCIP_VERSION_NO_DOTS)

    set(SCIP_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/SCIP${SCIP_VERSION_NO_DOTS}_defs.h)
    set(SCIP_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/SCIP${SCIP_VERSION_NO_DOTS}_maps.h)

    # Only add custom command if files do not exist
    if (NOT (EXISTS ${SCIP_DEFS_FILE} AND EXISTS ${SCIP_MAPS_FILE}))
        set(GENERATE_SCIP_HEADERS ON CACHE BOOL "Generate SCIP headers" FORCE)

        # ----- SCIP Solver tool -------------------------------------------- #
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

        # ----- Install instructions ---------------------------------------- #
        include(GNUInstallDirs)
        install(TARGETS scip_pars
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        set(GENERATE_SCIP_HEADERS OFF CACHE BOOL "Generate SCIP headers" FORCE)
    endif ()
endif ()

if (MILPSolver_USE_GUROBI)
    # Remove dots from GUROBI version
    remove_dots(${GUROBI_VERSION} GUROBI_VERSION_NO_DOTS)

    set(GUROBI_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/GRB${GUROBI_VERSION_NO_DOTS}_defs.h)
    set(GUROBI_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/GRB${GUROBI_VERSION_NO_DOTS}_maps.h)

    # Only add custom command if files do not exist
    if (NOT (EXISTS ${GUROBI_DEFS_FILE} AND EXISTS ${GUROBI_MAPS_FILE}))
        set(GENERATE_GUROBI_HEADERS ON CACHE BOOL "Generate GUROBI headers" FORCE)

        # ----- GRB Solver tool --------------------------------------------- #
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

        # ----- Install instructions ---------------------------------------- #
        include(GNUInstallDirs)
        install(TARGETS grb_pars
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        set(GENERATE_GUROBI_HEADERS OFF CACHE BOOL "Generate GUROBI headers" FORCE)
    endif ()
endif ()

if (MILPSolver_USE_HiGHS)
    # Remove dots from HiGHS version
    remove_dots(${HiGHS_VERSION} HiGHS_VERSION_NO_DOTS)

    set(HiGHS_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/HiGHS${HiGHS_VERSION_NO_DOTS}_defs.h)
    set(HiGHS_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/HiGHS${HiGHS_VERSION_NO_DOTS}_maps.h)

    # Only add custom command if files do not exist
    if (NOT (EXISTS ${HiGHS_DEFS_FILE} AND EXISTS ${HiGHS_MAPS_FILE}))
        set(GENERATE_HIGHS_HEADERS ON CACHE BOOL "Generate HiGHS headers" FORCE)

        # ----- HIGHS Solver tool ------------------------------------------- #
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

        # ----- Install instructions ---------------------------------------- #
        include(GNUInstallDirs)
        install(TARGETS highs_pars
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        set(GENERATE_HIGHS_HEADERS OFF CACHE BOOL "Generate HiGHS headers" FORCE)
    endif ()
endif ()

# --------------------------------------------------------------------------- #
