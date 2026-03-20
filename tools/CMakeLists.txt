# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/tools                                          #
# --------------------------------------------------------------------------- #

# Function to remove dots from version strings
function(remove_dots VERSION OUTPUT_VAR REPLACE_WITH_ZERO)
    if (REPLACE_WITH_ZERO)
        string(REPLACE "." "0" _VERSION_MODIFIED ${VERSION})
        set(_VERSION_MODIFIED "${_VERSION_MODIFIED}00")
    else()
        string(REPLACE "." "" _VERSION_MODIFIED ${VERSION})
    endif()
    set(${OUTPUT_VAR} ${_VERSION_MODIFIED} PARENT_SCOPE)
endfunction()

set(MILPSolver_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_custom_target(generate_headers)

# ----- Parameter generators ------------------------------------------------ #

if (MILPSolver_USE_CPLEX)
    # Remove dots from CPLEX version
    remove_dots(${CPLEX_VERSION} CPLEX_VERSION_NO_DOTS ON)

    set(CPLEX_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/CPX${CPLEX_VERSION_NO_DOTS}_defs.h)
    set(CPLEX_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/CPX${CPLEX_VERSION_NO_DOTS}_maps.h)

    add_executable(cpx_pars cpx_pars.cpp)
    target_compile_features(cpx_pars PRIVATE cxx_std_17)
    target_link_libraries(cpx_pars PRIVATE SMS++::SMS++ CPLEX::Cplex)

    add_custom_command(
            OUTPUT ${CPLEX_DEFS_FILE} ${CPLEX_MAPS_FILE}
            COMMAND "$<TARGET_FILE:cpx_pars>"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS cpx_pars
            COMMENT "Running cpx_pars to generate necessary headers"
            VERBATIM
    )

    add_custom_target(generate_cplex_headers
            DEPENDS ${CPLEX_DEFS_FILE} ${CPLEX_MAPS_FILE}
    )

    add_dependencies(generate_headers generate_cplex_headers)
endif ()

if (MILPSolver_USE_SCIP)
    # Remove dots from SCIP version
    remove_dots(${SCIP_VERSION} SCIP_VERSION_NO_DOTS OFF)

    set(SCIP_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/SCIP${SCIP_VERSION_NO_DOTS}_defs.h)
    set(SCIP_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/SCIP${SCIP_VERSION_NO_DOTS}_maps.h)

    add_executable(scip_pars scip_pars.cpp)
    target_compile_features(scip_pars PRIVATE cxx_std_17)
    target_link_libraries(scip_pars PRIVATE SMS++::SMS++ SCIP::SCIP)

    add_custom_command(
            OUTPUT ${SCIP_DEFS_FILE} ${SCIP_MAPS_FILE}
            COMMAND ${CMAKE_COMMAND} -E env
            DYLD_LIBRARY_PATH=$<TARGET_FILE_DIR:SCIP::SCIP>:$<TARGET_FILE_DIR:SMS++::SMS++>:$ENV{DYLD_LIBRARY_PATH}
            "$<TARGET_FILE:scip_pars>"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS scip_pars
            COMMENT "Running scip_pars to generate necessary headers"
            VERBATIM
    )

    add_custom_target(generate_scip_headers
            DEPENDS ${SCIP_DEFS_FILE} ${SCIP_MAPS_FILE}
    )

    add_dependencies(generate_headers generate_scip_headers)
endif ()

if (MILPSolver_USE_GUROBI)
    # Remove dots from GUROBI version
    remove_dots(${GUROBI_VERSION} GUROBI_VERSION_NO_DOTS OFF)

    set(GUROBI_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/GRB${GUROBI_VERSION_NO_DOTS}_defs.h)
    set(GUROBI_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/GRB${GUROBI_VERSION_NO_DOTS}_maps.h)

    add_executable(grb_pars grb_pars.cpp)
    target_compile_features(grb_pars PRIVATE cxx_std_17)
    target_link_libraries(grb_pars PRIVATE SMS++::SMS++ GUROBI::Gurobi)

    add_custom_command(
            OUTPUT ${GUROBI_DEFS_FILE} ${GUROBI_MAPS_FILE}
            COMMAND $<TARGET_FILE:grb_pars>
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS grb_pars
            COMMENT "Running grb_pars to generate necessary headers"
            VERBATIM
    )

    add_custom_target(generate_gurobi_headers
            DEPENDS ${GUROBI_DEFS_FILE} ${GUROBI_MAPS_FILE}
    )

    add_dependencies(generate_headers generate_gurobi_headers)
endif ()

if (MILPSolver_USE_HiGHS)
    # Remove dots from HiGHS version
    remove_dots(${HiGHS_VERSION} HiGHS_VERSION_NO_DOTS OFF)

    set(HiGHS_DEFS_FILE ${MILPSolver_INCLUDE_DIR}/HiGHS${HiGHS_VERSION_NO_DOTS}_defs.h)
    set(HiGHS_MAPS_FILE ${MILPSolver_INCLUDE_DIR}/HiGHS${HiGHS_VERSION_NO_DOTS}_maps.h)

    if (CMAKE_CROSSCOMPILING)
        add_custom_command(
                OUTPUT ${HiGHS_DEFS_FILE} ${HiGHS_MAPS_FILE}
                COMMAND /bin/sh -c
                        "'${CMAKE_MAKE_PROGRAM}' -C '${CMAKE_CURRENT_SOURCE_DIR}' \
                        HiGHS_ROOT='$$BUILD_PREFIX' highs_pars"
                COMMAND /bin/sh -c "'${CMAKE_CURRENT_SOURCE_DIR}/highs_pars'"
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                COMMENT "Running native highs_pars to generate necessary headers"
                VERBATIM
        )
    else ()
        add_executable(highs_pars highs_pars.cpp)
        target_compile_features(highs_pars PRIVATE cxx_std_17)
        target_link_libraries(highs_pars PRIVATE SMS++::SMS++ HiGHS::HiGHS)

        add_custom_command(
                OUTPUT ${HiGHS_DEFS_FILE} ${HiGHS_MAPS_FILE}
                COMMAND $<TARGET_FILE:highs_pars>
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DEPENDS highs_pars
                COMMENT "Running highs_pars to generate necessary headers"
                VERBATIM
        )
    endif ()

    add_custom_target(generate_highs_headers
            DEPENDS ${HiGHS_DEFS_FILE} ${HiGHS_MAPS_FILE}
    )

    add_dependencies(generate_headers generate_highs_headers)

endif ()

# --------------------------------------------------------------------------- #
