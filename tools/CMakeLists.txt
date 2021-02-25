# --------------------------------------------------------------------------- #
#    CMake file for tools                                                     #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else ()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- Parameter generators ------------------------------------------------ #
if (MILPSolver_USE_CPLEX)
    add_executable(cpx_pars cpx_pars.cpp)
    target_compile_features(cpx_pars PRIVATE cxx_std_17)
    target_link_libraries(cpx_pars PRIVATE CPLEX::Cplex)
endif ()

if (MILPSolver_USE_SCIP)
    add_executable(scip_pars scip_pars.cpp)
    target_compile_features(scip_pars PRIVATE cxx_std_17)
    target_link_libraries(scip_pars PUBLIC ${SCIP_LIBRARIES})
    target_include_directories(scip_pars PRIVATE ${SCIP_INCLUDE_DIRS})
endif ()

# ----- MILP Solver tool ---------------------------------------------------- #
add_executable(milp_solver solver.cpp ../test/SimpleMILPBlock.cpp)
target_link_libraries(milp_solver PRIVATE SMS++::MILPSolver)

# --------------------------------------------------------------------------- #
