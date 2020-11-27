# --------------------------------------------------------------------------- #
#    CMake file for tools                                                     #
# --------------------------------------------------------------------------- #

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
