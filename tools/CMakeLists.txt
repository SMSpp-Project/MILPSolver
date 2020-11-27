# --------------------------------------------------------------------------- #
#    CMake file for tools                                                     #
# --------------------------------------------------------------------------- #

if (MILPSolver_USE_CPLEX)
    add_executable(cpx_pars cpx_pars.cpp)
    target_compile_features(cpx_pars PRIVATE cxx_std_17)
    target_link_libraries(cpx_pars PRIVATE CPLEX::Cplex)
endif ()
