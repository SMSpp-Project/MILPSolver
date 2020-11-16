# --------------------------------------------------------------------------- #
#    CMake file for tools                                                     #
# --------------------------------------------------------------------------- #

if (MILPSolver_USE_CPLEX)
    add_executable(cpx_pars cpx_pars.cpp)
    target_compile_features(cpx_pars PUBLIC cxx_std_17)
    target_link_libraries(cpx_pars PUBLIC CPLEX::Cplex)
endif ()
