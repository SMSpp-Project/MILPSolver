# --------------------------------------------------------------------------- #
#    CMake file for tools                                                     #
# --------------------------------------------------------------------------- #

find_package(CPLEX REQUIRED)

add_executable(cpx_pars cpx_pars.cpp)
target_compile_features(cpx_pars PUBLIC cxx_std_17)
target_link_libraries(cpx_pars PUBLIC CPLEX::Cplex)
