# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test_dynamic                                   #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Enrico Calandrini                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

if (TARGET SMS++::MILPSolver)
    add_executable(DynamicLP_test test.cpp)

    target_link_libraries(DynamicLP_test PRIVATE
                          SMS++::MILPSolver)
endif ()

# --------------------------------------------------------------------------- #
