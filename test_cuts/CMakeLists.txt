# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test_dynamic                                   #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Antonio Frangioni                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

if (TARGET SMS++::MILPSolver)
    add_executable(test_cuts test.cpp)

    target_link_libraries(test_cuts PRIVATE
                          SMS++::MILPSolver)
endif ()

# --------------------------------------------------------------------------- #
