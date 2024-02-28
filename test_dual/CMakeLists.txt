# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test_dual                                      #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Enrico Calandrini                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# ----- MILPSolver_test_dual ------------------------------------------------ #
add_executable(dual_test test_dual.cpp)
target_link_libraries(dual_test PRIVATE
                      SMS++::MILPSolver)

add_test(NAME dual_test
         COMMAND dual_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# --------------------------------------------------------------------------- #