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

    # This version replicates the batch file in CMake.
    # This method was preferred because it allows to keep track
    # of the single tests via the CTest application and also
    # avoids copying all the auxiliary files.

    foreach (size 10 50)
          foreach (seed RANGE 0 19)
                add_test(NAME "test_cuts/[${seed},${size},3]"
                         COMMAND test_cuts ${seed} ${size} 3
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
          endforeach ()
    endforeach ()
endif ()

# --------------------------------------------------------------------------- #
