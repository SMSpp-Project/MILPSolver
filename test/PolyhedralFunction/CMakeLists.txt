# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test/PolyhedralFunction                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Enrico Calandrini                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

if ( TARGET SMS++::MILPSolver )
    add_executable(PolyhedralFunctionLP_test test.cpp)

    target_link_libraries(PolyhedralFunctionLP_test PRIVATE
                            SMS++::MILPSolver)

    set_target_properties(PolyhedralFunctionLP_test PROPERTIES COMPILE_FLAGS "-O3 -std=c++17 -DNDEBUG")

    # This version replicates the batch file in CMake.
    # This method was preferred because it allows to keep track
    # of the single tests via the CTest application and also
    # avoids copying all the auxiliary files.

    foreach(what 31 63 127)
        foreach (size 2 10 50)
            foreach (dens 3 1.1)
                foreach (seed RANGE 0 19)
                    add_test(NAME "PolyhedralFunctionLP_test/[${seed},${what},${size},${dens}]"
                            COMMAND PolyhedralFunctionLP_test ${seed} ${what} ${size} ${dens}
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
                endforeach ()
            endforeach ()
        endforeach ()
    endforeach()
endif ()

# --------------------------------------------------------------------------- #
