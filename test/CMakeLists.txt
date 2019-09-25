if (NOT SMS++_CONFIG_ALL)
    # Trying to find them while configuring all at once
    # gave some errors (probably race conditions).
    find_package(SMILPBLock REQUIRED)
    find_package(MCFBlock REQUIRED)
endif ()

# ----- MILPSolver_test -------------------------------------------------------
add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::SMILPBLock)

# ----- MCFSolver_test ---------------------------------------------------------
find_package(netCDF REQUIRED)
add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::MCFBlock
                      NetCDF::NetCDFCxx)
