# ----- MILPSolver_test --------------------------------------------------------
find_package(SMILPBLock REQUIRED)
add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::SMILPBLock
                      SMS++::MILPSolver)

# ----- MCFSolver_test ---------------------------------------------------------
# TODO: Do it in a more fashionable way and avoid duplication
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

find_package(MCFBlock REQUIRED)
find_package(netCDF REQUIRED)
add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test PRIVATE
                      SMS++::MCFBlock
                      SMS++::MILPSolver
                      NetCDF::NetCDFCxx)
