find_package(netCDF REQUIRED)

# ----- MILPSolver_test --------------------------------------------------------
#find_package(SMILPBLock REQUIRED)
add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_directories(MILPSolver_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::SMSpp
                      SMS++::SMILPBLock
                      SMS++::MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )

# ----- MCFSolver_test ---------------------------------------------------------
# TODO: Do it in a more fashionable way and avoid duplication
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

find_package(MCFClass REQUIRED)
find_package(MCFBlock REQUIRED)
add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_directories(MCFSolver_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MCFSolver_test
                      SMS++::SMSpp
                      SMS++::MCFBlock
                      SMS++::MILPSolver
                      MCFClass::MCFClass
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )
