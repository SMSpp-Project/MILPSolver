add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test
        SMS++
        SMILPBLock
        MILPSolver
        cplex-library
        netcdf
        netcdf-cxx4)

add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test
        SMS++
        MCFBlock
        MCFClass
        MILPSolver
        cplex-library
        netcdf
        netcdf-cxx4)
