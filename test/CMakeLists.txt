add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test
        SMS++
        SMILPBLock
        MILPSolver
        cplex
        ilocplex
        netcdf
        netcdf-cxx4)

add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test
        SMS++
        MCFBlock
        MCFClass
        MILPSolver
        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libcplex.a
        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libilocplex.a
        netcdf
        netcdf-cxx4)
