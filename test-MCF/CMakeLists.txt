add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)

add_executable(MCFSolver_test test.cpp)

target_link_libraries(MCFSolver_test
        SMS++
        MCFBlock
        MCFClass
        MILPSolver
        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libcplex.a
        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libilocplex.a
        netcdf
        netcdf-cxx4)
