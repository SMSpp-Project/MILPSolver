add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)

add_executable(MILPSolver_test main.cpp)

target_link_libraries(MILPSolver_test
        SMS++
        SMILPBLock
        MILPSolver
        cplex
        ilocplex
        netcdf
        netcdf-cxx4)
