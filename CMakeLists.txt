add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)

add_library(MILPSolver OBJECT src/MILPSolver.cpp src/CPXMILPSolver.cpp)
add_subdirectory(test)
