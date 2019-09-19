find_package(netCDF REQUIRED)

# ----- Custom includes --------------------------------------------------------
# TODO: When libraries will be installable, fix this
set(SMS++_MAIN_DIR ${PROJECT_SOURCE_DIR}/../..)
set(MCFBlock_INCLUDE_DIR ${SMS++_MAIN_DIR}/MCFBlock)
set(SMILPBlock_INCLUDE_DIR ${SMS++_MAIN_DIR}/SMILPBlock)
set(MCFClass_INCLUDE_DIR ${SMS++_MAIN_DIR}/MCFClass)
set(MCFClass_INCLUDE_DIRS
    ${MCFClass_INCLUDE_DIR}/MCFClass
    ${MCFClass_INCLUDE_DIR}/MCFClone
    ${MCFClass_INCLUDE_DIR}/MCFCplex
    ${MCFClass_INCLUDE_DIR}/MCFSimplex
    ${MCFClass_INCLUDE_DIR}/OPTUtils
    ${MCFClass_INCLUDE_DIR}/RelaxIV
    )

# ----- MILPSolver_test --------------------------------------------------------
add_executable(MILPSolver_test milpsolver_test.cpp)
target_include_directories(MILPSolver_test PRIVATE ${MCFBlock_INCLUDE_DIR})
target_include_directories(MILPSolver_test PRIVATE ${SMILPBlock_INCLUDE_DIR})
target_include_directories(MILPSolver_test PRIVATE ${MCFClass_INCLUDE_DIRS})
target_include_directories(MILPSolver_test PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(MILPSolver_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MILPSolver_test
                      SMS++
                      SMILPBLock
                      MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )

# ----- MCFSolver_test ---------------------------------------------------------
# TODO: Do it in a more fashionable way and avoid duplication
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

add_executable(MCFSolver_test mcfsolver_test.cpp)
target_include_directories(MILPSolver_test PRIVATE ${MCFBlock_INCLUDE_DIR})
target_include_directories(MILPSolver_test PRIVATE ${MCFClass_INCLUDE_DIRS})
target_include_directories(MCFSolver_test PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(MCFSolver_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MCFSolver_test
                      SMS++
                      MCFBlock
                      MCFClass
                      MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )
