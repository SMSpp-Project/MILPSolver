find_package(netCDF REQUIRED)
include_directories(${netCDF_INCLUDE_DIR})
link_directories(${netCDF_LIB_DIR})

# TODO: Write and use a FindSMILPBLock.cmake
set(SMS++_MAIN_DIR ${PROJECT_SOURCE_DIR}/../..)
include_directories(${SMS++_MAIN_DIR}/SMILPBlock)
include_directories(${SMS++_MAIN_DIR}/MCFBlock)
set(MCFClass_INCLUDE_DIR ${SMS++_MAIN_DIR}/MCFClass)
include_directories(${MCFClass_INCLUDE_DIR}/MCFClass
                    ${MCFClass_INCLUDE_DIR}/MCFClone
                    ${MCFClass_INCLUDE_DIR}/MCFCplex
                    ${MCFClass_INCLUDE_DIR}/MCFSimplex
                    ${MCFClass_INCLUDE_DIR}/OPTUtils
                    ${MCFClass_INCLUDE_DIR}/RelaxIV
                    )


add_executable(MILPSolver_test
               milpsolver_test.cpp
               )

# TODO: Linking doesn't work in standalone mode
target_link_libraries(MILPSolver_test
                      SMS++
                      SMILPBLock
                      MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )

# Options for choosing the MCF Solver
# TODO: Find how avoid duplication
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

add_executable(MCFSolver_test
               mcfsolver_test.cpp
               )

target_link_libraries(MCFSolver_test
                      SMS++
                      MCFBlock
                      MCFClass
                      MILPSolver
                      cplex-library
                      netcdf
                      netcdf-cxx4
                      )
