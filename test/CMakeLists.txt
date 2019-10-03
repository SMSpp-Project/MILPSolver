# ----- Fetch GoogleTest ---------------------------------------------------- #
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.8.1
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

# ----- Requirements -------------------------------------------------------- #
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Trying to find SMILPBLock while configuring all at once
    # gave some errors (probably race conditions).
    find_package(SMILPBLock REQUIRED)
    find_package(MCFBlock REQUIRED)
    find_package(netCDF REQUIRED)
endif ()

# ----- CPXMILPSolver_test -------------------------------------------------- #
add_executable(CPXMILPSolver_test tests_CPXMILPSolver.cpp)

target_link_libraries(CPXMILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::SMILPBLock
                      gtest)

add_test(
        NAME CPXMILPSolver_test
        COMMAND CPXMILPSolver_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# ----- MCFMILPSolvers_test ------------------------------------------------- #
add_executable(MCFMILPSolvers_test tests_MCFMILPSolvers.cpp)

target_link_libraries(MCFMILPSolvers_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::MCFBlock
                      NetCDF::NetCDFCxx
                      gtest)

add_test(
        NAME MCFMILPSolvers_test
        COMMAND MCFMILPSolvers_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


# ----- OLD TESTS ----------------------------------------------------------- #
add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::SMILPBLock)

add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::MCFBlock
                      NetCDF::NetCDFCxx)
