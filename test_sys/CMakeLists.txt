# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MILPSolver system tests                        #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11...3.15)
# TODO: set proper version
project(MILPSolver_test_sys VERSION 1.0.0 LANGUAGES CXX)

# TODO: Everything related to deeper deps was disabled to check if it works

# ----- Settings ------------------------------------------------------------ #
enable_testing()
set(CMAKE_BUILD_TYPE "Debug")
add_compile_options( "$<$<CONFIG:DEBUG>:-fno-inline>")
add_compile_definitions(HAVE_RELAX)

# ----- Requirements -------------------------------------------------------- #
find_package(SMILPBlock REQUIRED)
find_package(MCFBlock REQUIRED)
find_package(MILPSolver REQUIRED)
# TODO: Here
#find_package(netCDF REQUIRED)

# ----- Fetch GoogleTest ---------------------------------------------------- #
include(FetchContent)
fetchcontent_declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.8.1
)

fetchcontent_getproperties(googletest)
if (NOT googletest_POPULATED)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_MACOSX_RPATH 1)
    endif ()
    fetchcontent_populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

# ----- CPXMILPSolver_test -------------------------------------------------- #
add_executable(CPXMILPSolver_test tests_CPXMILPSolver.cpp)

target_link_libraries(CPXMILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::SMILPBlock
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
#                      NetCDF::NetCDFCxx # TODO: Here
                      gtest
                      )

add_test(
        NAME MCFMILPSolvers_test
        COMMAND MCFMILPSolvers_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


# ----- OLD TESTS ----------------------------------------------------------- #
add_executable(MILPSolver_test milpsolver_test.cpp)

target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::SMILPBlock
                      )

add_executable(MCFSolver_test mcfsolver_test.cpp)

target_link_libraries(MCFSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::MCFBlock
#                      NetCDF::NetCDFCxx # TODO: Here
                      )
