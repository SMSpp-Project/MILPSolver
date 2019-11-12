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

# ----- Settings ------------------------------------------------------------ #
include(CTest)
add_compile_definitions(HAVE_RELAX)

# Since we are using the block factory, objects from statically linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.

if (NOT BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(SMILPBlock REQUIRED)
find_package(MILPSolver REQUIRED)
find_package(MCFBlock REQUIRED)

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
                      SMS++::SMILPBlock)

add_executable(MCFSolver_test mcfsolver_test.cpp)
target_link_libraries(MCFSolver_test PRIVATE
                      SMS++::MILPSolver
                      SMS++::MCFBlock)
