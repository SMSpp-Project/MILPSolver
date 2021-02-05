# --------------------------------------------------------------------------- #
#    CMake file for unit tests                                                #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# ----- GoogleTest unit tests ----------------------------------------------- #
# The following tests are built using the GoogleTest framework,
# see: https://github.com/google/googletest
# GoogleTest is automatically fetched and built using CMake's FetchContent,
# see: https://cmake.org/cmake/help/latest/module/FetchContent.html
include(FetchContent)

# Declare GoogleTest repository
fetchcontent_declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.10.0)

# Fetch GoogleTest, if not already fetched
fetchcontent_getproperties(googletest)
if (NOT googletest_POPULATED)
    # This is so GoogleTest doesn't get installed alongside the library.
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

    # See: https://cmake.org/cmake/help/latest/prop_tgt/MACOSX_RPATH.html
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_MACOSX_RPATH 1)
    endif ()

    fetchcontent_populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

# ----- MILPSolver_test ----------------------------------------------------- #
add_executable(MILPSolver_test test.cpp SimpleMILPBlock.cpp)
target_link_libraries(MILPSolver_test PRIVATE
                      SMS++::MILPSolver
                      gtest)

add_test(NAME MILPSolver_test
         COMMAND MILPSolver_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# --------------------------------------------------------------------------- #
