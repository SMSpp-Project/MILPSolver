# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MILPSolver                                     #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
# TODO: set proper version
project(MILPSolver VERSION 1.0.0 LANGUAGES CXX)

# ----- Requirements -------------------------------------------------------- #
# find_package() is used to find libraries and automatically configure them.
# It requires either a package configuration module (***Config.cmake) or a
# find module (Find***.cmake), usually provided by libraries.

# We can provide a settings file with -DSETTINGS="filename"
if ((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND SETTINGS)
    message(STATUS "Using CMake custom settings in ${SETTINGS}")
    include(${SETTINGS} OPTIONAL)
endif ()

# This is to use the find modules we provide
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Boost 1.71.0 REQUIRED)
find_package(SMSpp REQUIRED)
find_package(CPLEX)

# ----- Library ------------------------------------------------------------- #
add_library(MILPSolver)
target_compile_features(MILPSolver PUBLIC cxx_std_17)

target_sources(MILPSolver PRIVATE src/MILPSolver.cpp)

target_include_directories(
        MILPSolver PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(MILPSolver PUBLIC
                      SMS++::SMSpp
                      Boost::boost)

if (CPLEX_FOUND)
    message(STATUS "CPLEX found: ${CPLEX_STUDIO_DIR}")
    target_sources(MILPSolver PRIVATE src/CPXMILPSolver.cpp)
    target_link_libraries(MILPSolver PUBLIC CPLEX::Cplex)
endif ()

# Alias for testing
add_library(SMS++::MILPSolver ALIAS MILPSolver)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS MILPSolver
        EXPORT MILPSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT MILPSolver-targets
        FILE MILPSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MILPSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Shipping custom find modules should not be encouraged because third
# party developers should provide them.
# We ship FindCPLEX because CPLEX does not on its own.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)


# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)

    export(EXPORT MILPSolver-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverTargets.cmake
           NAMESPACE SMS++::)

    # Export also the find modules for the local config file to use it
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
                   COPYONLY)

    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindnetCDF.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindnetCDF.cmake
                   COPYONLY)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE MILPSolver)
endif ()

# ----- Subdirectories ------------------------------------------------------ #

if (CPLEX_FOUND AND CMAKE_TESTING_ENABLED)
    # All tests use CPLEX
    add_subdirectory(test)
endif ()
