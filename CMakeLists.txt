# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MILPSolver                                     #
#                                                                             #
#    If Boost, Eigen3 and netCDF are installed properly,                      #
#    and SMS++ core library is either installed or already built              #
#    on this system using CMake, this should run just fine.                   #
#    Otherwise, edit CMakeConfig.txt to set the proper paths.                 #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8)
# ----- Project description and directories -----------------------------------
# TODO: set proper version
project(MILPSolver VERSION 1.0.0)

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

# This is to use the find modules we provide
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Requirements ----------------------------------------------------------
# This file contains HINTS supported by the find modules
include(CMakeConfig.txt OPTIONAL)

find_package(CPLEX)
find_package(SMSpp  REQUIRED)
find_package(Boost  REQUIRED)

# ----- Library ---------------------------------------------------------------
add_library(MILPSolver STATIC)
target_compile_features(MILPSolver PUBLIC cxx_std_17)

target_sources(MILPSolver PRIVATE
               ${PROJECT_SOURCE_DIR}/MILPSolver.cpp)

target_include_directories(MILPSolver PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>)

target_link_libraries(MILPSolver PUBLIC
                      SMS++::SMSpp
                      Boost::boost)

if (CPLEX_FOUND)
    message(STATUS "CPLEX found: ${CPLEX_STUDIO_DIR}")
    target_sources(MILPSolver PRIVATE ${PROJECT_SOURCE_DIR}/CPXMILPSolver.cpp)
    target_link_libraries(MILPSolver PUBLIC CPLEX::Cplex)
endif()

# Alias for testing
add_library(SMS++::MILPSolver ALIAS MILPSolver)

# ----- Install instructions --------------------------------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Install library
install(TARGETS MILPSolver
        EXPORT  MILPSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install export file *Targets.cmake
install(EXPORT MILPSolver-targets
        FILE   MILPSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${INSTALL_CONFIGDIR})

# Create *ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Create *Config.cmake
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MILPSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install the Config, ConfigVersion and FindCPLEX
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ----- Export from the build tree --------------------------------------------

export(EXPORT MILPSolver-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverTargets.cmake
       NAMESPACE SMS++::)

# Export the find module and the config file so other projects will use them
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
               COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindnetCDF.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindnetCDF.cmake
               COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeConfig.txt
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeConfig.txt
               COPYONLY)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MILPSolver)

# ----- Subdirectories --------------------------------------------------------
if (CPLEX_FOUND)
    # All tests use CPLEX
    add_subdirectory(test)
endif()
