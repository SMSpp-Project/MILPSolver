# ----- Project description and directories ------------------------------------
cmake_minimum_required(VERSION 3.8)
project(MILPSolver
        VERSION 1.0.0 # TODO: set proper version
        HOMEPAGE_URL "https://gitlab.com/niccolo/milpsolver"
        DESCRIPTION "MILPSolver for SMS++")

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(CPLEX REQUIRED)
find_package(SMSpp REQUIRED)

# ----- Library ----------------------------------------------------------------
add_library(MILPSolver OBJECT) # TODO: Find out why STATIC doesn't work
target_compile_features(MILPSolver PUBLIC cxx_std_17)

target_sources(MILPSolver PRIVATE
               ${PROJECT_SOURCE_DIR}/MILPSolver.cpp
               ${PROJECT_SOURCE_DIR}/CPXMILPSolver.cpp
               )

target_include_directories(MILPSolver PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>)

target_link_libraries(MILPSolver PUBLIC
                      SMS++::SMSpp
                      Boost::boost
                      Eigen3::Eigen
                      cplex-library)

# Alias for testing
add_library(SMS++::MILPSolver ALIAS MILPSolver)

# ----- Install instructions ---------------------------------------------------
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

# ----- Export from the build tree ---------------------------------------------
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
               COPYONLY)

export(EXPORT MILPSolver-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverTargets.cmake
       NAMESPACE SMS++::)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MILPSolver)

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
