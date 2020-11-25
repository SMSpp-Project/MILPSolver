# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MILPSolver                                     #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
project(MILPSolver VERSION 0.3.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# An option is an ON/OFF user-settable cache variable
option(MILPSolver_USE_CPLEX "Whether MILPSolver will use CPLEX" ON)
option(MILPSolver_USE_SCIP "Whether MILPSolver will use SCIP" ON)

# A lot of debug prints
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    #    add_compile_definitions(MILPSOLVER_DEBUG)
endif ()

# ----- Requirements -------------------------------------------------------- #
# This adds the cmake directory to the find module search paths,
# allowing us to use our find modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(SMSpp REQUIRED)

if (MILPSolver_USE_CPLEX)
    find_package(CPLEX 12.8) # Needed for CPXMILPSolver
endif ()

if (MILPSolver_USE_SCIP)
    find_package(SCIP) # Needed for SCIPMILPSolver
endif ()

# This removes from the search paths the cmake directory we added above.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# Update cache variables
# FIXME: Not sure if this is a best practice
if (NOT CPLEX_FOUND)
    set(MILPSolver_USE_CPLEX OFF)
endif ()

if (NOT SCIP_FOUND)
    set(MILPSolver_USE_SCIP OFF)
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/MILPSolverConfig.h.in MILPSolverConfig.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(MILPSolver)
target_compile_features(MILPSolver PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
target_sources(MILPSolver PRIVATE src/MILPSolver.cpp)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
target_include_directories(
        MILPSolver PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(MILPSolver PUBLIC
                      SMS++::SMSpp)

if (MILPSolver_USE_CPLEX)
    message(STATUS "MILPSolver: CPLEX found, configuring CPXMILPSolver")
    target_sources(MILPSolver PRIVATE src/CPXMILPSolver.cpp)
    target_link_libraries(MILPSolver PUBLIC CPLEX::Cplex)
endif ()

if (MILPSolver_USE_SCIP)
    message(STATUS "MILPSolver: SCIP found, configuring SCIPMILPSolver")
    target_sources(MILPSolver PRIVATE src/SCIPMILPSolver.cpp)
    # SCIP configuration file doesn't define a target
    target_link_libraries(MILPSolver PUBLIC ${SCIP_LIBRARIES})
    target_include_directories(MILPSolver PRIVATE ${SCIP_INCLUDE_DIRS})
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(SMS++::MILPSolver ALIAS MILPSolver)

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS MILPSolver
        EXPORT MILPSolver-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MILPSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

if (MILPSolver_USE_CPLEX)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++
            FILES_MATCHING PATTERN "CPX*")
endif ()

if (MILPSolver_USE_SCIP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++
            FILES_MATCHING PATTERN "SCIP*")
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/MILPSolverConfig.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT MILPSolver-targets
        FILE MILPSolverTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MILPSolver)

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MILPSolverConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MILPSolver)

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MILPSolver)

# Shipping custom find modules should not be encouraged because third
# party developers should provide them.
# We ship FindCPLEX because CPLEX does not on its own.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MILPSolver)

# ----- Export from the build tree ------------------------------------------ #
# The following commands are needed so other projects can find the library
# in a local build directory. They use CMake's User Package Registry
# (located in ~/.cmake/registry) to do so.
# See: https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html
# They don't affect the installation.

# CMAKE_EXPORT_PACKAGE_REGISTRY is OFF by default.
# It may be ON if the project is included by the SMS++ umbrella project.
if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    # Generate the target import file, that allows other
    # CMake projects to import the target.
    export(EXPORT MILPSolver-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MILPSolverTargets.cmake
           NAMESPACE SMS++::)

    # Export the find modules for the local configuration file to use them.
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCPLEX.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindCPLEX.cmake
                   COPYONLY)

    # Add the binary directory to RPATH, so executables that use
    # dynamic loading will look into it.
    get_directory_property(hasParent PARENT_DIRECTORY)
    if (hasParent)
        set(SMSPP_BUILD_RPATH "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
    endif ()

    # Export the package in user's package registry
    export(PACKAGE MILPSolver)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()
