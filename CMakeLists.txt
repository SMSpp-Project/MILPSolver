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
project(MILPSolver
        VERSION 0.3.0
        DESCRIPTION "SMS++ MILPSolver module"
        HOMEPAGE_URL https://gitlab.com/smspp/milpsolver
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace "SMS++")
set(modHeader "${modName}Config.h")
set(modTargets "${modName}Targets")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
#include(CTest)

# An option is an ON/OFF user-settable cache variable
option(${modName}_USE_CPLEX "Whether ${modName} will use CPLEX or not." ON)
option(${modName}_USE_SCIP "Whether ${modName} will use SCIP or not." ON)

# A lot of debug prints
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_compile_definitions(MILPSOLVER_DEBUG)
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
configure_file(cmake/${modHeader}.in ${modHeader})

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
               src/MILPSolver.cpp)

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
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(${modName} PUBLIC
                      ${modNamespace}::SMSpp)

if (${modName}_USE_CPLEX)
    message(STATUS "${modName}: CPLEX found, configuring CPX${modName}")

    target_sources(${modName} PRIVATE
                   src/CPX${modName}.cpp)

    target_link_libraries(${modName} PUBLIC
                          CPLEX::Cplex)
endif ()

if (${modName}_USE_SCIP)
    message(STATUS "${modName}: SCIP found, configuring SCIP${modName}")

    target_sources(${modName} PRIVATE
                   src/SCIP${modName}.cpp)

    # SCIP configuration file doesn't define a target
    target_link_libraries(${modName} PUBLIC ${SCIP_LIBRARIES})
    target_include_directories(${modName} PRIVATE ${SCIP_INCLUDE_DIRS})
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modTargets}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MILPSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

if (${modName}_USE_CPLEX)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "CPX*")
endif ()

if (${modName}_USE_SCIP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "SCIP*")
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modTargets}
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Shipping custom find modules should not be encouraged because third
# party developers should provide for them.
file(GLOB findMods ${CMAKE_CURRENT_LIST_DIR}/cmake/Find*)
install(FILES ${findMods}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

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
    export(EXPORT ${modTargets}
           FILE ${CMAKE_CURRENT_BINARY_DIR}/${modTargets}.cmake
           NAMESPACE ${modNamespace}::)

    # Export the find modules for the local configuration file to use them.
    foreach (findMod ${findMods})
        configure_file(${findMod} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
    endforeach ()

    # Add the binary directory to RPATH, so executables that use
    # dynamic loading will look into it.
    if (hasParent)
        set(SMSPP_BUILD_RPATH
            "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
    endif ()

    # Export the package in user's package registry
    export(PACKAGE ${modName})
endif ()

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html

# Using this with the umbrella project is too messy.
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Niccolò Iardella")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
        \\\\.git/
        \\\\.gitignore$
        \\\\.gitattributes$)
    include(CPack)
endif ()

# --------------------------------------------------------------------------- #
