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
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)
smspp_derive_version(SMSPP_MODULE_VERSION)

project(MILPSolver
        VERSION ${SMSPP_MODULE_VERSION}
        DESCRIPTION "SMS++ MILPSolver module"
        HOMEPAGE_URL https://gitlab.com/smspp/milpsolver
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})

# C-safe form of modName for the preprocessor macros in <modName>Config.h
string(MAKE_C_IDENTIFIER "${modName}" modMacro)
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# An option is an ON/OFF user-settable cache variable
option(MILPSolver_USE_CPLEX "Whether MILPSolver will use CPLEX or not." ON)
option(MILPSolver_USE_SCIP "Whether MILPSolver will use SCIP or not." ON)
option(MILPSolver_USE_GUROBI "Whether MILPSolver will use GUROBI or not." ON)
option(MILPSolver_USE_HiGHS "Whether MILPSolver will use HiGHS or not." ON)
option(MILPSolver_USE_PIPS "Whether MILPSolver will use PIPS-IPM++ or not." ON)
option(MILPSolver_USE_DAQP "Whether MILPSolver will use DAQP or not." ON)
option(MILPSolver_USE_qpOASES "Whether MILPSolver will use qpOASES or not." ON)
option(MILPSolver_USE_Clarabel
        "Whether MILPSolver will use Clarabel or not." ON)
option(MILPSolver_USE_MOSEK "Whether MILPSolver will use MOSEK or not." ON)
option(MILPSolver_USE_ProxQP "Whether MILPSolver will use ProxQP or not." ON)
option(MILPSolver_USE_PIQP "Whether MILPSolver will use PIQP or not." ON)
option(MILPSolver_USE_OSQP "Whether MILPSolver will use OSQP or not." ON)
option(MILPSolver_USE_qpmad "Whether MILPSolver will use qpmad or not." ON)
option(MILPSolver_USE_quadprog
        "Whether MILPSolver will use quadprog or not." ON)

# A lot of debug prints
#if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#    add_compile_definitions(MILPSolver_DEBUG)
#endif ()

# ----- Requirements -------------------------------------------------------- #
# If it's not being called by the umbrella, we need to
# look for the system-installed libraries.
if (NOT hasParent)
    find_package(SMS++ REQUIRED)
endif ()

if (MILPSolver_USE_CPLEX)
    # Needed for CPXMILPSolver
    find_package(CPLEX 12.8.0)
    if (NOT CPLEX_FOUND)
        set(MILPSolver_USE_CPLEX OFF)
    endif ()
endif ()

if (MILPSolver_USE_SCIP)
    # Needed for SCIPMILPSolver
    find_package(SCIP 7.0.0)
    if (NOT SCIP_FOUND)
        set(MILPSolver_USE_SCIP OFF)
    endif ()
endif ()

if (MILPSolver_USE_GUROBI)
    # Needed for GUROBIMILPSolver
    find_package(GUROBI 10.0.0)
    if (NOT GUROBI_FOUND)
        set(MILPSolver_USE_GUROBI OFF)
    endif ()
endif ()

if (MILPSolver_USE_HiGHS)
    # Needed for HiGHSMILPSolver
    # 1.8.0 is the first HiGHS whose C API has all that HiGHSMILPSolver
    # calls: the callbacks are there since 1.6.0, Highs_setSparseSolution
    # since 1.8.0
    find_package(HiGHS 1.8.0)
    if (NOT HiGHS_FOUND)
        set(MILPSolver_USE_HiGHS OFF)
    endif ()
endif ()

if (MILPSolver_USE_PIPS)
    # Needed for PIPSMILPSolver
    # PIPS-IPM++ exposes no version number, so none can be checked here
    find_package(PIPS)
    if (NOT PIPS_FOUND)
        set(MILPSolver_USE_PIPS OFF)
    endif ()
endif ()

if (MILPSolver_USE_DAQP)
    # Needed for DAQPMILPSolver
    # DAQP exposes no version number: FindDAQP accepts 0.8.2 or later only,
    # by a macro of its headers [see FindDAQP.cmake]
    find_package(DAQP)
    if (NOT DAQP_FOUND)
        set(MILPSolver_USE_DAQP OFF)
    endif ()
endif ()

if (MILPSolver_USE_qpOASES)
    # Needed for QPOASESMILPSolver
    find_package(qpOASES 3.2)
    if (NOT qpOASES_FOUND)
        set(MILPSolver_USE_qpOASES OFF)
    endif ()
endif ()

if (MILPSolver_USE_Clarabel)
    # Needed for ClarabelMILPSolver
    # Clarabel exposes no version number: FindClarabel looks for clarabel.h,
    # which Clarabel.cpp has from 0.11.1 [see FindClarabel.cmake]
    find_package(Clarabel)
    if (NOT Clarabel_FOUND)
        set(MILPSolver_USE_Clarabel OFF)
    endif ()
endif ()

if (MILPSolver_USE_MOSEK)
    # Needed for MOSEKMILPSolver
    # 11.0.1 is the first with the error codes of writing an LP file
    find_package(MOSEK 11.0.1)
    if (NOT MOSEK_FOUND)
        set(MILPSolver_USE_MOSEK OFF)
    endif ()
endif ()

if (MILPSolver_USE_ProxQP)
    # Needed for ProxQPMILPSolver
    # ProxSuite installs its CMake configuration, which finds Eigen (and
    # Simde, if it has been built with vectorization); proxsuite_ROOT is a
    # hint to its location. 0.4 is the first with the HessianType of dense
    # ProxQP, and its version file accepts any later 0.x
    find_package(proxsuite 0.4 CONFIG)
    if (NOT proxsuite_FOUND)
        set(MILPSolver_USE_ProxQP OFF)
    endif ()
endif ()

if (MILPSolver_USE_PIQP)
    # Needed for PIQPMILPSolver
    # PIQP installs its CMake configuration, which finds Eigen; piqp_ROOT is
    # a hint to its location. Its version file only accepts the exact
    # version, hence 0.6 (the two-sided inequalities) is checked here
    find_package(piqp CONFIG)
    if (NOT piqp_FOUND OR piqp_VERSION VERSION_LESS 0.6)
        set(MILPSolver_USE_PIQP OFF)
    endif ()
endif ()

if (MILPSolver_USE_OSQP)
    # Needed for OSQPMILPSolver
    # OSQP installs its CMake configuration; osqp_ROOT is a hint to its
    # location. Its version file says 0.0.0 unless OSQP_VERSION is given at
    # its build, hence the API of OSQP 1.0 is checked by its header
    # osqp_api_functions.h, which OSQP 0.6 does not have
    find_package(osqp CONFIG)
    if (osqp_FOUND AND TARGET osqp::osqpstatic)
        set(OSQP_TARGET osqp::osqpstatic)
    elseif (osqp_FOUND AND TARGET osqp::osqp)
        set(OSQP_TARGET osqp::osqp)
    endif ()
    if (OSQP_TARGET)
        get_target_property(OSQP_INCLUDE_DIRS ${OSQP_TARGET}
                            INTERFACE_INCLUDE_DIRECTORIES)
        find_file(OSQP_API_FUNCTIONS osqp_api_functions.h
                  PATHS ${OSQP_INCLUDE_DIRS} NO_DEFAULT_PATH)
    endif ()
    if (NOT OSQP_TARGET OR NOT OSQP_API_FUNCTIONS)
        set(MILPSolver_USE_OSQP OFF)
    endif ()
endif ()

if (MILPSolver_USE_qpmad)
    # Needed for QPMADMILPSolver
    # qpmad is header-only and installs its CMake configuration, which finds
    # Eigen; qpmad_ROOT is a hint to its location. 1.5 is the first that
    # exports the target qpmad::qpmad, and its version file accepts any 1.x
    find_package(qpmad 1.5 CONFIG)
    if (NOT qpmad_FOUND)
        set(MILPSolver_USE_qpmad OFF)
    endif ()
endif ()

if (MILPSolver_USE_quadprog)
    # Needed for QuadProgMILPSolver
    # quadprog is under the GPL, hence a MILPSolver built with it is too; it
    # has no version number, and its qpgen2_ has not changed for decades
    find_package(quadprog)
    if (NOT quadprog_FOUND)
        set(MILPSolver_USE_quadprog OFF)
    endif ()
endif ()

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- RPATH --------------------------------------------------------------- #
# The installed libraries and executables look for the libraries of the same
# installation relative to their own location, so that the installed tree can
# be moved as a whole; the directories of the external libraries are kept.
include(GNUInstallDirs)
if (NOT CMAKE_INSTALL_RPATH)
    if (APPLE)
        set(rpathOrigin "@loader_path")
    else ()
        set(rpathOrigin "$ORIGIN")
    endif ()
    file(RELATIVE_PATH binToLib ${CMAKE_INSTALL_FULL_BINDIR}
            ${CMAKE_INSTALL_FULL_LIBDIR})
    set(CMAKE_INSTALL_RPATH "${rpathOrigin}" "${rpathOrigin}/${binToLib}")
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# The name of the shared library carries the version of its binary interface,
# which any minor release may change while the major version is 0, and then
# any major release.
if (PROJECT_VERSION_MAJOR EQUAL 0)
    set(soVersion ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
else ()
    set(soVersion ${PROJECT_VERSION_MAJOR})
endif ()
set_target_properties(${modName} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${soVersion})

# The classes of the module register themselves in the factory from a static
# initialiser, and a linker that drops what looks unused takes the
# registration away with it: the target tells whoever links it to keep the
# symbol that forces the module in, or the library as a whole where naming
# the symbol is not enough.
if (BUILD_SHARED_LIBS)
    if (MSVC)
        target_link_options(${modName} INTERFACE
                "/INCLUDE:SMSpp_force_load_${modName}")
    elseif (UNIX AND (NOT APPLE))
        # the ELF linker decides by the symbols a regular object asks for, so
        # naming the symbol is not enough and the library has to be kept
        target_link_options(${modName} INTERFACE "LINKER:--no-as-needed")
    elseif (APPLE)
        # ld64 keeps a dylib that provides a symbol to the link, one asked
        # for by -u included, even under -dead_strip_dylibs
        target_link_options(${modName} INTERFACE
                "LINKER:-u,_SMSpp_force_load_${modName}")
    endif ()
endif ()

# Ensure that the executables from /tools are run before adding the /include
# directory to generate the necessary *_defs.h and *_maps.h files of the correct
# solver version installed on the system, if they are not already present
add_dependencies(${modName} generate_headers)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE src/MILPSolver.cpp)

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
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(${modName} PUBLIC ${modNamespace}::SMS++)

if (MILPSolver_USE_CPLEX)
    message(STATUS "MILPSolver: CPLEX found, configuring CPXMILPSolver")
    target_sources(${modName} PRIVATE src/CPXMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC CPLEX::Cplex)
endif ()

if (MILPSolver_USE_SCIP)
    message(STATUS "MILPSolver: SCIP found, configuring SCIPMILPSolver")
    target_sources(${modName} PRIVATE src/SCIPMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC SCIP::SCIP)
endif ()

if (MILPSolver_USE_GUROBI)
    message(STATUS "MILPSolver: GUROBI found, configuring GRBMILPSolver")
    target_sources(${modName} PRIVATE src/GRBMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC GUROBI::Gurobi)
endif ()

if (MILPSolver_USE_HiGHS)
    message(STATUS "MILPSolver: HiGHS found, configuring HiGHSMILPSolver")
    target_sources(${modName} PRIVATE src/HiGHSMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC HiGHS::HiGHS)
endif ()

if (MILPSolver_USE_PIPS)
    message(STATUS "MILPSolver: PIPS-IPM++ found, configuring PIPSMILPSolver")
    target_sources(${modName} PRIVATE src/PIPSMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC PIPS::PIPS)
endif ()

if (MILPSolver_USE_DAQP)
    message(STATUS "MILPSolver: DAQP found, configuring DAQPMILPSolver")
    target_sources(${modName} PRIVATE src/DAQPMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC DAQP::DAQP)
endif ()

if (MILPSolver_USE_qpOASES)
    message(STATUS "MILPSolver: qpOASES found, configuring QPOASESMILPSolver")
    target_sources(${modName} PRIVATE src/QPOASESMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC qpOASES::qpOASES)
endif ()

if (MILPSolver_USE_Clarabel)
    message(STATUS
            "MILPSolver: Clarabel found, configuring ClarabelMILPSolver")
    target_sources(${modName} PRIVATE src/ClarabelMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC Clarabel::Clarabel)
endif ()

if (MILPSolver_USE_MOSEK)
    message(STATUS "MILPSolver: MOSEK found, configuring MOSEKMILPSolver")
    target_sources(${modName} PRIVATE src/MOSEKMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC MOSEK::MOSEK)
endif ()

if (MILPSolver_USE_ProxQP)
    message(STATUS "MILPSolver: ProxQP found, configuring ProxQPMILPSolver")
    target_sources(${modName} PRIVATE src/ProxQPMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC proxsuite::proxsuite)
endif ()

if (MILPSolver_USE_PIQP)
    message(STATUS "MILPSolver: PIQP found, configuring PIQPMILPSolver")
    target_sources(${modName} PRIVATE src/PIQPMILPSolver.cpp)
    # the compiled library of PIQP, if it has been built, otherwise its
    # header-only target
    if (TARGET piqp::piqp)
        target_link_libraries(${modName} PUBLIC piqp::piqp)
    else ()
        target_link_libraries(${modName} PUBLIC piqp::piqp_header_only)
    endif ()
endif ()

if (MILPSolver_USE_OSQP)
    message(STATUS "MILPSolver: OSQP found, configuring OSQPMILPSolver")
    target_sources(${modName} PRIVATE src/OSQPMILPSolver.cpp)
    # the static library of OSQP, if it has been built, otherwise its shared
    # one. OSQP exports its internal functions, some of which have the name
    # of functions of other libraries (e.g., csc_copy, also in HiGHS), and
    # on Linux the first library loaded gives the function to both; hence
    # the symbols of the static library are kept local to the binary
    target_link_libraries(${modName} PUBLIC ${OSQP_TARGET})
    if (OSQP_TARGET STREQUAL "osqp::osqpstatic" AND
        CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(${modName} PUBLIC
                            "LINKER:--exclude-libs,libosqpstatic.a")
    endif ()
endif ()

if (MILPSolver_USE_qpmad)
    message(STATUS "MILPSolver: qpmad found, configuring QPMADMILPSolver")
    target_sources(${modName} PRIVATE src/QPMADMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC qpmad::qpmad)
endif ()

if (MILPSolver_USE_quadprog)
    message(STATUS
            "MILPSolver: quadprog found, configuring QuadProgMILPSolver")
    target_sources(${modName} PRIVATE src/QuadProgMILPSolver.cpp)
    target_link_libraries(${modName} PUBLIC quadprog::quadprog)
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- C++20 module -------------------------------------------------------- #
# The module smspp.milpsolver (src/MILPSolver.cppm) exports the names of the
# headers and re-exports smspp; like the module smspp of the core, it is built
# when SMSpp_BUILD_MODULE is ON, which needs CMake 3.28 or later, the Ninja or
# Visual Studio generator and a compiler that scans the dependencies of the
# modules.
# Where these are missing, or a module it imports is not there, the option
# is turned off with a warning [see CheckModule.cmake of SMS++].
smspp_check_module(${modNamespace}::SMS++_module)
if (SMSpp_BUILD_MODULE)
    # the scan of the module dependencies, and $<LINK_ONLY> honoured also for
    # the target itself, which cmake_policy(VERSION 3.15) above leaves off
    cmake_policy(SET CMP0155 NEW)
    cmake_policy(SET CMP0131 NEW)

    add_library(${modName}_module)
    target_sources(${modName}_module PUBLIC
            FILE_SET CXX_MODULES
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
            FILES src/${modName}.cppm)
    target_compile_features(${modName}_module PUBLIC cxx_std_20)
    # the interface takes the include paths and the definitions of the
    # library, but not its compile options (e.g., the -fopenmp of a source,
    # which no header needs), so that it agrees with the module it imports
    target_include_directories(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_INCLUDE_DIRECTORIES>)
    target_compile_definitions(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_COMPILE_DEFINITIONS>)
    target_link_libraries(${modName}_module PUBLIC $<LINK_ONLY:${modName}>
            ${modNamespace}::SMS++_module)
    # the interface takes the header and the names of a backend only if it
    # is there [see the #if in src/MILPSolver.cppm]
    foreach (backend CPLEX SCIP GUROBI HiGHS PIPS DAQP qpOASES Clarabel MOSEK
             ProxQP PIQP OSQP qpmad quadprog)
        if (MILPSolver_USE_${backend})
            target_compile_definitions(${modName}_module PRIVATE
                    MILPSolver_MODULE_${backend})
        endif ()
    endforeach ()
    add_library(${modNamespace}::${modName}_module ALIAS ${modName}_module)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .dll
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so .dylib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}) # .lib

# Install the module, whose interface is compiled again by each importer
if (SMSpp_BUILD_MODULE)
    install(TARGETS ${modName}_module
            EXPORT ${modName}Targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            FILE_SET CXX_MODULES
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}/modules)
endif ()

# Install the headers
install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/MILPSolver.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

if (MILPSolver_USE_CPLEX)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "CPX*")
endif ()

if (MILPSolver_USE_SCIP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "SCIP*")
endif ()

if (MILPSolver_USE_GUROBI)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "GRB*")
endif ()

if (MILPSolver_USE_HiGHS)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "HiGHS*")
endif ()

if (MILPSolver_USE_PIPS)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "PIPS*")
endif ()

if (MILPSolver_USE_DAQP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "DAQP*")
endif ()

if (MILPSolver_USE_qpOASES)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "QPOASES*")
endif ()

if (MILPSolver_USE_Clarabel)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "Clarabel*")
endif ()

if (MILPSolver_USE_MOSEK)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "MOSEK*")
endif ()

if (MILPSolver_USE_ProxQP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "ProxQP*")
endif ()

if (MILPSolver_USE_PIQP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "PIQP*")
endif ()

if (MILPSolver_USE_OSQP)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "OSQP*")
endif ()

if (MILPSolver_USE_qpmad)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "QPMAD*")
endif ()

if (MILPSolver_USE_quadprog)
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}
            FILES_MATCHING PATTERN "QuadProg*")
endif ()

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
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

# ----- Add the build tree to RPATH------------------------------------------ #
# Add the binary directory to RPATH, so executables that use
# dynamic loading will look into it.
if (hasParent)
    set(SMSPP_BUILD_RPATH
            "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
endif ()

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html

# Using this with the umbrella project is too messy.
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
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
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
