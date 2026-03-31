# --------------------------------------------------------------------------- #
#    CMake find module for SCIP                                               #
#                                                                             #
#    This module finds SCIP include directories and libraries.                #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(SCIP [version] [EXACT] [REQUIRED])                      #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        SCIP_FOUND         - True if headers are found                       #
#        SCIP_INCLUDE_DIRS  - Include directories                             #
#        SCIP_LIBRARIES     - Libraries to be linked                          #
#        SCIP_VERSION       - Version number                                  #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        SCIP_ROOT          - Custom path to SCIP                             #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        SCIP::SCIP                                                           #
#                                                                             #
#    This find module is provided because SCIP does not provide               #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                 Donato Meoli                                #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find SCIP directories and lib suffixes ----------------------------- #
# Based on the OS generate:
# - a list of possible SCIP directories
# - a list of possible lib suffixes to find the library

if (UNIX)
    if (APPLE)
        # macOS (usually /Library)
        set(SCIP_DIRS /Library)
    else ()
        # Other Unix-based systems (usually /opt)
        set(SCIP_DIRS /opt)
    endif ()
elseif (WIN32)
    # Windows (usually C:/Program Files)
    set(SCIP_DIRS "C:/Program Files")
    if (ARCH STREQUAL "x86")
        set(SCIP_DIRS "C:/Program Files (x86)" ${SCIP_DIRS})
    endif ()
endif ()
set(SCIP_LIB_PATH_SUFFIXES lib)

# ----- Find the path to SCIP ---------------------------------------------- #

foreach (dir ${SCIP_DIRS})
    file(GLOB SCIP_DIRS "${dir}/scip*")
    if (NOT SCIP_ROOT IN_LIST SCIP_DIRS)
        message(STATUS "Specified SCIP: ${SCIP_ROOT} not found")
        list(SORT SCIP_DIRS)
        list(REVERSE SCIP_DIRS)
        if (SCIP_DIRS)
            list(GET SCIP_DIRS 0 SCIP_ROOT)
            message(STATUS "Using SCIP: ${SCIP_ROOT}")
            break()
        else ()
            set(SCIP_ROOT SCIP_ROOT-NOTFOUND)
        endif ()
    else ()
        break()
    endif ()
endforeach ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

if (UNIX)
    find_package(TBB QUIET)
endif ()

# Check if already in cache
if (WIN32)
    if (SCIP_INCLUDE_DIR AND SCIP_LIBRARY AND SCIP_DLL AND TBB_DLL AND SCIP_VERSION)
        set(SCIP_FOUND TRUE)
    endif ()
else ()
    if (SCIP_INCLUDE_DIR AND SCIP_LIBRARY AND SCIP_VERSION)
        set(SCIP_FOUND TRUE)
    endif ()
endif ()

if (NOT SCIP_FOUND)

    # ----- Find the SCIP include directory --------------------------------- #
    find_path(SCIP_INCLUDE_DIR
            NAMES scip/scip.h
            PATHS ${SCIP_ROOT}/include
            DOC "SCIP include directory.")

    # ----- Find the SCIP library ------------------------------------------- #
    find_library(SCIP_LIBRARY
            NAMES scip
            PATHS ${SCIP_ROOT}
            PATH_SUFFIXES ${SCIP_LIB_PATH_SUFFIXES}
            DOC "SCIP library.")

    # ----- Find the SCIP runtime DLL on Windows ---------------------------- #
    if (WIN32)
        find_file(SCIP_DLL
                NAMES libscip.dll scip.dll
                PATHS ${SCIP_ROOT}
                PATH_SUFFIXES bin
                DOC "SCIP runtime DLL.")
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (SCIP_INCLUDE_DIR)
        file(STRINGS
                "${SCIP_INCLUDE_DIR}/scip/config.h"
                _SCIP_version_lines REGEX "#define SCIP_VERSION_(MAJOR|MINOR|PATCH)")

        string(REGEX REPLACE ".*SCIP_VERSION_MAJOR *\([0-9]*\).*" "\\1" _SCIP_version_major "${_SCIP_version_lines}")
        string(REGEX REPLACE ".*SCIP_VERSION_MINOR *\([0-9]*\).*" "\\1" _SCIP_version_minor "${_SCIP_version_lines}")
        string(REGEX REPLACE ".*SCIP_VERSION_PATCH *\([0-9]*\).*" "\\1" _SCIP_version_patch "${_SCIP_version_lines}")

        set(SCIP_VERSION "${_SCIP_version_major}.${_SCIP_version_minor}.${_SCIP_version_patch}")
        unset(_SCIP_version_lines)
        unset(_SCIP_version_major)
        unset(_SCIP_version_minor)
        unset(_SCIP_version_patch)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    if (WIN32)
        find_package_handle_standard_args(
                SCIP
                REQUIRED_VARS SCIP_LIBRARY SCIP_DLL SCIP_INCLUDE_DIR
                VERSION_VAR SCIP_VERSION)
    else ()
        find_package_handle_standard_args(
                SCIP
                REQUIRED_VARS SCIP_LIBRARY SCIP_INCLUDE_DIR
                VERSION_VAR SCIP_VERSION)
    endif ()
endif ()

# ----- Export the target --------------------------------------------------- #
if (SCIP_FOUND)
    set(SCIP_INCLUDE_DIRS ${SCIP_INCLUDE_DIR})
    set(SCIP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    if (TARGET TBB::tbb)
        set(SCIP_LIBRARIES ${SCIP_LIBRARIES} TBB::tbb)
    endif ()

    if (UNIX)
        set(SCIP_LIBRARIES ${SCIP_LIBRARIES} dl)
    endif ()

    if (NOT TARGET SCIP::SCIP)
        if (WIN32)
            add_library(SCIP::SCIP SHARED IMPORTED)
            set_target_properties(
                    SCIP::SCIP PROPERTIES
                    IMPORTED_IMPLIB "${SCIP_LIBRARY}"
                    IMPORTED_LOCATION "${SCIP_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${SCIP_LIBRARIES}")
        else ()
            add_library(SCIP::SCIP UNKNOWN IMPORTED)
            set_target_properties(
                    SCIP::SCIP PROPERTIES
                    IMPORTED_LOCATION "${SCIP_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${SCIP_LIBRARIES}")
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
if (WIN32)
    mark_as_advanced(SCIP_INCLUDE_DIR
            SCIP_LIBRARY
            SCIP_DLL
            SCIP_VERSION)
else ()
    mark_as_advanced(SCIP_INCLUDE_DIR
            SCIP_LIBRARY
            SCIP_VERSION)
endif ()

# --------------------------------------------------------------------------- #
