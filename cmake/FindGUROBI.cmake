# --------------------------------------------------------------------------- #
#    CMake find module for GUROBI                                             #
#                                                                             #
#    This module finds GUROBI include directories and libraries.              #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(GUROBI [version] [EXACT] [REQUIRED])                    #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        GUROBI_FOUND         - True if headers are found                     #
#        GUROBI_INCLUDE_DIRS  - Include directories                           #
#        GUROBI_LIBRARIES     - Libraries to be linked                        #
#        GUROBI_VERSION       - Version number                                #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        GUROBI_ROOT          - Custom path to GUROBI                         #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        GUROBI::Gurobi                                                       #
#                                                                             #
#    This find module is provided because GUROBI does not provide             #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find Gurobi directories and lib suffixes ---------------------------- #
# Based on the OS generate:
# - a list of possible Gurobi directories
# - a list of possible lib suffixes to find the library

if (UNIX)
    if (APPLE)
        # macOS (usually /Library)
        set(GUROBI_DIRS /Library)
    else ()
        # Other Unix-based systems (usually /opt)
        set(GUROBI_DIRS /opt)
    endif ()
else ()
    # Windows (usually C:)
    set(GUROBI_DIRS "C:")
endif ()
set(GUROBI_LIB_PATH_SUFFIXES lib)

# ----- Find the path to GUROBI --------------------------------------------- #

foreach (dir ${GUROBI_DIRS})
    file(GLOB GUROBI_DIRS "${dir}/gurobi*")
    if (NOT IS_DIRECTORY "${GUROBI_ROOT}")
        if (NOT "${GUROBI_ROOT}" STREQUAL "")
            message(STATUS "Specified Gurobi: ${GUROBI_ROOT} not found")
        endif ()
        list(SORT GUROBI_DIRS)
        list(REVERSE GUROBI_DIRS)
        if (GUROBI_DIRS)
            list(GET GUROBI_DIRS 0 GUROBI_ROOT)
            message(STATUS "Using Gurobi: ${GUROBI_ROOT}")
            break()
        else ()
            set(GUROBI_ROOT GUROBI_ROOT-NOTFOUND)
        endif ()
    else ()
        break()
    endif ()
endforeach ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

# Check if already in cache
if (GUROBI_INCLUDE_DIR AND GUROBI_LIBRARY AND GUROBI_LIBRARY_DEBUG)
    set(GUROBI_FOUND TRUE)
else ()

    if (UNIX)
        if (APPLE)
            set(GUROBI_DIR ${GUROBI_ROOT}/macos_universal2)
        else ()
            set(GUROBI_DIR ${GUROBI_ROOT}/linux64)
        endif ()
    else () # Windows
        if (ARCH STREQUAL "x64")
            set(GUROBI_DIR ${GUROBI_ROOT}/win64)
        elseif (ARCH STREQUAL "x86")
            set(GUROBI_DIR ${GUROBI_ROOT}/win32)
        endif ()
    endif ()

    # ----- Find the GUROBI include directory ------------------------------- #
    # Note that find_path() creates a cache entry
    find_path(GUROBI_INCLUDE_DIR
              NAMES gurobi_c.h
              PATHS ${GUROBI_DIR}/include
              DOC "GUROBI include directory.")

    if (UNIX)
        # ----- Find the GUROBI library ------------------------------------- #
        find_library(GUROBI_LIBRARY
                     NAMES gurobi gurobi100
                     PATHS ${GUROBI_DIR}
                     PATH_SUFFIXES ${GUROBI_LIB_PATH_SUFFIXES}
                     DOC "GUROBI library.")
        set(GUROBI_LIBRARY_DEBUG ${GUROBI_LIBRARY})
    elseif (NOT GUROBI_LIBRARY)

        # ----- Macro: find_win_gurobi_library ------------------------------ #
        # On Windows the version is appended to the library name which cannot be
        # handled by find_library, so here a macro to search manually.
        macro(find_win_gurobi_library var path_suffixes)
            foreach (s ${path_suffixes})
                file(GLOB GUROBI_LIBRARY_CANDIDATES "${GUROBI_DIR}/${s}/gurobi*.lib")
                if (GUROBI_LIBRARY_CANDIDATES)
                    list(GET GUROBI_LIBRARY_CANDIDATES 0 ${var})
                    break()
                endif ()
            endforeach ()
            if (NOT ${var})
                set(${var} NOTFOUND)
            endif ()
        endmacro ()

        # Library
        find_win_gurobi_library(GUROBI_LIB "${GUROBI_LIB_PATH_SUFFIXES}")
        set(GUROBI_LIBRARY ${GUROBI_LIB})

        # Debug library
        find_win_gurobi_library(GUROBI_LIB "${GUROBI_LIB_PATH_SUFFIXES_DEBUG}")
        set(GUROBI_LIBRARY_DEBUG ${GUROBI_LIB})
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (GUROBI_INCLUDE_DIR)
        file(STRINGS
                "${GUROBI_INCLUDE_DIR}/gurobi_c.h"
                _gurobi_version_lines REGEX "#define GRB_VERSION_(MAJOR|MINOR|TECHNICAL)")

        string(REGEX REPLACE ".*GRB_VERSION_MAJOR *\([0-9]*\).*" "\\1" _gurobi_version_major "${_gurobi_version_lines}")
        string(REGEX REPLACE ".*GRB_VERSION_MINOR *\([0-9]*\).*" "\\1" _gurobi_version_minor "${_gurobi_version_lines}")
        string(REGEX REPLACE ".*GRB_VERSION_TECHNICAL *\([0-9]*\).*" "\\1" _gurobi_version_technical "${_gurobi_version_lines}")

        set(GUROBI_VERSION "${_gurobi_version_major}.${_gurobi_version_minor}.${_gurobi_version_technical}")
        unset(_gurobi_version_lines)
        unset(_gurobi_version_major)
        unset(_gurobi_version_minor)
        unset(_gurobi_version_technical)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            GUROBI
            REQUIRED_VARS GUROBI_LIBRARY GUROBI_LIBRARY_DEBUG GUROBI_INCLUDE_DIR
            VERSION_VAR GUROBI_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (GUROBI_FOUND)
    set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
    set(GUROBI_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    # See: https://cmake.org/cmake/help/latest/module/CheckLibraryExists.html
    check_library_exists(m floor "" HAVE_LIBM)
    if (HAVE_LIBM)
        set(GUROBI_LINK_LIBRARIES ${GUROBI_LINK_LIBRARIES} m)
    endif ()

    if (UNIX)
        set(GUROBI_LINK_LIBRARIES ${GUROBI_LINK_LIBRARIES} dl)
    endif ()

    if (NOT TARGET GUROBI::Gurobi)
        add_library(GUROBI::Gurobi STATIC IMPORTED)
        set_target_properties(
                GUROBI::Gurobi PROPERTIES
                IMPORTED_LOCATION "${GUROBI_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${GUROBI_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${GUROBI_LINK_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(GUROBI_INCLUDE_DIR
                 GUROBI_LIBRARY
                 GUROBI_LIBRARY_DEBUG
                 GUROBI_VERSION)

# --------------------------------------------------------------------------- #
