# --------------------------------------------------------------------------- #
#    CMake find module for HiGHS                                              #
#                                                                             #
#    This module finds HiGHS include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(HiGHS [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        HiGHS_FOUND         - True if headers are found                      #
#        HiGHS_INCLUDE_DIRS  - Include directories                            #
#        HiGHS_LIBRARIES     - Libraries to be linked                         #
#        HiGHS_VERSION       - Version number                                 #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        HiGHS_ROOT          - Custom path to HiGHS                           #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        HiGHS::HiGHS                                                         #
#                                                                             #
#    This find module is provided because HiGHS does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                              Enrico Calandrini                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find HiGHS directories and lib suffixes ----------------------------- #
# Based on the OS generate:
# - a list of possible HiGHS directories
# - a list of possible lib suffixes to find the library

if (UNIX)
    if (APPLE)
        # macOS (usually /Library)
        set(HiGHS_DIRS /Library)
    else ()
        # Other Unix-based systems (usually /opt)
        set(HiGHS_DIRS /opt)
    endif ()
    set(HiGHS_LIB_PATH_SUFFIXES lib build/lib)
else ()
    # Windows (usually C:)
    set(HiGHS_DIRS "C:")
    set(HiGHS_LIB_PATH_SUFFIXES build/RELEASE/bin)
    set(HiGHS_LIB_PATH_SUFFIXES_DEBUG build/DEBUG/bin)
endif ()

# ----- Find the path to HiGHS ---------------------------------------------- #

foreach (dir ${HiGHS_DIRS})
    file(GLOB HiGHS_DIRS "${dir}/HiGHS")
    if (NOT HiGHS_ROOT IN_LIST HiGHS_DIRS)
        if (NOT "${HiGHS_ROOT}" STREQUAL "")
            message(STATUS "Specified HiGHS: ${HiGHS_ROOT} not found")
        endif ()
        list(SORT HiGHS_DIRS)
        list(REVERSE HiGHS_DIRS)
        if (HiGHS_DIRS)
            list(GET HiGHS_DIRS 0 HiGHS_ROOT)
            message(STATUS "Using HiGHS: ${HiGHS_ROOT}")
            break()
        else ()
            set(HiGHS_ROOT HiGHS_ROOT-NOTFOUND)
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
if (HiGHS_INCLUDE_DIR AND HiGHS_LIBRARY AND HiGHS_LIBRARY_DEBUG)
    set(HiGHS_FOUND TRUE)
else ()

    set(HiGHS_DIR ${HiGHS_ROOT})

    # ----- Find the HiGHS include directory -------------------------------- #
    # Note that find_path() creates a cache entry
    find_path(HiGHS_INCLUDE_DIR
              NAMES Highs.h interfaces/highs_c_api.h
              PATHS ${HiGHS_DIR}
              PATH_SUFFIXES include/highs src
              DOC "HiGHS include directory.")

    # ----- Find the HiGHS config include directory ------------------------- #
    # Note that find_path() creates a cache entry
    find_path(HiGHS_CONFIG_INCLUDE_DIR
              NAMES HConfig.h
              PATHS ${HiGHS_DIR}
              PATH_SUFFIXES include/highs build
              DOC "HiGHS config include directory.")

    if (UNIX)
        # ----- Find the HiGHS library -------------------------------------- #
        find_library(HiGHS_LIBRARY
                     NAMES highs
                     PATH_SUFFIXES ${HiGHS_LIB_PATH_SUFFIXES}
                     DOC "HiGHS library.")
        set(HiGHS_LIBRARY_DEBUG ${HiGHS_LIBRARY})
    elseif (NOT HiGHS_LIBRARY)

        # ----- Macro: find_win_HiGHS_library ------------------------------- #
        # On Windows the version is appended to the library name which cannot be
        # handled by find_library, so here a macro to search manually.
        macro(find_win_HiGHS_library var path_suffixes)
            foreach (s ${path_suffixes})
                file(GLOB HiGHS_LIBRARY_CANDIDATES "${HiGHS_DIR}/${s}/highs*.lib")
                if (HiGHS_LIBRARY_CANDIDATES)
                    list(GET HiGHS_LIBRARY_CANDIDATES 0 ${var})
                    break()
                endif ()
            endforeach ()
            if (NOT ${var})
                set(${var} NOTFOUND)
            endif ()
        endmacro ()

        # Library
        find_win_HiGHS_library(HiGHS_LIB "${HiGHS_LIB_PATH_SUFFIXES}")
        set(HiGHS_LIBRARY ${HiGHS_LIB})

        # Debug library
        find_win_HiGHS_library(HiGHS_LIB "${HiGHS_LIB_PATH_SUFFIXES_DEBUG}")
        set(HiGHS_LIBRARY_DEBUG ${HiGHS_LIB})
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (HiGHS_CONFIG_INCLUDE_DIR)
        file(STRINGS
                "${HiGHS_CONFIG_INCLUDE_DIR}/HConfig.h"
                _HiGHS_version_lines REGEX "#define HIGHS_VERSION_(MAJOR|MINOR|PATCH)")

        string(REGEX REPLACE ".*HIGHS_VERSION_MAJOR *\([0-9]*\).*" "\\1" _HiGHS_version_major "${_HiGHS_version_lines}")
        string(REGEX REPLACE ".*HIGHS_VERSION_MINOR *\([0-9]*\).*" "\\1" _HiGHS_version_minor "${_HiGHS_version_lines}")
        string(REGEX REPLACE ".*HIGHS_VERSION_PATCH *\([0-9]*\).*" "\\1" _HiGHS_version_patch "${_HiGHS_version_lines}")

        set(HiGHS_VERSION "${_HiGHS_version_major}.${_HiGHS_version_minor}.${_HiGHS_version_patch}")
        unset(_HiGHS_version_lines)
        unset(_HiGHS_version_major)
        unset(_HiGHS_version_minor)
        unset(_HiGHS_version_patch)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            HiGHS
            REQUIRED_VARS HiGHS_LIBRARY HiGHS_LIBRARY_DEBUG
                          HiGHS_INCLUDE_DIR HiGHS_CONFIG_INCLUDE_DIR
            VERSION_VAR HiGHS_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (HiGHS_FOUND)
    set(HiGHS_INCLUDE_DIRS "${HiGHS_INCLUDE_DIR}" "${HiGHS_CONFIG_INCLUDE_DIR}")
    set(HiGHS_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    # See: https://cmake.org/cmake/help/latest/module/CheckLibraryExists.html
    check_library_exists(m floor "" HAVE_LIBM)
    if (HAVE_LIBM)
        set(HiGHS_LINK_LIBRARIES ${HiGHS_LINK_LIBRARIES} m)
    endif ()

    if (UNIX)
        set(HiGHS_LINK_LIBRARIES ${HiGHS_LINK_LIBRARIES} dl)
    endif ()

    if (NOT TARGET HiGHS::HiGHS)
        add_library(HiGHS::HiGHS STATIC IMPORTED)
        set_target_properties(
                HiGHS::HiGHS PROPERTIES
                IMPORTED_LOCATION "${HiGHS_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${HiGHS_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${HiGHS_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${HiGHS_LINK_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(HiGHS_INCLUDE_DIR
                 HiGHS_CONFIG_INCLUDE_DIR
                 HiGHS_LIBRARY
                 HiGHS_LIBRARY_DEBUG
                 HiGHS_VERSION)

# --------------------------------------------------------------------------- #
