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
#                               Enrico Calandrini                             #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
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
        set(SCIP_DIRS /usr/local/include)
    endif ()
else ()
    # Windows (usually C:)
    set(SCIP_DIRS "C:")
endif ()
set(SCIP_LIB_PATH_SUFFIXES lib)

# ----- Find the path to SCIP ---------------------------------------------- #
foreach (dir ${SCIP_DIRS})
    file(GLOB SCIP_DIRS "${dir}")
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

# Check if already in cache
if (SCIP_INCLUDE_DIR AND SCIP_LIBRARY AND SCIP_LIBRARY_DEBUG)
    set(SCIP_FOUND TRUE)
else ()

    if (UNIX)
        set(SCIP_DIR ${SCIP_ROOT})
    else () # Windows
        if (ARCH MATCHES x64)
            set(SCIP_DIR ${SCIP_ROOT}/win64)
        elseif (ARCH MATCHES x86)
            set(SCIP_DIR ${SCIP_ROOT}/win32)
        endif ()
    endif ()

    # ----- Find the SCIP include directory -------------------------------- #
    # Note that find_path() creates a cache entry
    find_path(SCIP_INCLUDE_DIR
              NAMES scip/scip.h
              PATHS ${SCIP_DIR}
              PATH_SUFFIXES SCIP
              DOC "SCIP include directory.")
    set(SCIP_INCLUDE_DIR /usr/local/include)

    if (UNIX)
        # ----- Find the SCIP library -------------------------------------- #
        # Note that find_library() creates a cache entry
        find_library(SCIP_LIBRARY
                     NAMES scip
                     PATH_SUFFIXES ${SCIP_LIB_PATH_SUFFIXES}
                     DOC "SCIP library.")
        set(SCIP_LIBRARY_DEBUG ${SCIP_LIBRARY})
    elseif (NOT SCIP_LIBRARY)

        # ----- Macro: find_win_SCIP_library ------------------------------- #
        # On Windows the version is appended to the library name which cannot be
        # handled by find_library, so here a macro to search manually.
        macro(find_win_SCIP_library var path_suffixes)
            foreach (s ${path_suffixes})
                file(GLOB SCIP_LIBRARY_CANDIDATES "${SCIP_DIR}/${s}/SCIP*.lib")
                if (SCIP_LIBRARY_CANDIDATES)
                    list(GET SCIP_LIBRARY_CANDIDATES 0 ${var})
                    break()
                endif ()
            endforeach ()
            if (NOT ${var})
                set(${var} NOTFOUND)
            endif ()
        endmacro()

        # Library
        find_win_SCIP_library(SCIP_LIB "${SCIP_LIB_PATH_SUFFIXES}")
        set(SCIP_LIBRARY ${SCIP_LIB})

        # Debug library
        find_win_SCIP_library(SCIP_LIB "${SCIP_LIB_PATH_SUFFIXES_DEBUG}")
        set(SCIP_LIBRARY_DEBUG ${SCIP_LIB})

        # DLL
        if (SCIP_LIBRARY MATCHES ".*/(SCIP.*)\\.lib")
            file(GLOB SCIP_DLL_ "${SCIP_DIR}/bin/${CMAKE_MATCH_1}.dll")
            set(SCIP_DLL ${SCIP_DLL_})
        endif ()
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            SCIP
            REQUIRED_VARS SCIP_LIBRARY SCIP_LIBRARY_DEBUG SCIP_INCLUDE_DIR
            VERSION_VAR SCIP_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (SCIP_FOUND)
    set(SCIP_INCLUDE_DIRS "${SCIP_INCLUDE_DIR}")
    set(SCIP_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    # See: https://cmake.org/cmake/help/latest/module/CheckLibraryExists.html
    check_library_exists(m floor "" HAVE_LIBM)
    if (HAVE_LIBM)
        set(SCIP_LINK_LIBRARIES ${SCIP_LINK_LIBRARIES} m)
    endif ()

    if (UNIX)
        # Required under Unix since 12.8
        set(SCIP_LINK_LIBRARIES ${SCIP_LINK_LIBRARIES} dl)
    endif ()

    if (NOT TARGET SCIP::SCIP)
        add_library(SCIP::SCIP STATIC IMPORTED)
        set_target_properties(
                SCIP::SCIP PROPERTIES
                IMPORTED_LOCATION "${SCIP_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${SCIP_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${SCIP_LINK_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(SCIP_INCLUDE_DIR
                 SCIP_LIBRARY
                 SCIP_LIBRARY_DEBUG
                 SCIP_VERSION)

# --------------------------------------------------------------------------- #
