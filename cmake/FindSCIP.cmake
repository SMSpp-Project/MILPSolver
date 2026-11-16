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
#        SCIP_DLL           - The found runtime DLL, where there is one       #
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

# ----- Find the path to SCIP ----------------------------------------------- #
# This takes the greatest scip* found where the installers put it (/opt,
# /Library, C:/Program Files), of which only those of the system at hand
# are there

foreach (dir /opt /Library "C:/Program Files" "C:/Program Files (x86)")
    file(GLOB SCIP_DIRS "${dir}/scip*")
    if (NOT IS_DIRECTORY "${SCIP_ROOT}")
        if (NOT "${SCIP_ROOT}" STREQUAL "")
            message(STATUS "Specified SCIP: ${SCIP_ROOT} not found")
        endif ()
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

# A static SCIP built with PaPILO needs TBB, which is linked where found
find_package(TBB QUIET)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_SCIP_cached SCIP_INCLUDE_DIR SCIP_LIBRARY SCIP_DLL SCIP_VERSION)
foreach (_var ${_SCIP_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_SCIP_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_SCIP_stale)
    message(STATUS "Cached SCIP paths not found, searching again")
    foreach (_var ${_SCIP_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_SCIP_cached)
unset(_SCIP_stale)
unset(_var)
unset(_path)

# The cached version is the one of the include directory it was read from
if (NOT "${SCIP_INCLUDE_DIR}" STREQUAL
        "${_SCIP_VERSION_INCLUDE_DIR}")
    unset(SCIP_VERSION CACHE)
endif ()

# Check if already in cache
if (SCIP_INCLUDE_DIR AND SCIP_LIBRARY AND SCIP_VERSION)
    set(SCIP_FOUND TRUE)
endif ()

if (NOT SCIP_FOUND)

    # ----- Find the SCIP include directory --------------------------------- #
    find_path(SCIP_INCLUDE_DIR
            NAMES scip/scip.h
            PATHS ${SCIP_ROOT}/include
            DOC "SCIP include directory.")

    # ----- Find the SCIP library ------------------------------------------- #
    # where the library is a DLL, the one found is its import library, and
    # the DLL is in bin/
    find_library(SCIP_LIBRARY
            NAMES scip
            PATHS ${SCIP_ROOT}
            PATH_SUFFIXES lib
            DOC "SCIP library.")

    set(_SCIP_dll scip${CMAKE_SHARED_LIBRARY_SUFFIX})
    find_file(SCIP_DLL
            NAMES lib${_SCIP_dll} ${_SCIP_dll}
            PATHS ${SCIP_ROOT}
            PATH_SUFFIXES bin
            NO_DEFAULT_PATH
            DOC "SCIP runtime DLL.")
    unset(_SCIP_dll)

    # ----- Parse the version ----------------------------------------------- #
    if (SCIP_INCLUDE_DIR)
        file(STRINGS
                "${SCIP_INCLUDE_DIR}/scip/config.h"
                _SCIP_version_lines
                REGEX "#define SCIP_VERSION_(MAJOR|MINOR|PATCH)")

        string(REGEX REPLACE ".*SCIP_VERSION_MAJOR *\([0-9]*\).*" "\\1"
                _SCIP_version_major "${_SCIP_version_lines}")
        string(REGEX REPLACE ".*SCIP_VERSION_MINOR *\([0-9]*\).*" "\\1"
                _SCIP_version_minor "${_SCIP_version_lines}")
        string(REGEX REPLACE ".*SCIP_VERSION_PATCH *\([0-9]*\).*" "\\1"
                _SCIP_version_patch "${_SCIP_version_lines}")

        string(CONCAT SCIP_VERSION "${_SCIP_version_major}."
                "${_SCIP_version_minor}.${_SCIP_version_patch}")
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
    # https://cmake.org/cmake/help/latest/module/
    # FindPackageHandleStandardArgs.html
    # The DLL is not required: a static build has none, nor has a platform
    # whose shared libraries are not DLLs
    find_package_handle_standard_args(
            SCIP
            REQUIRED_VARS SCIP_LIBRARY SCIP_INCLUDE_DIR
            VERSION_VAR SCIP_VERSION)
endif ()

# A found package keeps its version in the cache, with the include
# directory it was read from
if (SCIP_FOUND)
    set(SCIP_VERSION "${SCIP_VERSION}" CACHE INTERNAL
            "SCIP version.")
    set(_SCIP_VERSION_INCLUDE_DIR "${SCIP_INCLUDE_DIR}"
            CACHE INTERNAL "Include directory of the SCIP version.")
endif ()

# ----- Export the target --------------------------------------------------- #
if (SCIP_FOUND)
    set(SCIP_INCLUDE_DIRS ${SCIP_INCLUDE_DIR})
    set(SCIP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    if (TARGET TBB::tbb)
        set(SCIP_LIBRARIES ${SCIP_LIBRARIES} TBB::tbb)
    endif ()

    set(SCIP_LIBRARIES ${SCIP_LIBRARIES} ${CMAKE_DL_LIBS})

    # A DLL is imported with its import library, any other library as it is
    if (NOT TARGET SCIP::SCIP)
        if (SCIP_DLL)
            add_library(SCIP::SCIP SHARED IMPORTED)
            set_target_properties(
                    SCIP::SCIP PROPERTIES
                    IMPORTED_IMPLIB "${SCIP_LIBRARY}"
                    IMPORTED_LOCATION "${SCIP_DLL}")
        else ()
            add_library(SCIP::SCIP UNKNOWN IMPORTED)
            set_target_properties(
                    SCIP::SCIP PROPERTIES
                    IMPORTED_LOCATION "${SCIP_LIBRARY}")
        endif ()
        set_target_properties(
                SCIP::SCIP PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${SCIP_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${SCIP_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(SCIP_INCLUDE_DIR
        SCIP_LIBRARY
        SCIP_DLL
        SCIP_VERSION)

# --------------------------------------------------------------------------- #
