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
#        GUROBI_DIR    - Custom path to GUROBI                                #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        GUROBI::Gurobi                                                       #
#                                                                             #
#    This find module is provided because GUROBI does not provide             #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Enrico Calandrini                              #
#                         Dipartimento di Matematica                          #
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
    set(GUROBI_LIB_PATH_SUFFIXES lib)
endif ()

# ----- Find the path to GUROBI --------------------------------------------- #

if (NOT GUROBI_DIR)
    foreach (dir ${GUROBI_DIRS})
        file(GLOB GUROBI_DIRS "${dir}/gurobi*")
        list(SORT GUROBI_DIRS)
        list(REVERSE GUROBI_DIRS)
        if (GUROBI_DIRS)
            list(GET GUROBI_DIRS 0 GUROBI_DIR_)
            message(STATUS "Found Gurobi: ${GUROBI_DIR_}")
            break()
        endif ()
    endforeach ()

    if (NOT GUROBI_DIR_)
        set(GUROBI_DIR_ GUROBI_DIR-NOTFOUND)
    endif ()
    # Set the path in the cache
    set(GUROBI_DIR ${GUROBI_DIR_})
endif ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

# Check if already in cache
if (GUROBI_INCLUDE_DIR AND GUROBI_LIBRARY AND GUROBI_LIBRARY_DEBUG)
    set(GUROBI_FOUND TRUE)
else ()

    # ----- Find the GUROBI include directory ------------------------------- #
    set(GUROBI_DIR ${GUROBI_DIR}/linux64)
    # Note that find_path() creates a cache entry
    find_path(GUROBI_INCLUDE_DIR gurobi_c.h
              PATHS ${GUROBI_DIR}/include
              DOC "GUROBI include directory.")

    # ----- Find the GUROBI library ----------------------------------------- #
    # Note that find_library() creates a cache entry
    find_library(GUROBI_LIBRARY
                 NAMES gurobi gurobi100 gurobi1002
                 PATHS ${GUROBI_DIR}
                 PATH_SUFFIXES ${GUROBI_LIB_PATH_SUFFIXES}
                 DOC "GUROBI library.")
    set(GUROBI_LIBRARY_DEBUG ${GUROBI_LIBRARY} CACHE FILEPATH "Debug GUROBI library.")

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
        unset(_gurobi_version_patch)
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
        # Required under Unix since 12.8
        set(GUROBI_LINK_LIBRARIES ${GUROBI_LINK_LIBRARIES} dl)
    endif ()

    if (NOT TARGET GUROBI::Gurobi)
        add_library(GUROBI::Gurobi STATIC IMPORTED)
        set_target_properties(
                GUROBI::Gurobi PROPERTIES
                IMPORTED_LOCATION "${GUROBI_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${GUROBI_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
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
