# --------------------------------------------------------------------------- #
#    CMake find module for MOSEK                                              #
#                                                                             #
#    This module finds MOSEK include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(MOSEK [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        MOSEK_FOUND         - True if headers are found                      #
#        MOSEK_LICENSE_FOUND - True if the mosek.lic file is found            #
#        MOSEK_INCLUDE_DIRS  - Include directories                            #
#        MOSEK_LIBRARIES     - Libraries to be linked                         #
#        MOSEK_DLL           - The found runtime DLL (Windows only)           #
#        MOSEK_VERSION       - Version number                                 #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        MOSEK_ROOT          - Custom path to MOSEK, either the directory of  #
#                              a version (the one holding tools/platform) or  #
#                              the one above it, holding a directory per      #
#                              version, of which the last one is taken        #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        MOSEK::MOSEK                                                         #
#                                                                             #
#    This find module is provided because MOSEK does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find the MOSEK directory and the platform --------------------------- #
# Based on the OS generate:
# - a list of possible MOSEK directories
# - the name of the platform directory of MOSEK

if (UNIX)
    if (APPLE)
        # macOS (usually /Library or the home directory)
        set(MOSEK_DIRS /Library/mosek "$ENV{HOME}/mosek")
        if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
            set(MOSEK_PLATFORM osxaarch64)
        else ()
            set(MOSEK_PLATFORM osx64x86)
        endif ()
    else ()
        # Other Unix-based systems (usually /opt or the home directory)
        set(MOSEK_DIRS /opt/mosek "$ENV{HOME}/mosek")
        if (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
            set(MOSEK_PLATFORM linuxaarch64)
        else ()
            set(MOSEK_PLATFORM linux64x86)
        endif ()
    endif ()
elseif (WIN32)
    # Windows (usually C:/Program Files/Mosek)
    set(MOSEK_DIRS "C:/Program Files/Mosek")
    set(MOSEK_PLATFORM win64x86)
endif ()

# ----- Find the path to MOSEK ---------------------------------------------- #
# The directory of a version is the one holding tools/platform: MOSEK_ROOT
# itself, or the last of its subdirectories that holds one.

macro(find_mosek_version_dir root var)
    if (IS_DIRECTORY "${root}/tools/platform")
        set(${var} "${root}")
    else ()
        file(GLOB _mosek_versions "${root}/*/tools/platform")
        if (_mosek_versions)
            list(SORT _mosek_versions)
            list(REVERSE _mosek_versions)
            list(GET _mosek_versions 0 _mosek_version)
            get_filename_component(_mosek_version "${_mosek_version}"
                    DIRECTORY)
            get_filename_component(${var} "${_mosek_version}" DIRECTORY)
        else ()
            set(${var} "")
        endif ()
        unset(_mosek_versions)
        unset(_mosek_version)
    endif ()
endmacro()

set(MOSEK_VERSION_DIR "")
if (NOT "${MOSEK_ROOT}" STREQUAL "")
    find_mosek_version_dir("${MOSEK_ROOT}" MOSEK_VERSION_DIR)
    if ("${MOSEK_VERSION_DIR}" STREQUAL "")
        message(STATUS "Specified MOSEK: ${MOSEK_ROOT} not found")
    endif ()
endif ()

if ("${MOSEK_VERSION_DIR}" STREQUAL "")
    foreach (dir ${MOSEK_DIRS})
        find_mosek_version_dir("${dir}" MOSEK_VERSION_DIR)
        if (NOT "${MOSEK_VERSION_DIR}" STREQUAL "")
            break()
        endif ()
    endforeach ()
endif ()

if (NOT "${MOSEK_VERSION_DIR}" STREQUAL "")
    message(STATUS "Using MOSEK: ${MOSEK_VERSION_DIR}")
endif ()

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

# Check if already in cache
if (WIN32)
    if (MOSEK_INCLUDE_DIR AND MOSEK_LIBRARY AND MOSEK_DLL AND MOSEK_VERSION)
        set(MOSEK_FOUND TRUE)
    endif ()
else ()
    if (MOSEK_INCLUDE_DIR AND MOSEK_LIBRARY AND MOSEK_VERSION)
        set(MOSEK_FOUND TRUE)
    endif ()
endif ()

if (NOT MOSEK_FOUND)

    set(MOSEK_DIR ${MOSEK_VERSION_DIR}/tools/platform/${MOSEK_PLATFORM})

    # ----- Find the MOSEK include directory -------------------------------- #
    find_path(MOSEK_INCLUDE_DIR
            NAMES mosek.h
            PATHS ${MOSEK_DIR}/h
            NO_DEFAULT_PATH
            DOC "MOSEK include directory.")

    # ----- Find the MOSEK library ------------------------------------------ #
    if (WIN32)
        file(GLOB _mosek_libs "${MOSEK_DIR}/bin/mosek64_*.lib")
    else ()
        file(GLOB _mosek_libs
                "${MOSEK_DIR}/bin/libmosek64${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif ()

    if (_mosek_libs)
        list(GET _mosek_libs 0 _mosek_lib)
        set(MOSEK_LIBRARY ${_mosek_lib}
                CACHE FILEPATH "MOSEK library." FORCE)
    else ()
        set(MOSEK_LIBRARY MOSEK_LIBRARY-NOTFOUND
                CACHE FILEPATH "MOSEK library." FORCE)
    endif ()
    unset(_mosek_libs)
    unset(_mosek_lib)

    if (WIN32)
        file(GLOB _mosek_dlls "${MOSEK_DIR}/bin/mosek64_*.dll")
        if (_mosek_dlls)
            list(GET _mosek_dlls 0 _mosek_dll)
            set(MOSEK_DLL ${_mosek_dll}
                    CACHE FILEPATH "MOSEK runtime DLL." FORCE)
        else ()
            set(MOSEK_DLL MOSEK_DLL-NOTFOUND
                    CACHE FILEPATH "MOSEK runtime DLL." FORCE)
        endif ()
        unset(_mosek_dlls)
        unset(_mosek_dll)
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (MOSEK_INCLUDE_DIR)
        file(STRINGS
                "${MOSEK_INCLUDE_DIR}/mosek.h"
                _mosek_version_lines
                REGEX "#define MSK_VERSION_(MAJOR|MINOR|REVISION)")

        string(REGEX REPLACE ".*MSK_VERSION_MAJOR *\([0-9]*\).*" "\\1"
                _mosek_version_major "${_mosek_version_lines}")
        string(REGEX REPLACE ".*MSK_VERSION_MINOR *\([0-9]*\).*" "\\1"
                _mosek_version_minor "${_mosek_version_lines}")
        string(REGEX REPLACE ".*MSK_VERSION_REVISION *\([0-9]*\).*" "\\1"
                _mosek_version_revision "${_mosek_version_lines}")

        set(MOSEK_VERSION "${_mosek_version_major}.${_mosek_version_minor}")
        set(MOSEK_VERSION "${MOSEK_VERSION}.${_mosek_version_revision}")
        unset(_mosek_version_lines)
        unset(_mosek_version_major)
        unset(_mosek_version_minor)
        unset(_mosek_version_revision)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    if (WIN32)
        find_package_handle_standard_args(
                MOSEK
                REQUIRED_VARS MOSEK_LIBRARY MOSEK_DLL MOSEK_INCLUDE_DIR
                VERSION_VAR MOSEK_VERSION)
    else ()
        find_package_handle_standard_args(
                MOSEK
                REQUIRED_VARS MOSEK_LIBRARY MOSEK_INCLUDE_DIR
                VERSION_VAR MOSEK_VERSION)
    endif ()

    # ----- Find the MOSEK license ------------------------------------------ #
    # MOSEK looks for it in $MOSEKLM_LICENSE_FILE, then in mosek/mosek.lic
    # under the home directory of the user.
    set(MOSEK_LICENSE_FOUND FALSE)

    if (DEFINED ENV{MOSEKLM_LICENSE_FILE})
        set(MOSEK_LICENSE_FOUND TRUE)
    else ()
        if (WIN32)
            set(_mosek_home "$ENV{USERPROFILE}")
        else ()
            set(_mosek_home "$ENV{HOME}")
        endif ()
        if (EXISTS "${_mosek_home}/mosek/mosek.lic")
            set(MOSEK_LICENSE_FOUND TRUE)
            message(STATUS "MOSEK license file found at: "
                    "${_mosek_home}/mosek/mosek.lic")
        endif ()
        unset(_mosek_home)
    endif ()

    if (MOSEK_FOUND AND NOT MOSEK_LICENSE_FOUND)
        message(WARNING "MOSEK license file not found: MOSEKMILPSolver is\
                         built, but it cannot solve anything until the file\
                         mosek.lic is put in the mosek directory of your home\
                         or MOSEKLM_LICENSE_FILE points to it.")
    endif ()
endif ()

# ----- Export the target --------------------------------------------------- #
if (MOSEK_FOUND)
    set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
    set(MOSEK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    if (NOT TARGET MOSEK::MOSEK)
        if (WIN32)
            add_library(MOSEK::MOSEK SHARED IMPORTED)
            set_target_properties(
                    MOSEK::MOSEK PROPERTIES
                    IMPORTED_IMPLIB "${MOSEK_LIBRARY}"
                    IMPORTED_LOCATION "${MOSEK_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${MOSEK_LIBRARIES}")
        else ()
            add_library(MOSEK::MOSEK UNKNOWN IMPORTED)
            set_target_properties(
                    MOSEK::MOSEK PROPERTIES
                    IMPORTED_LOCATION "${MOSEK_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${MOSEK_LIBRARIES}")
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
if (WIN32)
    mark_as_advanced(MOSEK_INCLUDE_DIR
            MOSEK_LIBRARY
            MOSEK_DLL
            MOSEK_VERSION)
else ()
    mark_as_advanced(MOSEK_INCLUDE_DIR
            MOSEK_LIBRARY
            MOSEK_VERSION)
endif ()

# --------------------------------------------------------------------------- #
