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
#        MOSEK_DLL           - The found runtime DLL, where there is one      #
#        MOSEK_RUNTIME_DIR   - The directory of the library, where it is not  #
#                              a DLL, for the run path of an installed target #
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

# ----- Find the MOSEK directories ------------------------------------------ #
# Where the installers put it (/opt, /Library, C:/Program Files) or the home
# directory, of which only those of the system at hand are there

# (in the form of CMake, with no backslash, whatever the environment gives)

file(TO_CMAKE_PATH "$ENV{HOME}" _mosek_home)
file(TO_CMAKE_PATH "$ENV{USERPROFILE}" _mosek_profile)
set(MOSEK_DIRS /opt/mosek /Library/mosek "C:/Program Files/Mosek"
        "${_mosek_home}/mosek" "${_mosek_profile}/mosek")
unset(_mosek_home)
unset(_mosek_profile)

# ----- Find the path to MOSEK ---------------------------------------------- #
# The directory of a version is the one holding tools/platform: MOSEK_ROOT
# itself, or the last of its subdirectories that holds one.

function(find_mosek_version_dir root var)
    file(TO_CMAKE_PATH "${root}" root)
    if (IS_DIRECTORY "${root}/tools/platform")
        set(${var} "${root}" PARENT_SCOPE)
    else ()
        file(GLOB _versions "${root}/*/tools/platform")
        if (_versions)
            list(SORT _versions)
            list(REVERSE _versions)
            list(GET _versions 0 _version)
            get_filename_component(_version "${_version}" DIRECTORY)
            get_filename_component(_version "${_version}" DIRECTORY)
            set(${var} "${_version}" PARENT_SCOPE)
        else ()
            set(${var} "" PARENT_SCOPE)
        endif ()
    endif ()
endfunction()

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

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_mosek_cached MOSEK_INCLUDE_DIR MOSEK_LIBRARY MOSEK_DLL MOSEK_VERSION)
foreach (_var ${_mosek_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_mosek_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_mosek_stale)
    message(STATUS "Cached MOSEK paths not found, searching again")
    foreach (_var ${_mosek_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_mosek_cached)
unset(_mosek_stale)
unset(_var)
unset(_path)

# The cached version is the one of the include directory it was read from
if (NOT "${MOSEK_INCLUDE_DIR}" STREQUAL
        "${_MOSEK_VERSION_INCLUDE_DIR}")
    unset(MOSEK_VERSION CACHE)
endif ()

# Check if already in cache
if (MOSEK_INCLUDE_DIR AND MOSEK_LIBRARY AND MOSEK_VERSION)
    set(MOSEK_FOUND TRUE)
endif ()

if (NOT MOSEK_FOUND)

    # ----- Find the MOSEK include directory -------------------------------- #
    # The directory of the platform is the one of those MOSEK is built for
    # that is there, and the library is taken from it, hence so is the
    # include directory
    unset(MOSEK_INCLUDE_DIR CACHE)
    find_path(MOSEK_INCLUDE_DIR
            NAMES mosek.h
            PATHS ${MOSEK_VERSION_DIR}/tools/platform
            PATH_SUFFIXES linux64x86/h linuxaarch64/h osx64x86/h
                          osxaarch64/h win64x86/h
            NO_DEFAULT_PATH
            DOC "MOSEK include directory.")
    get_filename_component(MOSEK_DIR "${MOSEK_INCLUDE_DIR}" DIRECTORY)

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
        set(_mosek_name
                mosek64_${_mosek_version_major}_${_mosek_version_minor})
        unset(_mosek_version_lines)
        unset(_mosek_version_major)
        unset(_mosek_version_minor)
        unset(_mosek_version_revision)
    endif ()

    # ----- Find the MOSEK library ------------------------------------------ #
    # The library is mosek64, to which CMake adds the prefix and the suffix
    # of the platform, or, where it is a DLL, the import library of the DLL,
    # whose name carries the major and minor version
    unset(MOSEK_LIBRARY CACHE)
    unset(MOSEK_DLL CACHE)
    if (_mosek_name)
        find_library(MOSEK_LIBRARY
                NAMES mosek64 ${_mosek_name}
                PATHS ${MOSEK_DIR}/bin
                NO_DEFAULT_PATH
                DOC "MOSEK library.")
        set(_mosek_dll ${_mosek_name}${CMAKE_SHARED_LIBRARY_SUFFIX})
        find_file(MOSEK_DLL
                NAMES ${_mosek_dll} ${CMAKE_SHARED_LIBRARY_PREFIX}${_mosek_dll}
                PATHS ${MOSEK_DIR}/bin
                NO_DEFAULT_PATH
                DOC "MOSEK runtime DLL.")
        unset(_mosek_dll)
        unset(_mosek_name)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/
    # FindPackageHandleStandardArgs.html
    # The DLL is not required, since only some platforms have one
    find_package_handle_standard_args(
            MOSEK
            REQUIRED_VARS MOSEK_LIBRARY MOSEK_INCLUDE_DIR
            VERSION_VAR MOSEK_VERSION)

    # ----- Find the MOSEK license ------------------------------------------ #
    # MOSEK looks for it in $MOSEKLM_LICENSE_FILE, then in mosek/mosek.lic
    # under the home directory of the user.
    set(MOSEK_LICENSE_FOUND FALSE)

    if (DEFINED ENV{MOSEKLM_LICENSE_FILE})
        set(MOSEK_LICENSE_FOUND TRUE)
    else ()
        foreach (_mosek_home "$ENV{HOME}" "$ENV{USERPROFILE}")
            if (NOT "${_mosek_home}" STREQUAL "" AND
                    EXISTS "${_mosek_home}/mosek/mosek.lic")
                set(MOSEK_LICENSE_FOUND TRUE)
                message(STATUS "MOSEK license file found at: "
                        "${_mosek_home}/mosek/mosek.lic")
                break()
            endif ()
        endforeach ()
        unset(_mosek_home)
    endif ()

    if (MOSEK_FOUND AND NOT MOSEK_LICENSE_FOUND)
        message(WARNING "MOSEK license file not found: MOSEKMILPSolver is\
                         built, but it cannot solve anything until the file\
                         mosek.lic is put in the mosek directory of your home\
                         or MOSEKLM_LICENSE_FILE points to it.")
    endif ()
endif ()

# A found package keeps its version in the cache, with the include
# directory it was read from
if (MOSEK_FOUND)
    set(MOSEK_VERSION "${MOSEK_VERSION}" CACHE INTERNAL
            "MOSEK version.")
    set(_MOSEK_VERSION_INCLUDE_DIR "${MOSEK_INCLUDE_DIR}"
            CACHE INTERNAL "Include directory of the MOSEK version.")
endif ()

# ----- Export the target --------------------------------------------------- #
if (MOSEK_FOUND)
    set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
    set(MOSEK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

    # A DLL is imported with its import library, any other library as it is
    if (NOT TARGET MOSEK::MOSEK)
        if (MOSEK_DLL)
            add_library(MOSEK::MOSEK SHARED IMPORTED)
            set_target_properties(
                    MOSEK::MOSEK PROPERTIES
                    IMPORTED_IMPLIB "${MOSEK_LIBRARY}"
                    IMPORTED_LOCATION "${MOSEK_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${MOSEK_LIBRARIES}")
        else ()
            # MOSEK keeps next to its library its own libtbb, older than the
            # one of the system that other packages (e.g., SCIP) link: with
            # the directory of MOSEK in the run path of a program, its libtbb
            # would be the one loaded for all of them, which CMake warns of
            # for every target. The library is taken from a directory of the
            # build tree that holds only links to it, and libmosek64 finds a
            # libtbb by itself, through its own RPATH ($ORIGIN), if none is
            # loaded already [see MOSEK_RUNTIME_DIR for the installation]
            get_filename_component(_mosek_real "${MOSEK_LIBRARY}" REALPATH)
            get_filename_component(MOSEK_RUNTIME_DIR "${_mosek_real}"
                    DIRECTORY)
            set(_mosek_links "${CMAKE_BINARY_DIR}/mosek-lib")
            file(MAKE_DIRECTORY "${_mosek_links}")
            file(GLOB _mosek_libs "${MOSEK_RUNTIME_DIR}/libmosek64*")
            foreach (_mosek_lib ${_mosek_libs})
                get_filename_component(_mosek_name "${_mosek_lib}" NAME)
                file(CREATE_LINK "${_mosek_lib}"
                        "${_mosek_links}/${_mosek_name}" SYMBOLIC)
            endforeach ()
            get_filename_component(_mosek_name "${MOSEK_LIBRARY}" NAME)
            add_library(MOSEK::MOSEK UNKNOWN IMPORTED)
            set_target_properties(
                    MOSEK::MOSEK PROPERTIES
                    IMPORTED_LOCATION "${_mosek_links}/${_mosek_name}"
                    INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIRS}"
                    INTERFACE_LINK_LIBRARIES "${MOSEK_LIBRARIES}")
            unset(_mosek_real)
            unset(_mosek_links)
            unset(_mosek_libs)
            unset(_mosek_lib)
            unset(_mosek_name)
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(MOSEK_INCLUDE_DIR
        MOSEK_LIBRARY
        MOSEK_DLL
        MOSEK_VERSION)

# --------------------------------------------------------------------------- #
