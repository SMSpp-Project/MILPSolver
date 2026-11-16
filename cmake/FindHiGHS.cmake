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
#        HiGHS_DLL           - The found runtime DLL, where there is one      #
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
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

find_package(ZLIB REQUIRED QUIET)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_HiGHS_cached HiGHS_INCLUDE_DIR HiGHS_LIBRARY HiGHS_LIBRARY_DEBUG
        HiGHS_DLL HiGHS_DLL_DEBUG HiGHS_EXTRAS_LIBRARY HiGHS_VERSION)
foreach (_var ${_HiGHS_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_HiGHS_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_HiGHS_stale)
    message(STATUS "Cached HiGHS paths not found, searching again")
    foreach (_var ${_HiGHS_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_HiGHS_cached)
unset(_HiGHS_stale)
unset(_var)
unset(_path)

# The cached version is the one of the include directory it was read from
if (NOT "${HiGHS_INCLUDE_DIR}" STREQUAL
        "${_HiGHS_VERSION_INCLUDE_DIR}")
    unset(HiGHS_VERSION CACHE)
endif ()

# Check if already in cache
if (HiGHS_INCLUDE_DIR AND HiGHS_LIBRARY AND HiGHS_LIBRARY_DEBUG
        AND HiGHS_VERSION)
    set(HiGHS_FOUND TRUE)
endif ()

if (NOT HiGHS_FOUND)

    # ----- Find the HiGHS include directory -------------------------------- #
    find_path(HiGHS_INCLUDE_DIR
            NAMES Highs.h interfaces/highs_c_api.h
            PATHS ${HiGHS_ROOT}
            PATH_SUFFIXES include/highs src
            DOC "HiGHS include directory.")

    # ----- Find the HiGHS library ------------------------------------------ #
    # also where a build of HiGHS, conda and vcpkg put it, with the debug
    # one apart where there is one, and the DLL where there is one: where
    # the library is a DLL, the one found is its import library
    set(_HiGHS_dirs ${HiGHS_ROOT}/lib ${HiGHS_ROOT}/build/lib/Release
            $ENV{LIBRARY_LIB})
    set(_HiGHS_debug_dirs ${HiGHS_ROOT}/debug/lib
            ${HiGHS_ROOT}/build/lib/Debug)
    set(_HiGHS_bins ${HiGHS_ROOT}/bin ${HiGHS_ROOT}/build/bin/Release
            $ENV{LIBRARY_BIN})
    set(_HiGHS_debug_bins ${HiGHS_ROOT}/debug/bin
            ${HiGHS_ROOT}/build/bin/Debug)
    if (VCPKG_INSTALLED_DIR)
        set(_HiGHS_vcpkg ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET})
        list(APPEND _HiGHS_dirs ${_HiGHS_vcpkg}/lib)
        list(APPEND _HiGHS_debug_dirs ${_HiGHS_vcpkg}/debug/lib)
        list(APPEND _HiGHS_bins ${_HiGHS_vcpkg}/bin)
        list(APPEND _HiGHS_debug_bins ${_HiGHS_vcpkg}/debug/bin)
        unset(_HiGHS_vcpkg)
    endif ()

    find_library(HiGHS_LIBRARY
            NAMES highs
            PATHS ${_HiGHS_dirs}
            DOC "HiGHS library.")

    find_library(HiGHS_LIBRARY_DEBUG
            NAMES highs
            PATHS ${_HiGHS_debug_dirs}
            NO_DEFAULT_PATH
            DOC "HiGHS debug library.")

    # Release-only distributions (e.g. conda-forge) ship no debug build:
    # fall back to the release library so a Release configure succeeds.
    if (NOT HiGHS_LIBRARY_DEBUG)
        set(HiGHS_LIBRARY_DEBUG ${HiGHS_LIBRARY}
                CACHE FILEPATH "HiGHS debug library." FORCE)
    endif ()

    set(_HiGHS_dll highs${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(_HiGHS_dlls ${_HiGHS_dll} ${CMAKE_SHARED_LIBRARY_PREFIX}${_HiGHS_dll})
    find_file(HiGHS_DLL
            NAMES ${_HiGHS_dlls}
            PATHS ${_HiGHS_bins}
            NO_DEFAULT_PATH
            DOC "HiGHS runtime DLL.")

    find_file(HiGHS_DLL_DEBUG
            NAMES ${_HiGHS_dlls}
            PATHS ${_HiGHS_debug_bins}
            NO_DEFAULT_PATH
            DOC "HiGHS debug runtime DLL.")

    if (NOT HiGHS_DLL_DEBUG AND HiGHS_DLL)
        set(HiGHS_DLL_DEBUG ${HiGHS_DLL}
                CACHE FILEPATH "HiGHS debug runtime DLL." FORCE)
    endif ()
    unset(_HiGHS_dirs)
    unset(_HiGHS_debug_dirs)
    unset(_HiGHS_bins)
    unset(_HiGHS_debug_bins)
    unset(_HiGHS_dll)
    unset(_HiGHS_dlls)

    # ----- Parse the version ----------------------------------------------- #
    if (HiGHS_INCLUDE_DIR)
        file(STRINGS
                "${HiGHS_INCLUDE_DIR}/HConfig.h"
                _HiGHS_version_lines
                REGEX "#define HIGHS_VERSION_(MAJOR|MINOR|PATCH)")

        string(REGEX REPLACE ".*HIGHS_VERSION_MAJOR *\([0-9]*\).*" "\\1"
                _HiGHS_version_major "${_HiGHS_version_lines}")
        string(REGEX REPLACE ".*HIGHS_VERSION_MINOR *\([0-9]*\).*" "\\1"
                _HiGHS_version_minor "${_HiGHS_version_lines}")
        string(REGEX REPLACE ".*HIGHS_VERSION_PATCH *\([0-9]*\).*" "\\1"
                _HiGHS_version_patch "${_HiGHS_version_lines}")

        string(CONCAT HiGHS_VERSION "${_HiGHS_version_major}."
                "${_HiGHS_version_minor}.${_HiGHS_version_patch}")
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
    # https://cmake.org/cmake/help/latest/module/
    # FindPackageHandleStandardArgs.html
    # The debug library and the DLLs are optional: the former falls back
    # to the release one, and a static build (e.g., a static vcpkg triplet)
    # has no DLL, nor has a platform whose shared libraries are not DLLs
    find_package_handle_standard_args(
            HiGHS
            REQUIRED_VARS HiGHS_LIBRARY HiGHS_INCLUDE_DIR
            VERSION_VAR HiGHS_VERSION)
endif ()

# A found package keeps its version in the cache, with the include
# directory it was read from
if (HiGHS_FOUND)
    set(HiGHS_VERSION "${HiGHS_VERSION}" CACHE INTERNAL
            "HiGHS version.")
    set(_HiGHS_VERSION_INCLUDE_DIR "${HiGHS_INCLUDE_DIR}"
            CACHE INTERNAL "Include directory of the HiGHS version.")
endif ()

# ----- Export the target --------------------------------------------------- #
if (HiGHS_FOUND)
    set(HiGHS_INCLUDE_DIRS ${HiGHS_INCLUDE_DIR})
    set(HiGHS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ZLIB::ZLIB ${CMAKE_DL_LIBS})

    # If HiGHS was built with the HiPO interior point solver (HIPO is then
    # #define-d in its HConfig.h) it depends on a BLAS, which must be linked
    # explicitly when libhighs is static; this mirrors the conditional
    # find_dependency(BLAS) in the upstream highs-config.cmake.
    file(STRINGS "${HiGHS_INCLUDE_DIR}/HConfig.h" _HiGHS_hipo_line
            REGEX "^#define HIPO[ \t]*$")
    if (_HiGHS_hipo_line)
        find_package(BLAS REQUIRED QUIET)
        set(HiGHS_LIBRARIES ${HiGHS_LIBRARIES} ${BLAS_LIBRARIES})
    endif ()
    unset(_HiGHS_hipo_line)

    # From version 1.15, a static libhighs calls the functions of the extras
    # library of HiGHS, which is static too in such a build (shared, and
    # loaded at run time, otherwise) and must be linked explicitly, as the
    # highs target of the upstream highs-config.cmake does; an import
    # library, which may have the suffix of a static one, comes with a DLL
    if (NOT HiGHS_DLL AND
            HiGHS_LIBRARY MATCHES "\\${CMAKE_STATIC_LIBRARY_SUFFIX}$")
        get_filename_component(_HiGHS_libdir "${HiGHS_LIBRARY}" DIRECTORY)
        set(_HiGHS_extras ${CMAKE_STATIC_LIBRARY_PREFIX}highs_extras)
        find_library(HiGHS_EXTRAS_LIBRARY
                NAMES ${_HiGHS_extras}${CMAKE_STATIC_LIBRARY_SUFFIX}
                PATHS ${_HiGHS_libdir}
                NO_DEFAULT_PATH
                DOC "HiGHS extras library.")
        if (HiGHS_EXTRAS_LIBRARY)
            set(HiGHS_LIBRARIES ${HiGHS_EXTRAS_LIBRARY} ${HiGHS_LIBRARIES})
        endif ()
        unset(_HiGHS_libdir)
        unset(_HiGHS_extras)
    endif ()

    # A DLL is imported with its import library, any other library as it is
    if (NOT TARGET HiGHS::HiGHS)
        if (HiGHS_DLL)
            add_library(HiGHS::HiGHS SHARED IMPORTED)
            set_target_properties(
                    HiGHS::HiGHS PROPERTIES
                    IMPORTED_IMPLIB "${HiGHS_LIBRARY}"
                    IMPORTED_IMPLIB_DEBUG "${HiGHS_LIBRARY_DEBUG}"
                    IMPORTED_LOCATION "${HiGHS_DLL}"
                    IMPORTED_LOCATION_DEBUG "${HiGHS_DLL_DEBUG}")
        else ()
            add_library(HiGHS::HiGHS UNKNOWN IMPORTED)
            set_target_properties(
                    HiGHS::HiGHS PROPERTIES
                    IMPORTED_LOCATION "${HiGHS_LIBRARY}"
                    IMPORTED_LOCATION_DEBUG "${HiGHS_LIBRARY_DEBUG}")
        endif ()
        set_target_properties(
                HiGHS::HiGHS PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${HiGHS_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${HiGHS_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(HiGHS_INCLUDE_DIR
        HiGHS_LIBRARY
        HiGHS_LIBRARY_DEBUG
        HiGHS_DLL
        HiGHS_DLL_DEBUG
        HiGHS_VERSION)

# --------------------------------------------------------------------------- #
