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
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
find_package(Threads QUIET)

find_package(ZLIB REQUIRED QUIET)

# Check if already in cache
if (HiGHS_INCLUDE_DIR AND HiGHS_LIBRARY AND HiGHS_LIBRARY_DEBUG AND HIGHS_VERSION)
    set(HiGHS_FOUND TRUE)
else ()

    # ----- Find the HiGHS include directory -------------------------------- #
    find_path(HiGHS_INCLUDE_DIR
              NAMES Highs.h interfaces/highs_c_api.h
              PATHS ${HiGHS_ROOT}
              PATH_SUFFIXES include/highs src
              DOC "HiGHS include directory.")

    # ----- Find the HiGHS library ------------------------------------------ #
    if (UNIX)
        find_library(HiGHS_LIBRARY
                NAMES highs
                PATHS ${HiGHS_ROOT}/lib
                DOC "HiGHS library.")

        set(HiGHS_LIBRARY_DEBUG ${HiGHS_LIBRARY}
                CACHE FILEPATH "HiGHS debug library." FORCE)
    elseif (WIN32)
        find_library(HiGHS_LIBRARY
                NAMES highs
                PATHS ${HiGHS_ROOT}/lib
                      ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
                      $ENV{LIBRARY_LIB}
                NO_DEFAULT_PATH
                DOC "HiGHS library.")

        find_library(HiGHS_LIBRARY_DEBUG
                     NAMES highs
                     PATHS ${HiGHS_ROOT}/debug/lib
                           ${HiGHS_ROOT}/build/lib/Debug
                           ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/debug/lib
                     NO_DEFAULT_PATH
                     DOC "HiGHS debug library.")
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (HiGHS_INCLUDE_DIR)
        file(STRINGS
                "${HiGHS_INCLUDE_DIR}/HConfig.h"
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
            REQUIRED_VARS HiGHS_LIBRARY HiGHS_INCLUDE_DIR
            VERSION_VAR HiGHS_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (HiGHS_FOUND)
    set(HiGHS_INCLUDE_DIRS ${HiGHS_INCLUDE_DIR})
    set(HiGHS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ZLIB::ZLIB)

    if (NOT TARGET HiGHS::HiGHS)
        add_library(HiGHS::HiGHS UNKNOWN IMPORTED)
        set_target_properties(
                HiGHS::HiGHS PROPERTIES
                IMPORTED_LOCATION "${HiGHS_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${HiGHS_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${HiGHS_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${HiGHS_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(HiGHS_INCLUDE_DIR
                 HiGHS_LIBRARY
                 HiGHS_LIBRARY_DEBUG
                 HiGHS_VERSION)

# --------------------------------------------------------------------------- #
