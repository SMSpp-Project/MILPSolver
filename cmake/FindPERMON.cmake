# --------------------------------------------------------------------------- #
#    CMake find module for PERMON                                             #
#                                                                             #
#    This module finds PERMON, and PETSc on which it is built, include        #
#    directories and libraries.                                               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(PERMON [REQUIRED])                                      #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        PERMON_FOUND        - True if headers are found                      #
#        PERMON_INCLUDE_DIRS - Include directories (of PERMON and PETSc)      #
#        PERMON_LIBRARIES    - Libraries to be linked (PERMON and PETSc)      #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        PERMON_ROOT         - Custom path to PERMON                          #
#        PETSC_DIR           - Custom path to PETSc (defaults to PERMON_ROOT) #
#                                                                             #
#    The following IMPORTED targets are also defined:                         #
#                                                                             #
#        PERMON::PETSc                                                        #
#        PERMON::PERMON                                                       #
#                                                                             #
#    This find module is provided because PERMON installs no CMake or         #
#    pkg-config file; PETSc is found by its pkg-config file, which has to be  #
#    3.22 or later (PETSC_CURRENT), and PERMON exposes no version number.     #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_PERMON_cached PERMON_INCLUDE_DIR PERMON_LIBRARY)
foreach (_var ${_PERMON_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_PERMON_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_PERMON_stale)
    message(STATUS "Cached PERMON paths not found, searching again")
    foreach (_var ${_PERMON_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_PERMON_cached)
unset(_PERMON_stale)
unset(_var)
unset(_path)

# ----- Find PETSc ---------------------------------------------------------- #
if (NOT PETSC_DIR)
    if (DEFINED ENV{PETSC_DIR})
        set(PETSC_DIR $ENV{PETSC_DIR})
    else ()
        set(PETSC_DIR ${PERMON_ROOT})
    endif ()
endif ()
find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
    set(_PERMON_pkg_path $ENV{PKG_CONFIG_PATH})
    if (PETSC_DIR)
        set(ENV{PKG_CONFIG_PATH}
            "${PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif ()
    pkg_check_modules(PERMON_PETSC QUIET IMPORTED_TARGET PETSc>=3.22)
    set(ENV{PKG_CONFIG_PATH} "${_PERMON_pkg_path}")
    unset(_PERMON_pkg_path)
endif ()

# ----- Find the PERMON include directory ----------------------------------- #
# The private headers (permon/private/qpimpl.h) are used as well
find_path(PERMON_INCLUDE_DIR
        NAMES permonqps.h permon/private/qpimpl.h
        PATHS ${PERMON_ROOT} $ENV{PERMON_DIR}
        PATH_SUFFIXES include
        DOC "PERMON include directory.")

# ----- Find the PERMON library --------------------------------------------- #
find_library(PERMON_LIBRARY
        NAMES permon
        PATHS ${PERMON_ROOT} $ENV{PERMON_DIR}
        PATH_SUFFIXES lib
        DOC "PERMON library.")

# ----- Handle the standard arguments --------------------------------------- #
# The following macro manages the QUIET and REQUIRED options passed to
# find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
# set. REQUIRED_VARS should be cache entries and not output variables.
# See:
# https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
find_package_handle_standard_args(
        PERMON
        REQUIRED_VARS PERMON_LIBRARY PERMON_INCLUDE_DIR PERMON_PETSC_FOUND)

# ----- Export the targets -------------------------------------------------- #
if (PERMON_FOUND)
    set(PERMON_INCLUDE_DIRS ${PERMON_INCLUDE_DIR}
        ${PERMON_PETSC_INCLUDE_DIRS})
    set(PERMON_LIBRARIES ${PERMON_LIBRARY} ${PERMON_PETSC_LINK_LIBRARIES})

    if (NOT TARGET PERMON::PETSc)
        add_library(PERMON::PETSc INTERFACE IMPORTED)
        target_link_libraries(PERMON::PETSc INTERFACE PkgConfig::PERMON_PETSC)
    endif ()

    if (NOT TARGET PERMON::PERMON)
        add_library(PERMON::PERMON UNKNOWN IMPORTED)
        set_target_properties(
                PERMON::PERMON PROPERTIES
                IMPORTED_LOCATION "${PERMON_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${PERMON_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES PERMON::PETSc)
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(PERMON_INCLUDE_DIR PERMON_LIBRARY)

# --------------------------------------------------------------------------- #
