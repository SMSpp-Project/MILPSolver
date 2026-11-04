# --------------------------------------------------------------------------- #
#    CMake find module for HPIPM                                              #
#                                                                             #
#    This module finds HPIPM, and BLASFEO on which it is built, include       #
#    directories and libraries. Use it by invoking find_package() with the    #
#    form:                                                                    #
#                                                                             #
#        find_package(HPIPM [REQUIRED])                                       #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        HPIPM_FOUND         - True if headers are found                      #
#        HPIPM_INCLUDE_DIRS  - Include directories (of HPIPM and BLASFEO)     #
#        HPIPM_LIBRARIES     - Libraries to be linked (HPIPM and BLASFEO)     #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        HPIPM_ROOT          - Custom path to HPIPM                           #
#        BLASFEO_ROOT        - Custom path to BLASFEO                         #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        HPIPM::HPIPM                                                         #
#                                                                             #
#    This find module is provided because the CMake configuration file of     #
#    HPIPM is installed in <prefix>/cmake, where find_package() does not look #
#    on Linux and macOS, and does not bring BLASFEO along; HPIPM exposes no   #
#    version number in its headers: 0.1.2 or later is told apart by the masks #
#    of the bounds, which HPIPMMILPSolver uses.                               #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_HPIPM_cached HPIPM_INCLUDE_DIR HPIPM_LIBRARY BLASFEO_INCLUDE_DIR
        BLASFEO_LIBRARY)
foreach (_var ${_HPIPM_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_HPIPM_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_HPIPM_stale)
    message(STATUS "Cached HPIPM paths not found, searching again")
    foreach (_var ${_HPIPM_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_HPIPM_cached)
unset(_HPIPM_stale)
unset(_var)
unset(_path)

# Check if already in cache
if (HPIPM_INCLUDE_DIR AND HPIPM_LIBRARY AND BLASFEO_INCLUDE_DIR AND
    BLASFEO_LIBRARY)
    set(HPIPM_FOUND TRUE)
endif ()

if (NOT HPIPM_FOUND)

    # ----- Find the include directories ------------------------------------ #
    find_path(HPIPM_INCLUDE_DIR
            NAMES hpipm_d_dense_qp_ipm.h
            PATHS ${HPIPM_ROOT}
            PATH_SUFFIXES include include/hpipm
            DOC "HPIPM include directory.")

    find_path(BLASFEO_INCLUDE_DIR
            NAMES blasfeo_target.h
            PATHS ${BLASFEO_ROOT} ${HPIPM_ROOT}
            PATH_SUFFIXES include include/blasfeo
            DOC "BLASFEO include directory.")

    # ----- Check that it is 0.1.2 or later --------------------------------- #
    # The headers carry no version number: the masks of the bounds, which
    # HPIPMMILPSolver uses for the infinite ones, are in those of 0.1.2 and
    # later only, hence HPIPM_BOUND_MASKS is required
    if (HPIPM_INCLUDE_DIR AND EXISTS "${HPIPM_INCLUDE_DIR}/hpipm_d_dense_qp.h")
        file(STRINGS "${HPIPM_INCLUDE_DIR}/hpipm_d_dense_qp.h"
             HPIPM_BOUND_MASKS REGEX "d_dense_qp_set_lb_mask")
    endif ()

    # ----- Find the libraries ---------------------------------------------- #
    find_library(HPIPM_LIBRARY
            NAMES hpipm
            PATHS ${HPIPM_ROOT}
            PATH_SUFFIXES lib
            DOC "HPIPM library.")

    find_library(BLASFEO_LIBRARY
            NAMES blasfeo
            PATHS ${BLASFEO_ROOT} ${HPIPM_ROOT}
            PATH_SUFFIXES lib
            DOC "BLASFEO library.")

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
    # set. REQUIRED_VARS should be cache entries and not output variables.
    # See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            HPIPM
            REQUIRED_VARS HPIPM_LIBRARY HPIPM_INCLUDE_DIR BLASFEO_LIBRARY
                          BLASFEO_INCLUDE_DIR HPIPM_BOUND_MASKS)
endif ()

# ----- Export the target --------------------------------------------------- #
if (HPIPM_FOUND)
    set(HPIPM_INCLUDE_DIRS ${HPIPM_INCLUDE_DIR} ${BLASFEO_INCLUDE_DIR})
    set(HPIPM_LIBRARIES ${HPIPM_LIBRARY} ${BLASFEO_LIBRARY})

    if (NOT TARGET HPIPM::BLASFEO)
        add_library(HPIPM::BLASFEO UNKNOWN IMPORTED)
        set_target_properties(
                HPIPM::BLASFEO PROPERTIES
                IMPORTED_LOCATION "${BLASFEO_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${BLASFEO_INCLUDE_DIR}")
    endif ()

    if (NOT TARGET HPIPM::HPIPM)
        add_library(HPIPM::HPIPM UNKNOWN IMPORTED)
        set_target_properties(
                HPIPM::HPIPM PROPERTIES
                IMPORTED_LOCATION "${HPIPM_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${HPIPM_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES HPIPM::BLASFEO)
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(HPIPM_INCLUDE_DIR HPIPM_LIBRARY BLASFEO_INCLUDE_DIR
                 BLASFEO_LIBRARY)

# --------------------------------------------------------------------------- #
