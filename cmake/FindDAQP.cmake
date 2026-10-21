# --------------------------------------------------------------------------- #
#    CMake find module for DAQP                                               #
#                                                                             #
#    This module finds DAQP include directories and libraries.                #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(DAQP [REQUIRED])                                        #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        DAQP_FOUND          - True if headers are found                      #
#        DAQP_INCLUDE_DIRS   - Include directories                            #
#        DAQP_LIBRARIES      - Libraries to be linked                         #
#        DAQP_DLL            - The runtime DLL, if any (Windows only)         #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        DAQP_ROOT           - Custom path to DAQP                            #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        DAQP::DAQP                                                           #
#                                                                             #
#    This find module is provided because the CMake configuration file of     #
#    DAQP is installed in <prefix>/cmake, where find_package() does not look  #
#    on Linux and macOS, and DAQP exposes no version number in its headers:   #
#    0.8.2 or later is told apart by a macro of its constants.h.              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Check if already in cache
if (DAQP_INCLUDE_DIR AND DAQP_LIBRARY)
    set(DAQP_FOUND TRUE)
endif ()

if (NOT DAQP_FOUND)

    # ----- Find the DAQP include directory --------------------------------- #
    # The headers are in <prefix>/include/daqp and are included as daqp/*.h
    find_path(DAQP_INCLUDE_DIR
            NAMES daqp/api.h
            PATHS ${DAQP_ROOT}
            PATH_SUFFIXES include
            DOC "DAQP include directory.")

    # ----- Check that it is 0.8.2 or later --------------------------------- #
    # The headers carry no version number: DAQP_DEFAULT_EPS_PROX, which
    # DAQPMILPSolver uses, is in those of 0.8.2 and later only, hence
    # DAQP_EPS_PROX is required
    if (DAQP_INCLUDE_DIR AND EXISTS "${DAQP_INCLUDE_DIR}/daqp/constants.h")
        file(STRINGS "${DAQP_INCLUDE_DIR}/daqp/constants.h" DAQP_EPS_PROX
             REGEX "define[ \t]+DAQP_DEFAULT_EPS_PROX")
    endif ()

    # ----- Find the DAQP library ------------------------------------------- #
    if (UNIX)
        find_library(DAQP_LIBRARY
                NAMES daqp
                PATHS ${DAQP_ROOT}/lib
                DOC "DAQP library.")
    elseif (WIN32)
        find_library(DAQP_LIBRARY
                NAMES daqp
                PATHS
                ${DAQP_ROOT}/lib
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
                $ENV{LIBRARY_LIB}
                NO_DEFAULT_PATH
                DOC "DAQP library.")

        find_file(DAQP_DLL
                NAMES daqp.dll libdaqp.dll
                PATHS
                ${DAQP_ROOT}/bin
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin
                $ENV{LIBRARY_BIN}
                NO_DEFAULT_PATH
                DOC "DAQP runtime DLL.")
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
    # set. REQUIRED_VARS should be cache entries and not output variables.
    # See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    # on Windows the DLL is not required: a static build (e.g., a static
    # vcpkg triplet) has none, and the library is then a static one
    find_package_handle_standard_args(
            DAQP
            REQUIRED_VARS DAQP_LIBRARY DAQP_INCLUDE_DIR DAQP_EPS_PROX)
endif ()

# ----- Export the target --------------------------------------------------- #
if (DAQP_FOUND)
    set(DAQP_INCLUDE_DIRS ${DAQP_INCLUDE_DIR})
    set(DAQP_LIBRARIES ${DAQP_LIBRARY})

    if (NOT TARGET DAQP::DAQP)
        if (WIN32 AND DAQP_DLL)
            add_library(DAQP::DAQP SHARED IMPORTED)
            set_target_properties(
                    DAQP::DAQP PROPERTIES
                    IMPORTED_IMPLIB "${DAQP_LIBRARY}"
                    IMPORTED_LOCATION "${DAQP_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${DAQP_INCLUDE_DIRS}")
        else ()
            add_library(DAQP::DAQP UNKNOWN IMPORTED)
            set_target_properties(
                    DAQP::DAQP PROPERTIES
                    IMPORTED_LOCATION "${DAQP_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${DAQP_INCLUDE_DIRS}")
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
if (WIN32)
    mark_as_advanced(DAQP_INCLUDE_DIR DAQP_LIBRARY DAQP_DLL)
else ()
    mark_as_advanced(DAQP_INCLUDE_DIR DAQP_LIBRARY)
endif ()

# --------------------------------------------------------------------------- #
