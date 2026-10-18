# --------------------------------------------------------------------------- #
#    CMake find module for qpOASES                                            #
#                                                                             #
#    This module finds qpOASES include directories and libraries.             #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(qpOASES [version] [REQUIRED])                           #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        qpOASES_FOUND         - True if headers are found                    #
#        qpOASES_INCLUDE_DIRS  - Include directories                          #
#        qpOASES_LIBRARIES     - Libraries to be linked                       #
#        qpOASES_VERSION       - The version, if known                        #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        qpOASES_ROOT          - Custom path to qpOASES                       #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        qpOASES::qpOASES                                                     #
#                                                                             #
#    This find module is provided because the CMake configuration file that   #
#    qpOASES installs defines no target, only variables, while the one of     #
#    its vcpkg port defines qpOASES::qpOASES: the configuration file is       #
#    used if there is one, and the headers and the library are looked for     #
#    otherwise. On Windows qpOASES is a static library only.                  #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Check if already in cache
if (qpOASES_INCLUDE_DIR AND qpOASES_LIBRARY)
    set(qpOASES_FOUND TRUE)
endif ()

if (NOT qpOASES_FOUND)

    # ----- Try the CMake configuration file of qpOASES --------------------- #
    # The headers carry no version number, which is only known from it
    find_package(qpOASES CONFIG QUIET HINTS ${qpOASES_ROOT})

    if (TARGET qpOASES::qpOASES)
        get_target_property(qpOASES_INCLUDE_DIR qpOASES::qpOASES
                INTERFACE_INCLUDE_DIRECTORIES)
        set(qpOASES_LIBRARY qpOASES::qpOASES)
    endif ()

    # ----- Find the qpOASES include directory ------------------------------ #
    # The main header is <prefix>/include/qpOASES.hpp, the others being in
    # <prefix>/include/qpOASES
    find_path(qpOASES_INCLUDE_DIR
            NAMES qpOASES.hpp
            PATHS ${qpOASES_ROOT}
            PATH_SUFFIXES include
            DOC "qpOASES include directory.")

    # ----- Find the qpOASES library ---------------------------------------- #
    if (UNIX)
        find_library(qpOASES_LIBRARY
                NAMES qpOASES
                PATHS ${qpOASES_ROOT}
                PATH_SUFFIXES lib lib64
                DOC "qpOASES library.")
    elseif (WIN32)
        find_library(qpOASES_LIBRARY
                NAMES qpOASES libqpOASES
                PATHS
                ${qpOASES_ROOT}/lib
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
                $ENV{LIBRARY_LIB}
                NO_DEFAULT_PATH
                DOC "qpOASES library.")
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
    # set. REQUIRED_VARS should be cache entries and not output variables
    # (see the documentation of the FindPackageHandleStandardArgs module).
    find_package_handle_standard_args(
            qpOASES
            REQUIRED_VARS qpOASES_LIBRARY qpOASES_INCLUDE_DIR
            VERSION_VAR qpOASES_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (qpOASES_FOUND)
    set(qpOASES_INCLUDE_DIRS ${qpOASES_INCLUDE_DIR})
    set(qpOASES_LIBRARIES ${qpOASES_LIBRARY})

    if (NOT TARGET qpOASES::qpOASES)
        add_library(qpOASES::qpOASES UNKNOWN IMPORTED)
        set_target_properties(
                qpOASES::qpOASES PROPERTIES
                IMPORTED_LOCATION "${qpOASES_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${qpOASES_INCLUDE_DIRS}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(qpOASES_INCLUDE_DIR qpOASES_LIBRARY)

# --------------------------------------------------------------------------- #
