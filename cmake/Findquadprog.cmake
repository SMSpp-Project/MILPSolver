# --------------------------------------------------------------------------- #
#    CMake find module for quadprog                                           #
#                                                                             #
#    This module finds the quadprog library.                                  #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(quadprog [REQUIRED])                                    #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        quadprog_FOUND        - True if the library is found                 #
#        quadprog_LIBRARIES    - Libraries to be linked                       #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        quadprog_ROOT         - Custom path to quadprog                      #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        quadprog::quadprog                                                   #
#                                                                             #
#    quadprog is distributed as the C sources of a Python package, with no    #
#    header, no library and no CMake configuration: the library is the one    #
#    that INSTALL.sh compiles from those sources in <quadprog_ROOT>/lib, and  #
#    QuadProgMILPSolver declares the function it calls.                       #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Check if already in cache
if (quadprog_LIBRARY)
    set(quadprog_FOUND TRUE)
endif ()

if (NOT quadprog_FOUND)

    # ----- Find the quadprog library --------------------------------------- #
    find_library(quadprog_LIBRARY
            NAMES quadprog
            PATHS ${quadprog_ROOT}
            PATH_SUFFIXES lib lib64
            DOC "quadprog library.")

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
    # set. REQUIRED_VARS should be cache entries and not output variables
    # (see the documentation of the FindPackageHandleStandardArgs module).
    find_package_handle_standard_args(
            quadprog
            REQUIRED_VARS quadprog_LIBRARY)
endif ()

# ----- Export the target --------------------------------------------------- #
if (quadprog_FOUND)
    set(quadprog_LIBRARIES ${quadprog_LIBRARY})

    if (NOT TARGET quadprog::quadprog)
        add_library(quadprog::quadprog UNKNOWN IMPORTED)
        set_target_properties(
                quadprog::quadprog PROPERTIES
                IMPORTED_LOCATION "${quadprog_LIBRARY}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(quadprog_LIBRARY)

# --------------------------------------------------------------------------- #
