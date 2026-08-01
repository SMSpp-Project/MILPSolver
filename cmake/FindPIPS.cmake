# --------------------------------------------------------------------------- #
#    CMake find module for PIPS-IPM++                                         #
#                                                                             #
#    This module finds PIPS-IPM++ include directories and libraries.          #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(PIPS [REQUIRED])                                        #
#                                                                             #
#    This module mirrors the dependency information in makefile-libPIPS.      #
#    PIPS-IPM++ currently installs its library and a subset of its headers,   #
#    but does not install an exported CMake package.  PIPSMILPSolver also     #
#    needs headers that are only available in the PIPS-IPM++ source tree.     #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        PIPSIPMPP_ROOT             PIPS-IPM++ checkout root                  #
#        PIPS_ROOT                  Alternative checkout-root hint            #
#        PIPS_BUILD_DIR             PIPS-IPM++ CMake build directory          #
#        PIPS_PARDISO_DIR           Directory containing libpardiso           #
#                                                                             #
#    Optional solver switches:                                                #
#        PIPS_WITH_MUMPS            Link the MUMPS adapter (default ON)       #
#        PIPS_WITH_PARDISO          Link the Pardiso adapter (auto-detected)  #
#        PIPS_WITH_MA57             Link the MA57 adapter (auto-detected)     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        PIPS_FOUND                 True when all requested components exist  #
#        PIPS_INCLUDE_DIRS          PIPS-IPM++ source include directories     #
#        PIPS_LIBRARIES             PIPS-IPM++ and all link dependencies      #
#        PIPS_SOURCE_ROOT           Path to the PIPS-IPM subdirectory         #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        PIPS::PIPS                                                           #
#                                                                             #
#    This find module is provided because PIPS-IPM++ does not provide         #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
# External dependencies that PIPS-IPM++ needs at link time.
find_package(MPI QUIET COMPONENTS CXX)
find_package(OpenMP QUIET COMPONENTS CXX)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(Threads QUIET)

# Check if already in cache
if (PIPS_SOURCE_ROOT AND PIPS_LIBRARY)
    set(PIPS_FOUND TRUE)
endif ()

if (NOT PIPS_FOUND)

    # ----- Find the PIPS-IPM++ source tree --------------------------------- #
    # PIPSMILPSolver uses internal PIPS-IPM++ headers, so locate the source
    # tree rather than only an installation prefix.  Both the checkout root
    # and its PIPS-IPM subdirectory are accepted as hints.
    find_path(PIPS_SOURCE_ROOT
            NAMES Core/Interface/PIPSIPMppInterface.hpp
            HINTS
            "${PIPSIPMPP_ROOT}/PIPS-IPM"
            "${PIPSIPMPP_ROOT}"
            "${PIPS_ROOT}/PIPS-IPM"
            "${PIPS_ROOT}"
            "$ENV{PIPSIPMPP_ROOT}/PIPS-IPM"
            "$ENV{PIPSIPMPP_ROOT}"
            "$ENV{PIPS_ROOT}/PIPS-IPM"
            "$ENV{PIPS_ROOT}"
            "$ENV{PIPSIPM_ROOT}/PIPS-IPM"
            "$ENV{PIPSIPM_ROOT}"
            PATHS
            "/opt/pips-ipmpp/PIPS-IPM"
            "$ENV{HOME}/pips-ipmpp/PIPS-IPM"
            DOC "PIPS-IPM++ PIPS-IPM source directory.")

    if (PIPS_SOURCE_ROOT)
        get_filename_component(_PIPSIPMPP_ROOT "${PIPS_SOURCE_ROOT}" DIRECTORY)
        set(PIPSIPMPP_ROOT "${_PIPSIPMPP_ROOT}")
        set(PIPSIPMPP_ROOT "${_PIPSIPMPP_ROOT}" CACHE PATH
                "PIPS-IPM++ checkout root.")
    endif ()

    if (NOT PIPS_BUILD_DIR AND PIPSIPMPP_ROOT)
        set(PIPS_BUILD_DIR "${PIPSIPMPP_ROOT}/build" CACHE PATH
                "PIPS-IPM++ CMake build directory.")
    endif ()

    if (NOT PIPS_CORE_LIB_DIR AND PIPS_BUILD_DIR)
        set(PIPS_CORE_LIB_DIR "${PIPS_BUILD_DIR}/PIPS-IPM/Core" CACHE PATH
                "Directory containing the PIPS-IPM++ core library.")
    endif ()

    if (NOT PIPS_LINEAR_SOLVERS_LIB_DIR AND PIPS_CORE_LIB_DIR)
        set(PIPS_LINEAR_SOLVERS_LIB_DIR
                "${PIPS_CORE_LIB_DIR}/LinearSolvers" CACHE PATH
                "Directory with the PIPS-IPM++ solver adapter libraries.")
    endif ()

    if (NOT PIPS_PARDISO_DIR AND PIPSIPMPP_ROOT)
        if (DEFINED ENV{PARDISO_DIR} AND NOT "$ENV{PARDISO_DIR}" STREQUAL "")
            set(_PIPS_PARDISO_DIR_DEFAULT "$ENV{PARDISO_DIR}")
        else ()
            set(_PIPS_PARDISO_DIR_DEFAULT
                    "${PIPSIPMPP_ROOT}/ThirdPartyLibs/PARDISO/src")
        endif ()
        set(PIPS_PARDISO_DIR "${_PIPS_PARDISO_DIR_DEFAULT}" CACHE PATH
                "Directory containing the Pardiso library used by PIPS-IPM++.")
    endif ()

    # ----- Find the PIPS-IPM++ and adapter libraries ----------------------- #
    find_library(PIPS_LIBRARY
            NAMES pips-ipmpp
            HINTS "${PIPS_CORE_LIB_DIR}"
            DOC "PIPS-IPM++ core library.")

    find_library(PIPS_MUMPS_SOLVER_LIBRARY
            NAMES mumps_solver
            HINTS "${PIPS_LINEAR_SOLVERS_LIB_DIR}"
            DOC "PIPS-IPM++ MUMPS adapter library.")

    find_library(PIPS_PARDISO_SOLVER_LIBRARY
            NAMES pardiso_solver
            HINTS "${PIPS_LINEAR_SOLVERS_LIB_DIR}"
            DOC "PIPS-IPM++ Pardiso adapter library.")

    find_library(PIPS_PARDISO_ABSTRACT_LIBRARY
            NAMES pardiso_abstract
            HINTS "${PIPS_LINEAR_SOLVERS_LIB_DIR}"
            DOC "PIPS-IPM++ abstract Pardiso adapter library.")

    find_library(PIPS_MA57_LIBRARY
            NAMES ma57
            HINTS "${PIPS_LINEAR_SOLVERS_LIB_DIR}"
            DOC "PIPS-IPM++ MA57 adapter library.")

    find_library(PIPS_METIS_ADAPTER_LIBRARY
            NAMES metis_adapter
            HINTS "${PIPS_LINEAR_SOLVERS_LIB_DIR}"
            DOC "PIPS-IPM++ METIS adapter library.")

    # ----- Find the external solver libraries ------------------------------ #
    find_library(PIPS_MUMPS_D_LIBRARY NAMES dmumps
            DOC "MUMPS double-precision library.")
    find_library(PIPS_MUMPS_COMMON_LIBRARY NAMES mumps_common
            DOC "MUMPS common library.")
    find_library(PIPS_PORD_LIBRARY NAMES pord
            DOC "MUMPS PORD library.")
    find_library(PIPS_PARDISO_LIBRARY NAMES pardiso
            HINTS "${PIPS_PARDISO_DIR}"
            DOC "Pardiso library used by PIPS-IPM++.")
    find_library(PIPS_METIS_LIBRARY NAMES metis
            DOC "METIS library used by the PIPS-IPM++ MA57 adapter.")
    find_library(PIPS_MPI_MPIFH_LIBRARY NAMES mpi_mpifh
            DOC "MPI Fortran mpif.h compatibility library.")
    find_library(PIPS_GFORTRAN_LIBRARY
            NAMES gfortran libgfortran.so.5
            HINTS ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}
            DOC "GNU Fortran runtime library.")

    # ----- Solver adapter switches ----------------------------------------- #
    option(PIPS_WITH_MUMPS "Link the PIPS-IPM++ MUMPS solver adapter." ON)

    if (PIPS_PARDISO_SOLVER_LIBRARY AND PIPS_PARDISO_ABSTRACT_LIBRARY
            AND PIPS_PARDISO_LIBRARY)
        set(_PIPS_WITH_PARDISO_DEFAULT ON)
    else ()
        set(_PIPS_WITH_PARDISO_DEFAULT OFF)
    endif ()
    option(PIPS_WITH_PARDISO "Link the PIPS-IPM++ Pardiso solver adapter."
            ${_PIPS_WITH_PARDISO_DEFAULT})

    if (PIPS_MA57_LIBRARY AND PIPS_METIS_ADAPTER_LIBRARY)
        set(_PIPS_WITH_MA57_DEFAULT ON)
    else ()
        set(_PIPS_WITH_MA57_DEFAULT OFF)
    endif ()
    option(PIPS_WITH_MA57 "Link the PIPS-IPM++ MA57 solver adapter."
            ${_PIPS_WITH_MA57_DEFAULT})

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS
    # are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    set(_PIPS_REQUIRED_VARS
            PIPS_SOURCE_ROOT
            PIPS_LIBRARY
            MPI_CXX_FOUND
            OpenMP_CXX_FOUND
            BLAS_FOUND
            LAPACK_FOUND
            Threads_FOUND)

    if (PIPS_WITH_MUMPS)
        list(APPEND _PIPS_REQUIRED_VARS
                PIPS_MUMPS_SOLVER_LIBRARY
                PIPS_MUMPS_D_LIBRARY
                PIPS_MUMPS_COMMON_LIBRARY
                PIPS_PORD_LIBRARY)
    endif ()

    if (PIPS_WITH_PARDISO)
        list(APPEND _PIPS_REQUIRED_VARS
                PIPS_PARDISO_SOLVER_LIBRARY
                PIPS_PARDISO_ABSTRACT_LIBRARY
                PIPS_PARDISO_LIBRARY)
    endif ()

    if (PIPS_WITH_MA57)
        list(APPEND _PIPS_REQUIRED_VARS
                PIPS_MA57_LIBRARY
                PIPS_METIS_ADAPTER_LIBRARY
                PIPS_METIS_LIBRARY)
    endif ()

    find_package_handle_standard_args(PIPS
            REQUIRED_VARS ${_PIPS_REQUIRED_VARS})
endif ()

# ----- Export the target --------------------------------------------------- #
if (PIPS_FOUND)
    set(PIPS_INCLUDE_DIRS
            "${PIPS_SOURCE_ROOT}/Core/Interface"
            "${PIPS_SOURCE_ROOT}/Core/Options"
            "${PIPS_SOURCE_ROOT}/Core/Utilities"
            "${PIPS_SOURCE_ROOT}/Core/InteriorPointMethod"
            "${PIPS_SOURCE_ROOT}/Core/LinearAlgebra"
            "${PIPS_SOURCE_ROOT}/Core"
            "${PIPS_SOURCE_ROOT}/Drivers/CallbackExample"
            "${PIPS_SOURCE_ROOT}/Core/Readers/Distributed")

    set(_PIPS_LINK_LIBRARIES)

    # Keep the static adapter libraries before their external dependencies,
    # matching makefile-libPIPS and the PIPS-IPM++ executable link line.
    if (PIPS_WITH_MUMPS)
        list(APPEND _PIPS_LINK_LIBRARIES
                "${PIPS_MUMPS_SOLVER_LIBRARY}")
    endif ()

    if (PIPS_WITH_PARDISO)
        list(APPEND _PIPS_LINK_LIBRARIES
                "${PIPS_PARDISO_SOLVER_LIBRARY}"
                "${PIPS_PARDISO_ABSTRACT_LIBRARY}")
    endif ()

    if (PIPS_WITH_MA57)
        list(APPEND _PIPS_LINK_LIBRARIES
                "${PIPS_MA57_LIBRARY}"
                "${PIPS_METIS_ADAPTER_LIBRARY}")
    endif ()

    if (PIPS_WITH_MUMPS)
        list(APPEND _PIPS_LINK_LIBRARIES
                "${PIPS_MUMPS_D_LIBRARY}"
                "${PIPS_MUMPS_COMMON_LIBRARY}"
                "${PIPS_PORD_LIBRARY}")
    endif ()

    if (PIPS_WITH_PARDISO)
        list(APPEND _PIPS_LINK_LIBRARIES "${PIPS_PARDISO_LIBRARY}")
    endif ()

    if (PIPS_WITH_MA57)
        list(APPEND _PIPS_LINK_LIBRARIES "${PIPS_METIS_LIBRARY}")
    endif ()

    if (PIPS_MPI_MPIFH_LIBRARY AND UNIX AND NOT APPLE)
        # This is the CMake equivalent of the --no-as-needed section in
        # makefile-libPIPS, needed by some shared MUMPS/OpenMPI combinations.
        list(APPEND _PIPS_LINK_LIBRARIES
                -Wl,--no-as-needed
                "${PIPS_MPI_MPIFH_LIBRARY}"
                MPI::MPI_CXX
                -Wl,--as-needed)
    else ()
        if (PIPS_MPI_MPIFH_LIBRARY)
            list(APPEND _PIPS_LINK_LIBRARIES "${PIPS_MPI_MPIFH_LIBRARY}")
        endif ()
        list(APPEND _PIPS_LINK_LIBRARIES MPI::MPI_CXX)
    endif ()

    list(APPEND _PIPS_LINK_LIBRARIES
            OpenMP::OpenMP_CXX
            LAPACK::LAPACK
            BLAS::BLAS)

    if (PIPS_GFORTRAN_LIBRARY)
        list(APPEND _PIPS_LINK_LIBRARIES "${PIPS_GFORTRAN_LIBRARY}")
    else ()
        # Compilers normally know their own Fortran runtime search directory.
        list(APPEND _PIPS_LINK_LIBRARIES gfortran)
    endif ()

    list(APPEND _PIPS_LINK_LIBRARIES Threads::Threads m ${CMAKE_DL_LIBS})

    set(PIPS_LIBRARIES "${PIPS_LIBRARY};${_PIPS_LINK_LIBRARIES}")

    if (NOT TARGET PIPS::PIPS)
        add_library(PIPS::PIPS UNKNOWN IMPORTED)
        set_target_properties(PIPS::PIPS PROPERTIES
                IMPORTED_LOCATION "${PIPS_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${PIPS_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${_PIPS_LINK_LIBRARIES}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(
        PIPS_SOURCE_ROOT
        PIPS_CORE_LIB_DIR
        PIPS_LINEAR_SOLVERS_LIB_DIR
        PIPS_LIBRARY
        PIPS_MUMPS_SOLVER_LIBRARY
        PIPS_PARDISO_SOLVER_LIBRARY
        PIPS_PARDISO_ABSTRACT_LIBRARY
        PIPS_MA57_LIBRARY
        PIPS_METIS_ADAPTER_LIBRARY
        PIPS_MUMPS_D_LIBRARY
        PIPS_MUMPS_COMMON_LIBRARY
        PIPS_PORD_LIBRARY
        PIPS_PARDISO_LIBRARY
        PIPS_METIS_LIBRARY
        PIPS_MPI_MPIFH_LIBRARY
        PIPS_GFORTRAN_LIBRARY)

unset(_PIPSIPMPP_ROOT)
unset(_PIPS_PARDISO_DIR_DEFAULT)
unset(_PIPS_REQUIRED_VARS)
unset(_PIPS_LINK_LIBRARIES)
unset(_PIPS_WITH_PARDISO_DEFAULT)
unset(_PIPS_WITH_MA57_DEFAULT)

# --------------------------------------------------------------------------- #
