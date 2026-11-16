# --------------------------------------------------------------------------- #
#    CMake find module for Clarabel                                           #
#                                                                             #
#    This module finds the include directories and the library of the C       #
#    interface of Clarabel (Clarabel.cpp).                                    #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(Clarabel [REQUIRED])                                    #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        Clarabel_FOUND        - True if headers are found                    #
#        Clarabel_INCLUDE_DIRS - Include directories                          #
#        Clarabel_LIBRARIES    - Libraries to be linked                       #
#        Clarabel_DLL          - The runtime DLL, if any (Windows only)       #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        Clarabel_ROOT         - Custom path to Clarabel                      #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        Clarabel::Clarabel                                                   #
#                                                                             #
#    This find module is provided because Clarabel.cpp installs nothing: its  #
#    CMake builds the Rust library libclarabel_c with cargo inside the        #
#    source tree, and the headers stay in its include directory. Clarabel is  #
#    looked for in a prefix where clarabel.h and its c/ subdirectory are in   #
#    include and libclarabel_c is in lib; a static libclarabel_c also needs   #
#    the system libraries of the Rust standard library.                       #
#    Clarabel.cpp carries no version number, and the clarabel.h looked for    #
#    is there from 0.11.1 (Clarabel.h before), the first version accepted.    #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_Clarabel_cached Clarabel_INCLUDE_DIR Clarabel_LIBRARY Clarabel_DLL)
foreach (_var ${_Clarabel_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_Clarabel_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
if (_Clarabel_stale)
    message(STATUS "Cached Clarabel paths not found, searching again")
    foreach (_var ${_Clarabel_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_Clarabel_cached)
unset(_Clarabel_stale)
unset(_var)
unset(_path)

# Check if already in cache
if (Clarabel_INCLUDE_DIR AND Clarabel_LIBRARY)
    set(Clarabel_FOUND TRUE)
endif ()

if (NOT Clarabel_FOUND)

    # ----- Find the Clarabel include directory ----------------------------- #
    # clarabel.h includes the headers in its c/ subdirectory
    find_path(Clarabel_INCLUDE_DIR
            NAMES clarabel.h
            PATHS ${Clarabel_ROOT}
            PATH_SUFFIXES include include/clarabel
            DOC "Clarabel include directory.")

    # ----- Find the Clarabel library --------------------------------------- #
    # also where conda and vcpkg put it, and the DLL where there is one
    set(_Clarabel_dirs ${Clarabel_ROOT}/lib $ENV{LIBRARY_LIB})
    set(_Clarabel_bins ${Clarabel_ROOT}/bin ${Clarabel_ROOT}/lib
            $ENV{LIBRARY_BIN})
    if (VCPKG_INSTALLED_DIR)
        list(APPEND _Clarabel_dirs
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib)
        list(APPEND _Clarabel_bins
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin)
    endif ()
    find_library(Clarabel_LIBRARY
            NAMES clarabel_c clarabel_c${CMAKE_SHARED_LIBRARY_SUFFIX}
            PATHS ${_Clarabel_dirs}
            DOC "Clarabel library.")
    # the DLL, which a link with an import library has next to it
    if (NOT Clarabel_LIBRARY MATCHES
            "${CMAKE_SHARED_LIBRARY_SUFFIX}$")
        find_file(Clarabel_DLL
                NAMES clarabel_c${CMAKE_SHARED_LIBRARY_SUFFIX}
                PATHS ${_Clarabel_bins}
                NO_DEFAULT_PATH
                DOC "Clarabel runtime DLL.")
    endif ()
    unset(_Clarabel_dirs)
    unset(_Clarabel_bins)

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET and REQUIRED options passed to
    # find_package(). It also sets <PackageName>_FOUND if REQUIRED_VARS are
    # set. REQUIRED_VARS should be cache entries and not output variables.
    # See the FindPackageHandleStandardArgs module in the CMake manual.
    # the DLL is not required: a static build has none, and the library is
    # then a static one
    find_package_handle_standard_args(
            Clarabel
            REQUIRED_VARS Clarabel_LIBRARY Clarabel_INCLUDE_DIR)
endif ()

# ----- Export the target --------------------------------------------------- #
if (Clarabel_FOUND)
    set(Clarabel_INCLUDE_DIRS ${Clarabel_INCLUDE_DIR})
    set(Clarabel_LIBRARIES ${Clarabel_LIBRARY})

    if (NOT TARGET Clarabel::Clarabel)
        if (Clarabel_DLL)
            add_library(Clarabel::Clarabel SHARED IMPORTED)
            set_target_properties(
                    Clarabel::Clarabel PROPERTIES
                    IMPORTED_IMPLIB "${Clarabel_LIBRARY}"
                    IMPORTED_LOCATION "${Clarabel_DLL}"
                    INTERFACE_INCLUDE_DIRECTORIES "${Clarabel_INCLUDE_DIRS}")
        else ()
            add_library(Clarabel::Clarabel UNKNOWN IMPORTED)
            set_target_properties(
                    Clarabel::Clarabel PROPERTIES
                    IMPORTED_LOCATION "${Clarabel_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${Clarabel_INCLUDE_DIRS}")
        endif ()

        # the Rust standard library in a static libclarabel_c needs these:
        # the threads, the loader and those of the system libraries of Rust
        # that the system has (m, or ws2_32, userenv, ntdll and bcrypt)
        if (Clarabel_LIBRARY MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}$")
            find_package(Threads REQUIRED)
            set(clarabelSysLibs Threads::Threads ${CMAKE_DL_LIBS})
            foreach (_lib m ws2_32 userenv ntdll bcrypt)
                find_library(Clarabel_SYS_${_lib} ${_lib})
                mark_as_advanced(Clarabel_SYS_${_lib})
                if (Clarabel_SYS_${_lib})
                    list(APPEND clarabelSysLibs ${_lib})
                endif ()
            endforeach ()
            set_property(TARGET Clarabel::Clarabel APPEND PROPERTY
                    INTERFACE_LINK_LIBRARIES ${clarabelSysLibs})
        endif ()
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(Clarabel_INCLUDE_DIR Clarabel_LIBRARY Clarabel_DLL)

# --------------------------------------------------------------------------- #
