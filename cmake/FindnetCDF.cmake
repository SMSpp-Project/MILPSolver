# --------------------------------------------------------------------------- #
#    CMake find module for NetCDF                                             #
#                                                                             #
#    This find module is provided because NetCDF and NetCDF-C++ support for   #
#    both CMake and PkgConfig configuration files is inconsistent between     #
#    original distribution, Homebrew formula and .deb packages.               #
#                                                                             #
#    Accepts the following HINTS:                                             #
#                                                                             #
#    - NETCDF_INC - Custom path to NetCDF headers                             #
#    - NETCDF_LIB - Custom path to NetCDF libraries                           #
#                                                                             #
#    Provides (at least) the following variables:                             #
#                                                                             #
#    - netCDF_FOUND - Whether NetCDF was found or not                         #
#    - netCDF_INCLUDE_DIRS - Include directories for NetCDF-C                 #
#    - netCDF_LIBRARIES - Libraries necessary to use NetCDF-C                 #
#    - netCDF_VERSION - The version of NetCDF-C found                         #
#                                                                             #
#    - netCDFCxx_FOUND - Whether NetCDF-C++ was found or not                  #
#    - netCDFCxx_INCLUDE_DIRS - Include directories for NetCDF-C++            #
#    - netCDFCxx_LIBRARIES - Libraries necessary to use NetCDF-C++            #
#    - netCDFCxx_VERSION - The version of NetCDF-C++ found                    #
#                                                                             #
#    - NetCDF::NetCDF - A target to use with target_link_libraries()          #
#    - NetCDF::NetCDFCxx - A target to use with target_link_libraries()       #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# ----- NetCDF-C ------------------------------------------------------------ #

# Try to find a CMake configuration file
find_package(netCDF CONFIG QUIET)
if (netCDF_FOUND)
    # Forward the variables in a consistent way
    set(netCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
    set(netCDF_VERSION "${NetCDFVersion}")

    # The original CMake config file exports a target without a namespace
    if (NOT TARGET NetCDF::NetCDF)
        add_library(NetCDF::NetCDF INTERFACE IMPORTED)
        set_target_properties(
                NetCDF::NetCDF PROPERTIES
                INTERFACE_LINK_LIBRARIES "${netCDF_LIBRARIES}")
    endif ()
endif ()

# Try to find a PkgConfig module
if (NOT netCDF_FOUND)

    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(_NetCDF QUIET netcdf IMPORTED_TARGET)
        if (_NetCDF_FOUND)
            # Forward the variables in a consistent way
            set(netCDF_FOUND "${_NetCDF_FOUND}")
            set(netCDF_INCLUDE_DIRS "${_NetCDF_INCLUDE_DIRS}")
            set(netCDF_INCLUDE_DIR "${_NetCDF_INCLUDE_DIRS}")
            set(netCDF_LIBRARIES "${_NetCDF_LIBRARIES}")
            set(netCDF_VERSION "${_NetCDF_VERSION}")

            # Export a target
            if (NOT TARGET NetCDF::NetCDF)
                add_library(NetCDF::NetCDF INTERFACE IMPORTED)
                set_target_properties(
                        NetCDF::NetCDF PROPERTIES
                        INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDF")
            endif ()
        endif ()
    endif ()
endif ()


# Try to find it manually
if (NOT netCDF_FOUND)

    find_path(netCDF_INCLUDE_DIR
              NAMES netcdf.h
              HINTS ${NETCDF_INC}
              DOC "NetCDF include directories")
    mark_as_advanced(netCDF_INCLUDE_DIR)

    find_library(netCDF_LIBRARY
                 NAMES netcdf
                 HINTS ${NETCDF_LIB}
                 DOC "NetCDF library")
    mark_as_advanced(netCDF_LIBRARY)

    # Parse version
    if (netCDF_INCLUDE_DIR)
        file(STRINGS "${netCDF_INCLUDE_DIR}/netcdf_meta.h" _netcdf_version_lines
             REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
        string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
        string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
        set(netCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
        unset(_netcdf_version_major)
        unset(_netcdf_version_minor)
        unset(_netcdf_version_patch)
        unset(_netcdf_version_note)
        unset(_netcdf_version_lines)
    endif ()

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(
            NetCDF
            REQUIRED_VARS netCDF_LIBRARY netCDF_INCLUDE_DIR
            VERSION_VAR netCDF_VERSION)

    if (netCDF_FOUND)
        set(netCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
        set(netCDF_LIBRARIES "${netCDF_LIBRARY}")

        if (NOT TARGET NetCDF::NetCDF)
            add_library(NetCDF::NetCDF UNKNOWN IMPORTED)
            set_target_properties(
                    NetCDF::NetCDF PROPERTIES
                    IMPORTED_LOCATION "${netCDF_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${netCDF_INCLUDE_DIR}")
        endif ()
    endif ()
endif ()

# ----- NetCDF-C++ ---------------------------------------------------------- #

# Try to find a PkgConfig module
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
    pkg_check_modules(_NetCDFCxx QUIET netcdf-cxx4 IMPORTED_TARGET)
    if (_NetCDFCxx_FOUND)
        # Forward the variables in a consistent way
        set(netCDFCxx_FOUND "${_NetCDFCxx_FOUND}")
        set(netCDFCxx_INCLUDE_DIRS "${_NetCDFCxx_INCLUDE_DIRS}")
        set(netCDFCxx_INCLUDE_DIR "${_NetCDFCxx_INCLUDE_DIRS}")
        set(netCDFCxx_LIBRARIES "${_NetCDFCxx_LIBRARIES}")
        set(netCDFCxx_VERSION "${_NetCDFCxx_VERSION}")

        # Export a target
        if (NOT TARGET NetCDF::NetCDFCxx)
            add_library(NetCDF::NetCDFCxx INTERFACE IMPORTED)
            set_target_properties(
                    NetCDF::NetCDFCxx PROPERTIES
                    INTERFACE_LINK_LIBRARIES "PkgConfig::_NetCDFCxx")
        endif ()
    endif ()
endif ()

# Try to find it manually
if (NOT netCDFCxx_FOUND)

    find_path(netCDFCxx_INCLUDE_DIR
              NAMES netcdf
              HINTS ${NETCDF_INC}
              DOC "NetCDF-C++ include directories")
    mark_as_advanced(netCDFCxx_INCLUDE_DIR)

    find_library(netCDFCxx_LIBRARY
                 NAMES netcdf-cxx4 netcdf_c++4
                 HINTS ${NETCDF_LIB}
                 DOC "NetCDF-C++ library")
    mark_as_advanced(netCDFCxx_LIBRARY)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(
            netCDFCxx
            REQUIRED_VARS netCDFCxx_LIBRARY netCDFCxx_INCLUDE_DIR
            VERSION_VAR netCDF_VERSION)

    if (netCDFCxx_FOUND)
        set(netCDFCxx_INCLUDE_DIRS "${netCDFCxx_INCLUDE_DIR}")
        set(netCDFCxx_LIBRARIES "${netCDFCxx_LIBRARY}")

        if (NOT TARGET NetCDF::NetCDFCxx)
            add_library(NetCDF::NetCDFCxx UNKNOWN IMPORTED)
            set_target_properties(
                    NetCDF::NetCDFCxx PROPERTIES
                    IMPORTED_LOCATION "${netCDFCxx_LIBRARY}"
                    INTERFACE_INCLUDE_DIRECTORIES "${netCDFCxx_INCLUDE_DIR}")
        endif ()
    endif ()
endif ()
