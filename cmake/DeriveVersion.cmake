# --------------------------------------------------------------------------- #
# DeriveVersion.cmake                                                          #
#                                                                              #
# Derive the module version from git, so that cutting a release is just        #
# `git tag x.y.z` with NO version number written by hand anywhere: the version #
# in project(), in the generated <Module>Config.h and in the tools' --version  #
# all follow the tag.                                                          #
#                                                                              #
# Resolution order:                                                            #
#   1. a working git clone: `git describe --tags --abbrev=0`;                  #
#   2. a source tarball: the VERSION file, which git fills in at `git archive` #
#      time via the `VERSION export-subst` .gitattributes entry (its content   #
#      is "$Format:%(describe:tags)$", substituted with the tag by git);       #
#   3. otherwise the sentinel "0.0.0" (no version info available).             #
#                                                                              #
# Call it BEFORE project():                                                    #
#     include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)           #
#     smspp_derive_version(MOD_VERSION)                                        #
#     project(<Module> VERSION ${MOD_VERSION} ...)                            #
# --------------------------------------------------------------------------- #

function(smspp_derive_version out_var)
    set(_version "0.0.0")
    set(_resolved FALSE)

    # 1) working git clone
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(
                COMMAND "${GIT_EXECUTABLE}" describe --tags --abbrev=0
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE _tag
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _result)
        if(_result EQUAL 0 AND _tag)
            string(REGEX REPLACE "^v" "" _tag "${_tag}")
            if(_tag MATCHES "^[0-9]+(\\.[0-9]+)*$")
                set(_version "${_tag}")
                set(_resolved TRUE)
            endif()
        endif()
    endif()

    # 2) tarball: the VERSION file filled in by `git archive` (export-subst)
    if(NOT _resolved AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _file_version)
        string(STRIP "${_file_version}" _file_version)
        string(REGEX REPLACE "^v" "" _file_version "${_file_version}")
        # only accept it if git actually substituted the placeholder
        if(_file_version MATCHES "^[0-9]+(\\.[0-9]+)*$")
            set(_version "${_file_version}")
        endif()
    endif()

    set(${out_var} "${_version}" PARENT_SCOPE)
endfunction()
