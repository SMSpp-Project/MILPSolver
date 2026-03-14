# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test_dynamic                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- MILPSolver_test_dynamic --------------------------------------------- #
add_executable(test_dynamic test_dynamic.cpp)
target_link_libraries(test_dynamic PRIVATE SMS++::MILPSolver)

if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
    # copy runtime DLLs next to the exe after build
    add_custom_command(TARGET test_dynamic POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy_if_different
                       $<TARGET_RUNTIME_DLLS:test_dynamic>
                       $<TARGET_FILE_DIR:test_dynamic>
                       COMMAND_EXPAND_LISTS
                       COMMENT "Copying runtime DLLs next to test_dynamic.exe")
endif()

# --------------------------------------------------------------------------- #
