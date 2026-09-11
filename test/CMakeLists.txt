# --------------------------------------------------------------------------- #
#    CMake file for MILPSolver/test                                           #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
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

# Utility: build a test executable linked against MILPSolver, copying the
# runtime DLLs next to it on Windows.
function(add_milp_test target_name source)
    add_executable(${target_name} ${source})
    target_link_libraries(${target_name} PRIVATE SMS++::MILPSolver)

    if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
        # On Windows with shared libraries enabled copies all runtime DLL
        # dependencies next to the executable in the build directory
        # (POST_BUILD), using CMake's TARGET_RUNTIME_DLLS generator expression
        add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${target_name}>
                $<TARGET_FILE_DIR:${target_name}>
                COMMAND_EXPAND_LISTS
                COMMENT "Copying runtime DLLs next to ${target_name}.exe")
    endif ()
endfunction()

# Utility: register a batch file as a ctest driving the given executable.
function(add_batch_test exe batch)
    if (UNIX)
        add_test(NAME ${exe}/${batch}
                 COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/${batch}
                         ${CMAKE_CURRENT_BINARY_DIR}/${exe}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    elseif (WIN32)
        add_test(NAME ${exe}/${batch}
                 # using Git for Windows that includes Git Bash
                 COMMAND "C:/Program Files/Git/bin/bash.exe"
                         "${CMAKE_CURRENT_SOURCE_DIR}/${batch}"
                         "$<TARGET_FILE:${exe}>"
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif ()
endfunction()

# ----- test_cuts ----------------------------------------------------------- #

add_milp_test(test_cuts test_cuts.cpp)
add_batch_test(test_cuts batch-cuts)

# ----- test_dual ----------------------------------------------------------- #

add_milp_test(test_dual test_dual.cpp)
add_test(NAME test_dual
         COMMAND test_dual
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_dynamic -------------------------------------------------------- #

add_milp_test(test_dynamic test_dynamic.cpp)
add_batch_test(test_dynamic batch-dynamic)
add_batch_test(test_dynamic batch-dynamic-L)

# --------------------------------------------------------------------------- #
