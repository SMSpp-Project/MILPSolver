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

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# Utility: build a test executable linked against MILPSolver, copying the
# runtime DLLs next to it on Windows.
function(add_milp_test target_name source)
    add_executable(${target_name} ${source})
    target_link_libraries(${target_name} PRIVATE SMS++::MILPSolver)

    smspp_copy_runtime(${target_name})
endfunction()

# Utility: register a batch file as a ctest driving the given executable.
function(add_batch_test exe batch)
    add_test(NAME ${exe}/${batch}
             COMMAND ${SMSpp_BASH}
                     "${CMAKE_CURRENT_SOURCE_DIR}/${batch}"
                     "$<TARGET_FILE:${exe}>"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

# ----- test_cuts ----------------------------------------------------------- #

add_milp_test(test_cuts test_cuts.cpp)
add_batch_test(test_cuts batch-cuts)

# ----- test_batch_fix ------------------------------------------------------ #

add_milp_test(test_batch_fix test_batch_fix.cpp)
add_test(NAME test_batch_fix
         COMMAND test_batch_fix
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_dual ----------------------------------------------------------- #

add_milp_test(test_dual test_dual.cpp)
add_test(NAME test_dual
         COMMAND test_dual
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_quad_modify ---------------------------------------------------- #

add_milp_test(test_quad_modify test_quad_modify.cpp)
add_test(NAME test_quad_modify
         COMMAND test_quad_modify
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_quad_row_modify ------------------------------------------------ #

add_milp_test(test_quad_row_modify test_quad_row_modify.cpp)
add_test(NAME test_quad_row_modify
         COMMAND test_quad_row_modify
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_qp ------------------------------------------------------------- #

add_milp_test(test_qp test_qp.cpp)
add_test(NAME test_qp
         COMMAND test_qp
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_farkas --------------------------------------------------------- #

add_milp_test(test_farkas test_farkas.cpp)
add_test(NAME test_farkas
         COMMAND test_farkas
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_groups --------------------------------------------------------- #

add_milp_test(test_groups test_groups.cpp)
add_test(NAME test_groups
         COMMAND test_groups
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_reduced_costs -------------------------------------------------- #

add_milp_test(test_reduced_costs test_reduced_costs.cpp)
add_test(NAME test_reduced_costs
         COMMAND test_reduced_costs
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_stopped_bounds ------------------------------------------------- #

add_milp_test(test_stopped_bounds test_stopped_bounds.cpp)
add_test(NAME test_stopped_bounds
         COMMAND test_stopped_bounds
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_free_rows ------------------------------------------------------ #

add_milp_test(test_free_rows test_free_rows.cpp)
add_test(NAME test_free_rows
         COMMAND test_free_rows
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_ipm_free_rows -------------------------------------------------- #

add_milp_test(test_ipm_free_rows test_ipm_free_rows.cpp)
add_test(NAME test_ipm_free_rows
         COMMAND test_ipm_free_rows
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_ipm_fallback --------------------------------------------------- #

add_milp_test(test_ipm_fallback test_ipm_fallback.cpp)
add_test(NAME test_ipm_fallback
         COMMAND test_ipm_fallback
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_cutoff_status -------------------------------------------------- #

add_milp_test(test_cutoff_status test_cutoff_status.cpp)
add_test(NAME test_cutoff_status
         COMMAND test_cutoff_status
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_free_workspace ------------------------------------------------- #

add_milp_test(test_free_workspace test_free_workspace.cpp)
add_test(NAME test_free_workspace
         COMMAND test_free_workspace
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_partial_solution ----------------------------------------------- #

add_milp_test(test_partial_solution test_partial_solution.cpp)
add_test(NAME test_partial_solution
         COMMAND test_partial_solution
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_qcp_dual ------------------------------------------------------- #

add_milp_test(test_qcp_dual test_qcp_dual.cpp)
add_test(NAME test_qcp_dual
         COMMAND test_qcp_dual
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_fail_feas ------------------------------------------------------ #

add_milp_test(test_fail_feas test_fail_feas.cpp)
add_test(NAME test_fail_feas
         COMMAND test_fail_feas
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ----- test_dynamic -------------------------------------------------------- #

add_milp_test(test_dynamic test_dynamic.cpp)
add_batch_test(test_dynamic batch-dynamic)
add_batch_test(test_dynamic batch-dynamic-L)

# ----- labels and the status of a test that has nothing to do -------------- #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
#
# A test whose configuration names more :MILPSolver than this build has gives
# up with 77 rather than failing [see keep_available_Solvers() of
# test_common.h], which is what a machine without CPLEX does to test_dynamic,
# and ctest reads that status as "skipped".
get_property(_tests DIRECTORY PROPERTY TESTS)
set_tests_properties(${_tests} PROPERTIES LABELS "${modName}"
                                          SKIP_RETURN_CODE 77)

# --------------------------------------------------------------------------- #
