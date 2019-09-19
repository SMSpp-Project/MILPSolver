cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MILPSolver DESCRIPTION "MILPSolver for SMS++")

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(CPLEX REQUIRED)

# ----- Custom includes --------------------------------------------------------
# TODO: When libraries will be installable, fix this
set(SMS++_INCLUDE_DIR ../SMS++/include)

# ----- Library ----------------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    ${PROJECT_SOURCE_DIR}/MILPSolver.cpp
    ${PROJECT_SOURCE_DIR}/CPXMILPSolver.cpp
    )

# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE_DIR})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${SMS++_INCLUDE_DIR})

# FIXME: Workaround for CPLEX
if (OVERRIDE_CPLEX_INC)
    target_include_directories(${PROJECT_NAME} PUBLIC ${OVERRIDE_CPLEX_INC})
else ()
    find_package(CPLEX REQUIRED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIR})
endif ()

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
