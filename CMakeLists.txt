cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MILPSolver DESCRIPTION "MILPSolver for SMS++")

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Requirements -----------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(CPLEX REQUIRED)

# ----- Includes ---------------------------------------------------------------
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${CPLEX_INCLUDE_DIR})
include_directories("${PROJECT_INCLUDE_DIR}")
# TODO: Write and use a FindSMS++.cmake
include_directories(../SMS++/include)

# ----- SOURCE FILES -----------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    ${PROJECT_SOURCE_DIR}/MILPSolver.cpp
    ${PROJECT_SOURCE_DIR}/CPXMILPSolver.cpp
    )

# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
