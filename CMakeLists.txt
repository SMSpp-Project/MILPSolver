cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(MILPSolver DESCRIPTION "MILPSolver for SMS++")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(MILPSolver_SOURCES
    ${PROJECT_SOURCE_DIR}/MILPSolver.cpp
    ${PROJECT_SOURCE_DIR}/CPXMILPSolver.cpp
    )

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(CPLEX)
include_directories(${CPLEX_INCLUDE_DIR})

# TODO: Write and use a FindSMS++.cmake
include_directories(../SMS++/include)

include_directories("${PROJECT_INCLUDE_DIR}")
# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${MILPSolver_SOURCES})

add_subdirectory(test)
