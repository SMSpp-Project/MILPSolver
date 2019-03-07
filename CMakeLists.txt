project(MILPSolver DESCRIPTION "MILP Solver")
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)

# Local include directory
include_directories(include)

# TODO: Figure out how to install SMS++ library globally on macOS
include_directories("../sms_plus_plus_project/SMS++/include")
include_directories("../sms_plus_plus_project/SMILPBlock")

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories()

add_library(MILPSolver OBJECT src/MILPSolver.cpp)
