cmake_minimum_required(VERSION 3.19)
project(test_package)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()
#conan_basic_setup(KEEP_RPATHS)

add_executable(${PROJECT_NAME} test.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})
