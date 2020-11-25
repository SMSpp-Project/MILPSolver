cmake_minimum_required(VERSION 2.8.12)
project(test_package)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()
#conan_basic_setup(KEEP_RPATHS)

if (NOT BUILD_SHARED_LIBS)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_link_options("-Wl,-all_load")
    else()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

add_executable(${PROJECT_NAME} test.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})
