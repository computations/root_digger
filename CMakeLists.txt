cmake_minimum_required(VERSION 3.14)
project(root_digger)

set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)
add_subdirectory(lib)
include_directories(${LOCALDEPS_DIR}/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
add_subdirectory(src)
add_subdirectory(test)
