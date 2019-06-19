cmake_minimum_required(VERSION 3.00)
project(root_digger)

set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g -fno-omit-frame-pointer" CACHE INTERNAL "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -fno-omit-frame-pointer" CACHE INTERNAL "")
set(CMAKE_C_FLAGS_DEBUG "-O3 -g -fno-omit-frame-pointer" CACHE INTERNAL "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -g -fno-omit-frame-pointer" CACHE INTERNAL "")

find_package(Git)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY "${local_dir}"
        OUTPUT_VARIABLE RD_BUILD_VERSION
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY "${local_dir}"
        OUTPUT_VARIABLE RD_BUILD_COMMIT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

string(TIMESTAMP RD_BUILD_DATE "%Y-%m-%dT%H:%M:%SZ")

set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)
add_subdirectory(lib)
include_directories(${LOCALDEPS_DIR}/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
add_subdirectory(src)
add_subdirectory(test)
