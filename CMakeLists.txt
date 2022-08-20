cmake_minimum_required(VERSION 3.00)
project(root_digger)

option(MPI_BUILD "Build with MPI" OFF)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer" CACHE INTERNAL "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE INTERNAL "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer " CACHE INTERNAL "")
set(CMAKE_C_FLAGS_RELEASE "-O3" CACHE INTERNAL "")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

find_package(Git)
if(GIT_FOUND)
    if(EXISTS ${CMAKE_SOURCE_DIR})
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags
            WORKING_DIRECTORY "${local_dir}"
            OUTPUT_VARIABLE RD_BUILD_VERSION
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (${RD_BUILD_VERSION})
            string(REGEX REPLACE "-g.+" "" RD_BUILD_VERSION ${RD_BUILD_VERSION})
            if(STATIC_BUILD)
                set(RD_BUILD_VERSION "${RD_BUILD_VERSION}-static")
            endif()
        endif()
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY "${local_dir}"
            OUTPUT_VARIABLE RD_BUILD_COMMIT
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if(NOT (RD_BUILD_VERSION AND RD_BUILD_COMMIT))
    set(RD_BUILD_VERSION "1.7.0")
    set(RD_BUILD_COMMIT "From Release")
endif()

if(MPI_BUILD)
    find_package(MPI)
    set(RD_BUILD_VERSION "${RD_BUILD_VERSION}-mpi")
endif()

string(TIMESTAMP RD_BUILD_DATE "%Y-%m-%d %H:%M:%S")

add_subdirectory(src)
set(LOCALDEPS_DIR ${PROJECT_BINARY_DIR}/localdeps)
add_subdirectory(lib)
include_directories(${LOCALDEPS_DIR}/include)


find_library(ASAN_FOUND asan)
if(ASAN_FOUND AND NOT MPI_BUILD)
    add_subdirectory(test)
    add_subdirectory(benchmark)
endif()
