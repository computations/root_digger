find_package(benchmark)
if(benchmark_FOUND)
    add_subdirectory(src)
endif()
