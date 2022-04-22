add_executable(rd_bench
    EXCLUDE_FROM_ALL
    ${RD_SOURCES}
    main.cpp
    model.cpp
    tree.cpp
    ${CMAKE_SOURCE_DIR}/test/src/data.cpp
)

set(TEST_LINK_LIBS ${PLLMODULES_LIBRARIES} ${PLLMODULES_LIBRARIES}
    benchmark pll_static ${GSL_LIBS})
target_compile_definitions(rd_bench PRIVATE
    DATA_DIRECTORY_DNA_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../../test/data/dna/
    DATA_DIRECTORY_TREE_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../../test/data/tree/)
target_compile_options(rd_bench PRIVATE -Wall -Wextra -pedantic)
target_include_directories(rd_bench PRIVATE ${CMAKE_SOURCE_DIR}/lib/libpll/src/)
target_include_directories(rd_bench PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(rd_bench PRIVATE ${CMAKE_SOURCE_DIR}/test/src/)
target_link_libraries(rd_bench ${TEST_LINK_LIBS})
