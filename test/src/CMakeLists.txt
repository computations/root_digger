set(CMAKE_BUILD_TYPE Debug)
add_executable(rd_test
    main.cpp
    data.cpp
    msa.cpp
    model.cpp
    tree.cpp
    ${RD_SOURCES}
)

set(TEST_LINK_LIBS ${PLLMODULES_LIBRARIES} ${PLLMODULES_LIBRARIES} Catch2::Catch2 pll_static)
target_compile_definitions(rd_test PRIVATE
    DATA_DIRECTORY_DNA_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../data/dna/
    DATA_DIRECTORY_TREE_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../data/tree/)

target_compile_options(rd_test PRIVATE -Wall -Wextra -pedantic
    -fno-omit-frame-pointer
    -fsanitize=address)
target_include_directories(rd_test PRIVATE ${CMAKE_SOURCE_DIR}/lib/libpll/src/)
target_include_directories(rd_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(rd_test ${TEST_LINK_LIBS})
