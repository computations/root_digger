set(CMAKE_BUILD_TYPE Debug)
add_executable(rd_test
    EXCLUDE_FROM_ALL
    main.cpp
    data.cpp
    msa.cpp
    model.cpp
    tree.cpp
    checkpoint.cpp
    util.cpp
    test_util.cpp
    ${RD_SOURCES}
)

set(TEST_LINK_LIBS corax_obj Catch2::Catch2 lbfgs_lib lapacke)
target_compile_definitions(rd_test PRIVATE
    DATA_DIRECTORY_DNA_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../data/dna/
    DATA_DIRECTORY_TREE_ABS=${CMAKE_CURRENT_SOURCE_DIR}/../data/tree/)
target_compile_options(rd_test PRIVATE -Wall -Wextra -pedantic
    -fno-omit-frame-pointer
    -fsanitize=address)
target_include_directories(rd_test PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/lib/lbfgsb/ ${CMAKE_SOURCE_DIR}/lib/coraxlib/src)
target_link_libraries(rd_test ${TEST_LINK_LIBS})
