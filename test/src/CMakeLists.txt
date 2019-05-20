set(CMAKE_BUILD_TYPE Debug)
add_executable(rd_test
    main.cpp
    msa.cpp
    ${RD_SOURCES}
)

set(TEST_LINK_LIBS ${LINK_LIBS} Catch2::Catch2)

target_compile_options(rd_test PRIVATE -Wall -Wextra -pedantic
    -fno-omit-frame-pointer
    -fsanitize=address)
target_include_directories(rd_test PRIVATE ${CMAKE_SOURCE_DIR}/lib/libpll/src/)
target_include_directories(rd_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(rd_test ${TEST_LINK_LIBS})
