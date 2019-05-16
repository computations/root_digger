set(CMAKE_BUILD_TYPE Debug)
add_executable(rd_test
    main.cpp
)

set(TEST_LINK_LIBS ${LINK_LIBS} Catch2::Catch2)

target_link_libraries(rd_test ${TEST_LINK_LIBS})
target_compile_options(rd_test PRIVATE -Wall -Wextra -pedantic)
