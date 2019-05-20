file(GLOB RD_SOURCES *.cpp)
list(REMOVE_ITEM RD_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

set(RD_SOURCES ${RD_SOURCES} PARENT_SCOPE)

add_executable(rd
    main.cpp
    ${RD_SOURCES}
)

set_target_properties(rd PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)

set(LINK_LIBS ${PLLMODULES_LIBRARIES} ${PLL_LIBRARIES})

target_link_libraries(rd 
    ${LINK_LIBS})
target_compile_options(rd PRIVATE -Wall -Wextra -Wshadow
    -Wdouble-promotion -Wmissing-include-dirs -Wtrampolines -pedantic)
