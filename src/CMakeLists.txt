add_executable(rd
    main.cpp
    model.cpp
    msa.cpp
    tree.cpp
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
