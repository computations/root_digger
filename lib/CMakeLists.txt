add_subdirectory(pll-modules)

set(PLLMOD_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/src")
set(LIBPLL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/libs/libpll/src")

file(COPY 
    ${PLLMOD_SRC}/pllmod_common.h
    ${PLLMOD_SRC}/tree/pll_tree.h
    ${LIBPLL_SRC}/pll.h
    DESTINATION ${LOCALDEPS_DIR}/include/libpll)
