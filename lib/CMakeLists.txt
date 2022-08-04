set(CORAX_NONREV ON CACHE BOOL "Enable nonrev corax")
add_subdirectory(coraxlib EXCLUDE_FROM_ALL)

include_directories(AFTER ${LOCALDEPS_DIR}/include)

set(PLLMOD_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/src")
set(LIBPLL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/libs/libpll/src")

file(COPY 
    ${PLLMOD_SRC}/pllmod_common.h
    ${PLLMOD_SRC}/msa/pll_msa.h
    ${PLLMOD_SRC}/tree/pll_tree.h
    ${PLLMOD_SRC}/util/pllmod_util.h
    ${LIBPLL_SRC}/pll.h
    DESTINATION ${LOCALDEPS_DIR}/include/libpll)
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/lbfgsb/lbfgsb.h
    DESTINATION ${LOCALDEPS_DIR}/include)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LOCALDEPS_DIR}/lib)

add_subdirectory(pll-modules)
add_subdirectory(lbfgsb)
target_include_directories(rd PRIVATE lbfgsb)
