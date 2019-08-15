set(LBFGS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/lbfgsb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/linesearch.c
    ${CMAKE_CURRENT_SOURCE_DIR}/linpack.c
    ${CMAKE_CURRENT_SOURCE_DIR}/miniCBLAS.c
    ${CMAKE_CURRENT_SOURCE_DIR}/subalgorithms.c
)

add_library(lbfgs_lib
    ${LBFGS_SOURCES})
