add_subdirectory(pll-modules)

set(PLLMOD_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/src")
set(LIBPLL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/pll-modules/libs/libpll/src")

file(COPY 
    ${PLLMOD_SRC}/pllmod_common.h
    ${PLLMOD_SRC}/msa/pll_msa.h
    ${PLLMOD_SRC}/tree/pll_tree.h
    ${PLLMOD_SRC}/util/pllmod_util.h
    ${LIBPLL_SRC}/pll.h
    DESTINATION ${LOCALDEPS_DIR}/include/libpll)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LOCALDEPS_DIR}/lib)

find_package(GSL)
if(GSL_FOUND)
    set(GSL_LIBS gsl gslcblas PARENT_SCOPE)
else()
    set(GSL_PREFIX "${CMAKE_BINARY_DIR}/localdeps/")
    file(MAKE_DIRECTORY ${GSL_PREFIX})
    include(ExternalProject)
    ExternalProject_Add(gsl
        URL ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
        URL_HASH SHA 9b3a32ba4f4e79fc6164fdc62e56003ae6e991e5
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gsl
        CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gsl/autogen.sh &&
        ${CMAKE_CURRENT_BINARY_DIR}/gsl/configure --prefix=${GSL_PREFIX}
        BUILD_COMMAND ${MAKE_INSTALL}
    )
    set(GSL_LIBS ${LOCALDEPS_DIR}/lib/libgsl.a
        ${LOCALDEPS_DIR}/lib/libgslcblas.a PARENT_SCOPE)
endif()
