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

set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

# We duplicate the libararies here to solve a dumb link time error
set(LINK_LIBS ${PLLMODULES_LIBRARIES} ${PLL_LIBRARIES} ${PLLMODULES_LIBRARIES}
    ${GSL_LIBS} lbfgs_lib )

if(STATIC_BUILD)
    #target_link_options(rd PRIVATE -static -pie)
    #target_compile_options(rd PRIVATE -fPIE -static)
endif()

set(RD_COMPILE_DEFS GIT_REV=${RD_BUILD_VERSION} GIT_COMMIT=${RD_BUILD_COMMIT} BUILD_DATE=${RD_BUILD_DATE})

if(${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(RD_COMPILE_DEFS ${RD_COMPILE_DEFS} RD_DEBUG_FLAG)
endif()

if(MPI_BUILD)
    set(RD_COMPILE_DEFS ${RD_COMPILE_DEFS} MPI_BUILD)
    set(LINK_LIBS ${LINK_LIBS} MPI::MPI_CXX)
endif()

find_package(Threads)
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(rd PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(rd PRIVATE ${LINK_LIBS})
#target_link_options(rd PRIVATE -fno-omit-frame-pointer -fsanitize=address)
target_compile_options(rd PRIVATE -Wall -Wextra -Wshadow
    -Wdouble-promotion -Wmissing-include-dirs -Wtrampolines
    -pedantic -Wsign-conversion -Wnarrowing)
    #-fsanitize=address -fno-omit-frame-pointer)
target_compile_definitions(rd PRIVATE
    ${RD_COMPILE_DEFS})
