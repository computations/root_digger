set(CORAX_NONREV ON CACHE BOOL "Enable nonrev corax")
add_subdirectory(coraxlib EXCLUDE_FROM_ALL)
add_subdirectory(lbfgsb EXCLUDE_FROM_ALL)

target_include_directories(rd PRIVATE lbfgsb)
