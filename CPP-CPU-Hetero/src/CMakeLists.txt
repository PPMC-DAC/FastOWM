
# Include all the files in src/ and include/
file(GLOB_RECURSE LIB_CODE
        "${PRJ_SRC}/envi.cpp"
        "${PRJ_SRC}/environment.c"
        "${PRJ_INC}/environment.h"
        "${PRJ_INC}/envi.h" )

add_library(src SHARED ${LIB_CODE})

# Declare include directories
target_include_directories(src PUBLIC $<BUILD_INTERFACE:${PRJ_INC}>)

# Linked libraries
target_link_libraries(src PUBLIC tbb rt m)
