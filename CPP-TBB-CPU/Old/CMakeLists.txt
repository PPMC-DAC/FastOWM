cmake_minimum_required(VERSION 3.10)
project(owm LANGUAGES C CXX)

# Using C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # Required

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON) # Required

set(PRJ_SRC "${PROJECT_SOURCE_DIR}/src")
set(PRJ_INC "${PROJECT_SOURCE_DIR}/include")

# Flags
set(CMAKE_CXX_FLAGS "-fopenmp -O3")

set(CMAKE_C_FLAGS "-fopenmp -O3")

# add the src library
add_subdirectory(src)
add_subdirectory(application)
