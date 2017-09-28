cmake_minimum_required (VERSION 2.6)
project (raytrace_cpu)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/obj)

include_directories("${CMAKE_SOURCE_DIR}/src")
add_subdirectory (src)
