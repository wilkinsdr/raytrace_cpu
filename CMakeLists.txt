cmake_minimum_required (VERSION 2.6)
project (raytrace_cpu)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/obj)

include_directories("${CMAKE_SOURCE_DIR}/src")
add_subdirectory (src)

# HDF5
# HDF5_ROOT and HDF5_INCLUDE_DIRS need to come from CMakeLocal.txt
#set(HDF5_USE_STATIC_LIBRARIES ON)
#find_package(HDF5 REQUIRED CXX)
#set(hdf5 hdf5)
#set(hdf5_cpp hdf5_cpp)
#set(hdf5_libs hdf5_cpp hdf5)

