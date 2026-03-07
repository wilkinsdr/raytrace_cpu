add_executable(caustic_imageplane caustic_imageplane.cpp)
target_link_libraries(caustic_imageplane imageplane raytracer ${cfitsio_libs})
