add_executable(imageplane_disc_image imageplane_disc_image.cpp)
target_link_libraries(imageplane_disc_image imageplane raytracer ${cfitsio_libs})

