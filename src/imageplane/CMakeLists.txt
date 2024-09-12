add_executable(imageplane_disc_image imageplane_disc_image.cpp ../raytracer/log_imageplane.cpp ../raytracer/log_imageplane.h)
target_link_libraries(imageplane_disc_image PUBLIC imageplane raytracer ${cfitsio_libs})

add_executable(imageplane_disc_image_geometry imageplane_disc_image_geometry.cpp ../raytracer/log_imageplane.cpp ../raytracer/log_imageplane.h)
target_link_libraries(imageplane_disc_image_geometry PUBLIC imageplane raytracer ${cfitsio_libs})