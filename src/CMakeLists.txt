include_directories("raytracer")
add_subdirectory (raytracer)

add_executable(trace_rays trace_rays.cpp)
target_link_libraries(trace_rays pointsource raytracer)

add_executable(trace_rays_imageplane trace_rays_imageplane.cpp)
target_link_libraries(trace_rays_imageplane imageplane raytracer)

add_executable(outflow outflow.cpp)
target_link_libraries(outflow sourcetracer_imageplane sourcetracer raytracer)

add_executable(pcyg pcyg.cpp)

add_executable(pointsource_mapper pointsource_mapper.cpp raytracer/mapper_pointsource.cpp raytracer/mapper.cpp)
target_link_libraries(pointsource_mapper raytracer hdf5_cpp hdf5)

add_executable(healpix_to_disc healpix_to_disc.cpp)
target_link_libraries(healpix_to_disc healpix_pointsource raytracer)
