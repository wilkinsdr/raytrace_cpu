
add_executable(outflow outflow.cpp)
target_link_libraries(outflow sourcetracer_imageplane source_tracer raytracer)

add_executable(outflow_spectrum outflow_spectrum.cpp)
target_link_libraries(outflow_spectrum sourcetracer_imageplane source_tracer raytracer)

add_executable(outflow_emis_bin outflow_emis_bin.cpp)
target_link_libraries(outflow_emis_bin mapper_imageplane mapper raytracer)

add_executable(pcyg pcyg.cpp)

add_executable(pointsource_mapper pointsource_mapper.cpp)
target_link_libraries(pointsource_mapper mapper_pointsource mapper raytracer hdf5_cpp hdf5)

