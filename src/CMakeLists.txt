include_directories("raytracer")
add_subdirectory (raytracer)

add_executable(trace_rays trace_rays.cpp)
target_link_libraries(trace_rays pointsource raytracer)

add_executable(trace_rays_imageplane trace_rays_imageplane.cpp)
target_link_libraries(trace_rays_imageplane imageplane raytracer)

add_executable(outflow outflow.cpp)
target_link_libraries(outflow sourcetracer_imageplane sourcetracer raytracer)
