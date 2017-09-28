include_directories("common")
add_subdirectory (common)

add_executable(trace_rays trace_rays.cpp)
target_link_libraries(trace_rays pointsource raytracer)