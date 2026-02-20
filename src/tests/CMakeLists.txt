add_executable(gsl_test gsl_test.cpp)
target_link_libraries(gsl_test GSL::gsl)

add_executable(raytrace_rk4_test raytrace_rk4_test.cpp)
target_link_libraries(raytrace_rk4_test pointsource raytracer)
