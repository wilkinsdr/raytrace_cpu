include_directories("raytracer")
add_subdirectory (raytracer)

include_directories("source_tracer")
add_subdirectory (source_tracer)

include_directories("mapper")
add_subdirectory (mapper)

include_directories("outflow")
add_subdirectory (outflow)

add_executable(trace_rays trace_rays.cpp)
target_link_libraries(trace_rays pointsource raytracer)

add_executable(trace_rays_imageplane trace_rays_imageplane.cpp)
target_link_libraries(trace_rays_imageplane imageplane raytracer)

add_executable(healpix_to_disc healpix_to_disc.cpp)
target_link_libraries(healpix_to_disc healpix_pointsource raytracer)

add_executable(pointsource_to_disc pointsource_to_disc.cpp)
target_link_libraries(pointsource_to_disc pointsource raytracer)

add_executable(disc_source_photonfrac disc_source_photonfrac.cpp)
target_link_libraries(disc_source_photonfrac pointsource raytracer)

add_executable(disc_source_photonfrac_r disc_source_photonfrac_r.cpp)
target_link_libraries(disc_source_photonfrac_r pointsource raytracer)

add_executable(disc_source_return_angdist disc_source_return_angdist.cpp)
target_link_libraries(disc_source_return_angdist pointsource raytracer)

add_executable(healpix_disc_source_photonfrac healpix_disc_source_photonfrac.cpp)
target_link_libraries(healpix_disc_source_photonfrac healpix_pointsource raytracer)

add_executable(emissivity emissivity.cpp)
target_link_libraries(emissivity pointsource raytracer)

add_executable(gsl_test gsl_test.cpp)
target_link_libraries(gsl_test GSL::gsl)

add_executable(pcyg_sei pcyg_sei.cpp)
target_link_libraries(pcyg_sei GSL::gsl)

add_executable(pcyg_rel pcyg_rel.cpp)
target_link_libraries(pcyg_rel GSL::gsl)

add_executable(disc_wind disc_wind.cpp)
target_link_libraries(disc_wind GSL::gsl)