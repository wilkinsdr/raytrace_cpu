

add_library(raytracer raytracer.cpp)
add_library(pointsource pointsource.cpp)
add_library(imageplane imageplane.cpp)
add_library(healpix_pointsource healpix_pointsource.cpp stop_functs.cpp)

#if(OpenMP_CXX_FOUND)
#    add_library(raytracer raytracer.cpp)
#    target_link_libraries(raytracer PUBLIC OpenMP::OpenMP_CXX)
#
#    add_library(pointsource pointsource.cpp)
#    target_link_libraries(pointsource PUBLIC OpenMP::OpenMP_CXX)
#
#    add_library(imageplane imageplane.cpp)
#    target_link_libraries(imageplane PUBLIC OpenMP::OpenMP_CXX)
#
#    add_library(healpix_pointsource healpix_pointsource.cpp stop_functs.cpp)
#    target_link_libraries(healpix_pointsource PUBLIC OpenMP::OpenMP_CXX)
#endif()