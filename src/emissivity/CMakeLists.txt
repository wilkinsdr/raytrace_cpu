add_executable(emissivity emissivity.cpp)
target_link_libraries(emissivity pointsource raytracer)

add_executable(emissivity_rd emissivity_rd.cpp)
target_link_libraries(emissivity_rd pointsource raytracer)
