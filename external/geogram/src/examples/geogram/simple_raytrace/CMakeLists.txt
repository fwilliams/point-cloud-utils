aux_source_directories(SOURCES "" .)
vor_add_executable(simple_raytrace ${SOURCES})
target_link_libraries(simple_raytrace geogram)

set_target_properties(simple_raytrace PROPERTIES FOLDER "GEOGRAM")
