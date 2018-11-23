aux_source_directories(SOURCES "" .)
vor_add_executable(compute_delaunay ${SOURCES})
target_link_libraries(compute_delaunay geogram)

set_target_properties(compute_delaunay PROPERTIES FOLDER "GEOGRAM")
