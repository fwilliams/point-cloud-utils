aux_source_directories(SOURCES "" .)
vor_add_executable(manifold_harmonics ${SOURCES})
target_link_libraries(manifold_harmonics geogram)

set_target_properties(manifold_harmonics PROPERTIES FOLDER "GEOGRAM")
