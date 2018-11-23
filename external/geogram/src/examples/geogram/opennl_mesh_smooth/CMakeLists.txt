aux_source_directories(SOURCES "" .)
vor_add_executable(opennl_mesh_smooth ${SOURCES})
target_link_libraries(opennl_mesh_smooth geogram)

set_target_properties(opennl_mesh_smooth PROPERTIES FOLDER "GEOGRAM")
