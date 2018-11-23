aux_source_directories(SOURCES "" .)
vor_add_executable(opennl_LSCM ${SOURCES})
target_link_libraries(opennl_LSCM geogram)

set_target_properties(opennl_LSCM PROPERTIES FOLDER "GEOGRAM")
