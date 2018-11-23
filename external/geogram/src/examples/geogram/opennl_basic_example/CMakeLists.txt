aux_source_directories(SOURCES "" .)
vor_add_executable(opennl_basic_example ${SOURCES})
target_link_libraries(opennl_basic_example geogram)

set_target_properties(opennl_basic_example PROPERTIES FOLDER "GEOGRAM")


