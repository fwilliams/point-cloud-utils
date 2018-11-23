aux_source_directories(SOURCES "" .)
vor_add_executable(vorpacomp ${SOURCES})
target_link_libraries(vorpacomp geogram)

set_target_properties(vorpacomp PROPERTIES FOLDER "GEOGRAM")

