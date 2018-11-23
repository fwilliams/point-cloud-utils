vor_reset_warning_level()

aux_source_directories(SOURCES "" .)
include_directories(.)
vor_add_executable(mcc ${SOURCES})
target_link_libraries(mcc geogram)

set_target_properties(mcc PROPERTIES FOLDER "GEOGRAM")
