aux_source_directories(SOURCES "" .)
vor_add_executable(compute_OTM ${SOURCES})

if(GEOGRAM_WITH_VORPALINE)
  add_definitions(-DGEOGRAM_WITH_VORPALINE)
  target_link_libraries(compute_OTM exploragram geogram vorpalib)
else()
  target_link_libraries(compute_OTM exploragram geogram)
endif()

set_target_properties(compute_OTM PROPERTIES FOLDER "GEOGRAM")
