include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(vorpaview ${SOURCES})
target_link_libraries(vorpaview geogram_gfx geogram ${GLFW_LIBRARIES})
install_runtime_targets(vorpaview)

set_target_properties(vorpaview PROPERTIES FOLDER "GEOGRAM")
