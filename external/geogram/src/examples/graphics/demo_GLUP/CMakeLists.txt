include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(geogram_demo_GLUP ${SOURCES})
target_link_libraries(geogram_demo_GLUP geogram_gfx geogram ${GLFW_LIBRARIES})
install_runtime_targets(geogram_demo_GLUP)

set_target_properties(geogram_demo_GLUP PROPERTIES FOLDER "GEOGRAM")
