include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(GLUP_basic_example ${SOURCES})

target_link_libraries(GLUP_basic_example geogram_gfx geogram ${GLFW_LIBRARIES})

set_target_properties(GLUP_basic_example PROPERTIES FOLDER "GEOGRAM")
