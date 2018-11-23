include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
aux_source_directories(SOURCES "lua" lua)
add_executable(geocod ${SOURCES})
target_link_libraries(geocod geogram_gfx geogram ${GLFW_LIBRARIES})
install_runtime_targets(geocod)

set_target_properties(geocod PROPERTIES FOLDER "GEOGRAM")
