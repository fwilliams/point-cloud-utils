cmake_minimum_required(VERSION 3.2)
project(LibiglPythonTests)


file(GLOB PYTHON_TESTS
    "*.py"
)
#message("${PYTHON_TESTS}")
foreach(test ${PYTHON_TESTS})
get_filename_component(name ${test} NAME_WE)
add_test(NAME ${name} COMMAND ${PYTHON_EXECUTABLE} ${test} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endforeach()
