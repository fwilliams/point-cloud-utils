cmake_minimum_required(VERSION 3.2)
project(py-sample-mesh)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

option(NPE_WITH_EIGEN "Whether to use the bundled version of Eigen" ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/external/numpyeigen/cmake)
include(numpyeigen)
option(NPE_WITH_EIGEN "Whether to use the bundled version of Eigen" ON)

#add_library(Eigen3::Eigen INTERFACE IMPORTED)
#set_target_properties(Eigen3::Eigen PROPERTIES
#  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib/eigenlib")

npe_add_module(py_sample_mesh
  BINDING_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample.cpp)
target_include_directories(py_sample_mesh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib)
#target_link_libraries(py_sample_mesh PUBLIC Eigen3::Eigen)
set_target_properties(py_sample_mesh PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -msse2")
set_target_properties(py_sample_mesh PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sample_mesh")
