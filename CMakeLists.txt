cmake_minimum_required(VERSION 3.2)
project(py-sample-mesh)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

option(NPE_WITH_EIGEN "Whether to use the bundled version of Eigen" ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/external/numpyeigen/cmake)
include(numpyeigen)
option(NPE_WITH_EIGEN "Whether to use the bundled version of Eigen" ON)

npe_add_module(sample_mesh
  BINDING_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/nanoflann.cpp)
target_include_directories(sample_mesh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib ${CMAKE_CURRENT_SOURCE_DIR}/src/)
set_target_properties(sample_mesh PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -msse2")
#set_target_properties(sample_mesh PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sample_mesh")
