cmake_minimum_required(VERSION 3.2)
project(point-cloud-utils)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

option(GEOGRAM_WITH_HLBFGS "Non-linear solver (Yang Liu's HLBFGS)"       ON)
option(NPE_WITH_EIGEN      "Whether to use the bundled version of Eigen" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/numpyeigen/cmake)

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
include(geogram)
include(numpyeigen)

npe_add_module(point_cloud_utils
  BINDING_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/point_cloud_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lloyd.cpp)

target_include_directories(point_cloud_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib ${CMAKE_CURRENT_SOURCE_DIR}/src/)
target_include_directories(point_cloud_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann)
target_link_libraries(point_cloud_utils PRIVATE geogram Eigen3::Eigen -fopenmp)
set_target_properties(point_cloud_utils PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -msse2")
#set_target_properties(sample_mesh PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sample_mesh")
