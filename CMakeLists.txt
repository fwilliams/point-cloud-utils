cmake_minimum_required(VERSION 3.2)
project(point-cloud-utils)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

option(GEOGRAM_WITH_HLBFGS "Non-linear solver (Yang Liu's HLBFGS)"             ON)
option(NPE_WITH_EIGEN      "Whether to use the bundled version of Eigen"       ON)
option(EIGEN_WITH_MKL      "Whether to build Eigen with intel MKL or not"      OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/numpyeigen/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)

# This will define Eigen3::Eigen if we enabled NPE_WITH_EIGEN
include(numpyeigen)

if (${EIGEN_WITH_MKL})
  find_package(MKL)
  target_link_libraries(Eigen3::Eigen INTERFACE ${MKL_LIBRARIES} -fopenmp)
  target_include_directories(Eigen3::Eigen INTERFACE ${MKL_INCLUDE_DIRS})
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
endif(${EIGEN_WITH_MKL})

include(geogram)

npe_add_module(point_cloud_utils
  BINDING_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/point_cloud_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lloyd.cpp)
target_include_directories(point_cloud_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib ${CMAKE_CURRENT_SOURCE_DIR}/src/)
target_include_directories(point_cloud_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann)
target_link_libraries(point_cloud_utils PRIVATE geogram Eigen3::Eigen -fopenmp)
set_target_properties(point_cloud_utils PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -msse3")


