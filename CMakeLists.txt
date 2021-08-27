cmake_minimum_required(VERSION 3.10)
project(point-cloud-utils)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

option(GEOGRAM_WITH_HLBFGS "Non-linear solver (Yang Liu's HLBFGS)"             ON)
#option(NPE_WITH_EIGEN      "Whether to use the bundled version of Eigen"       ON)
option(EIGEN_WITH_MKL      "Whether to build Eigen with intel MKL or not"      OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Download NumpyEigen if we haven't already
set(EXTERNAL_DEP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
include(DownloadExternalDeps)
download_dep(numpyeigen
  GIT_REPOSITORY https://github.com/fwilliams/numpyeigen.git
  GIT_TAG 7252ee0d8ae51414bbca3e0445ecfc285930b285
)
list(APPEND CMAKE_MODULE_PATH ${EXTERNAL_DEP_DIR}/numpyeigen/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(numpyeigen) # This will define Eigen3::Eigen if we enabled NPE_WITH_EIGEN

# FIXME: MKL Support is totally broken
# Optionally build Eigen with the MKL backend
#if (${EIGEN_WITH_MKL})
#  find_package(MKL REQUIRED)
#  target_link_libraries(Eigen3::Eigen INTERFACE ${MKL_LIBRARIES})
#  target_include_directories(Eigen3::Eigen INTERFACE ${MKL_INCLUDE_DIRS})
#  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
#endif(${EIGEN_WITH_MKL})

# Embree
set(EMBREE_DIR "${EXTERNAL_DEP_DIR}/embree")
if(NOT TARGET embree)
  download_dep(embree
    GIT_REPOSITORY https://github.com/embree/embree.git
    GIT_TAG        v3.12.1
  )

  # Note: On macOS, building embree as a static lib can only be done with a single ISA target.
  set(EMBREE_MAX_ISA "DEFAULT" CACHE STRING "Selects highest ISA to support.")
  set(EMBREE_TESTING_INTENSITY 0 CACHE STRING "")
  set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL " ")
  set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE BOOL " ")
  set(EMBREE_TUTORIALS OFF CACHE BOOL " ")
  set(EMBREE_STATIC_LIB ON CACHE BOOL " ")
  if(MSVC)
    set(EMBREE_STATIC_RUNTIME ON CACHE BOOL "Use the static version of the C/C++ runtime library.")
  endif()

  add_subdirectory("${EMBREE_DIR}" "embree" EXCLUDE_FROM_ALL)
endif()
#compile_igl_module("embree")
#target_link_libraries(igl_embree INTERFACE embree)
#target_include_directories(igl_embree INTERFACE ${EMBREE_DIR}/include)
#target_compile_definitions(igl_embree INTERFACE -DEMBREE_STATIC_LIB)


# Build geogram
set(THIRD_PARTY_DIR ${EXTERNAL_DEP_DIR})
include(geogram)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET_TARGET_PROPERTIES(geogram PROPERTIES COMPILE_FLAGS -fopenmp LINK_FLAGS -fopenmp)
  target_compile_options(geogram PUBLIC    -fopenmp)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

  find_package(OpenMP REQUIRED)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
  endif()
endif()



# Create the python module
npe_add_module(_pcu_internal
  BINDING_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sample_point_cloud.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/point_cloud_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lloyd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/meshio.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/normals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/misc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/morton.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/remove_duplicates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/octree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/signed_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/closest_point_on_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/connected_components.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ray_mesh_intersection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smooth.cpp
  EXTRA_MODULE_FUNCTIONS
  hack_extra_bindings
  )
target_sources(_pcu_internal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/geogram_utils.cpp)
target_sources(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib/wrap/ply/plylib.cpp)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external)
target_include_directories(_pcu_internal PRIVATE ${EMBREE_DIR}/include)
target_link_libraries(_pcu_internal PRIVATE embree)
target_link_libraries(_pcu_internal PRIVATE geogram)
set_target_properties(_pcu_internal PROPERTIES COMPILE_FLAGS "-fvisibility=hidden -msse3")

if(MSVC)
  target_compile_options(_pcu_internal PRIVATE
    # Type conversion warnings. These can be fixed with some effort and possibly more verbose code.
    /wd4267 # conversion from 'size_t' to 'type', possible loss of data
    /wd4244 # conversion from 'type1' to 'type2', possible loss of data
    /wd4018 # signed/unsigned mismatch
    /wd4305 # truncation from 'double' to 'float'
    # This one is from template instantiations generated by autoexplicit.sh:
    /wd4667 # no function template defined that matches forced instantiation ()
    # This one is easy to fix, just need to switch to safe version of C functions
    /wd4996 # this function or variable may be unsafe
    # This one is when using bools in adjacency matrices
    /wd4804 #'+=': unsafe use of type 'bool' in operation
  )
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(_pcu_internal PUBLIC OpenMP::OpenMP_CXX)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries(_pcu_internal PUBLIC OpenMP::OpenMP_CXX)
endif()

enable_testing()
add_subdirectory(tests)

