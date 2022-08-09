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
  GIT_TAG 4916d926aa2b939bd8f625c7537563a1575dafe9
)
list(APPEND CMAKE_MODULE_PATH ${EXTERNAL_DEP_DIR}/numpyeigen/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(numpyeigen) # This will define Eigen3::Eigen if we enabled NPE_WITH_EIGEN

download_dep(tinyply
        GIT_REPOSITORY https://github.com/fwilliams/tinyply.git
        GIT_TAG point_cloud_utils
)

# Build manifold as a static lib
download_dep(manifold
  GIT_REPOSITORY https://github.com/hjwdzh/Manifold.git
  GIT_TAG 81fd342e578e29fc15fb75d2b4d1e3c821fe33cb
)
set(MANIFOLD_SRC_DIR ${EXTERNAL_DEP_DIR}/manifold/src)
set(
    manifold_SRC
    ${MANIFOLD_SRC_DIR}/BVH.cpp
    ${MANIFOLD_SRC_DIR}/BVH.h
    ${MANIFOLD_SRC_DIR}/Intersection.cpp
    ${MANIFOLD_SRC_DIR}/Intersection.h
    ${MANIFOLD_SRC_DIR}/main.cpp
    ${MANIFOLD_SRC_DIR}/Model_OBJ.cpp
    ${MANIFOLD_SRC_DIR}/Model_OBJ.h
    ${MANIFOLD_SRC_DIR}/Octree.h)
add_library(manifold STATIC ${manifold_SRC})
set_target_properties(manifold PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(manifold PUBLIC ${MANIFOLD_SRC_DIR})
target_include_directories(manifold PRIVATE ${EXTERNAL_DEP_DIR}/manifold/3rd)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(manifold PRIVATE "${CMAKE_CXX_FLAGS} -Wall -pthread -DWITH_OMP  -Wno-int-in-bool-context -Wsign-compare -fsanitize=address")
else()
    target_compile_options(manifold PRIVATE "${CMAKE_CXX_FLAGS} -DWITH_OMP ")
endif()



# set(
#     mls_utils_SRC
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/mls_utils/balltree.cpp
# )
# add_library(mls_utils STATIC ${mls_utils_SRC})
# set_target_properties(mls_utils PROPERTIES LINKER_LANGUAGE CXX)
# target_include_directories(mls_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external)
# target_include_directories(mls_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib)
# target_include_directories(mls_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
# target_include_directories(mls_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/numpyeigen/external/eigen)



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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_normals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/point_cloud_normals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/misc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/morton.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/remove_duplicates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/octree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/signed_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/closest_point_on_mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/connected_components.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ray_mesh_intersection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ray_point_cloud_intersection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smooth.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/manifold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/curvature.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/sparse_voxel_grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/marching_cubes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_decimate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/remove_unreferenced_mesh_vertices.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/face_areas.cpp
  EXTRA_MODULE_FUNCTIONS
  hack_extra_bindings
  hack_extra_ray_mesh_bindings
  )

# TODO: Make common into its own library that we link statically
target_sources(_pcu_internal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common/morton_code.cpp)
target_sources(_pcu_internal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common/geogram_utils.cpp)
target_sources(_pcu_internal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common/embree_intersector.cpp)
target_sources(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib/wrap/ply/plylib.cpp)
# target_sources(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyply/source/tinyply.cpp)
target_include_directories(_pcu_internal PRIVATE ${EXTERNAL_DEP_DIR}/tinyply/source)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcglib)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/nanoflann)
target_include_directories(_pcu_internal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external)
target_include_directories(_pcu_internal PRIVATE ${EMBREE_DIR}/include)
target_link_libraries(_pcu_internal PRIVATE embree)
target_link_libraries(_pcu_internal PRIVATE geogram)
target_link_libraries(_pcu_internal PRIVATE manifold)
# target_link_libraries(_pcu_internal PRIVATE mls_utils)
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

