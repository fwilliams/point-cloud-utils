#include <Eigen/Core>

#include <vector>
#include <mutex>
#include <sstream>
#include <chrono>


#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/mesh_preprocessing.h>
#include <geogram/mesh/mesh_degree3_vertices.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/numerics/predicates.h>

#include <npe.h>
#include <pybind11/stl.h>

#include "common.h"
#include "geogram_utils.h"


namespace {

/*
 * Sample a mesh using lloyd relaxation
 */
template <typename DerivedV, typename DerivedF, typename DerivedP>
void sample_tri_mesh_lloyd(const Eigen::MatrixBase<DerivedV> &V,
                       const Eigen::MatrixBase<DerivedF> &F,
                       int num_samples, int num_lloyd, int num_newton,
                       Eigen::PlainObjectBase<DerivedP> &P) {
    assert(num_samples > 3);
    GEO::Mesh M;
    vf_to_geogram_mesh(V, F, M);

    GEO::CentroidalVoronoiTesselation CVT(&M);
    bool was_quiet = GEO::Logger::instance()->is_quiet();
    GEO::Logger::instance()->set_quiet(false);
    CVT.compute_initial_sampling(num_samples);
    GEO::Logger::instance()->set_quiet(was_quiet);

    if (num_lloyd > 0) {
        CVT.Lloyd_iterations(num_lloyd);
    }

    if (num_newton > 0) {
        CVT.Newton_iterations(num_newton);
    }

    P.resize(3, num_samples);
    std::copy_n(CVT.embedding(0), 3*num_samples, P.data());
    P.transposeInPlace();
}


/*
 * Sample a mesh using lloyd relaxation
 */
template <typename DerivedV, typename DerivedF, typename DerivedT, typename DerivedP>
void sample_tet_mesh_lloyd(const Eigen::MatrixBase<DerivedV> &V,
                          const Eigen::MatrixBase<DerivedF> &F,
                          const Eigen::MatrixBase<DerivedT> &T,
                          int num_samples, int num_lloyd, int num_newton,
                          Eigen::PlainObjectBase<DerivedP> &P) {
    assert(num_samples > 3);
    GEO::Mesh M;
    vft_to_geogram_tet_mesh(V, F, T, M);

    GEO::CentroidalVoronoiTesselation CVT(&M);
    CVT.set_volumetric(true);
    bool was_quiet = GEO::Logger::instance()->is_quiet();
    GEO::Logger::instance()->set_quiet(true);
    CVT.compute_initial_sampling(num_samples);
    GEO::Logger::instance()->set_quiet(was_quiet);

    if (num_lloyd > 0) {
        CVT.Lloyd_iterations(num_lloyd);
    }

    if (num_newton > 0) {
        CVT.Newton_iterations(num_newton);
    }

    P.resize(3, num_samples);
    std::copy_n(CVT.embedding(0), 3*num_samples, P.data());
    P.transposeInPlace();
}

/*
 * Sample a mesh using lloyd relaxation
 */
template <typename DerivedV, typename DerivedF, typename DerivedT, typename DerivedP>
void voronoi_centroids(const Eigen::MatrixBase<DerivedV> &V,
                     const Eigen::MatrixBase<DerivedF> &F,
                     const Eigen::MatrixBase<DerivedT> &T,
                     const Eigen::MatrixBase<DerivedP> &centers,
                     Eigen::PlainObjectBase<DerivedP> &P) {

    std::vector<double> centers_vec;
    for (GEO::index_t i = 0; i < centers.rows(); i++) {
        for (GEO::index_t j = 0; j < centers.cols(); j++) {
            centers_vec.push_back(centers(i, j));
        }
    }
    GEO::Mesh M;
    vft_to_geogram_tet_mesh(V, F, T, M);

    GEO::coord_index_t dimension_ = GEO::coord_index_t(M.vertices.dimension());
    const std::string& delaunay = "default";
    GEO::Delaunay_var delaunay_ = GEO::Delaunay::create(dimension_, delaunay);
    GEO::RestrictedVoronoiDiagram_var RVD_ = GEO::RestrictedVoronoiDiagram::create(delaunay_, &M);
    RVD_->set_volumetric(true);

    std::vector<double> mg;
    std::vector<double> m;
    RVD_->set_check_SR(false);

    P.resize(centers.rows(), centers.cols());

    GEO::index_t nb_points = centers.rows();
    mg.assign(nb_points * dimension_, 0.0);
    m.assign(nb_points, 0.0);
    delaunay_->set_vertices(nb_points, centers_vec.data());
    RVD_->compute_centroids(mg.data(), m.data());
    GEO::index_t cur = 0;
    for(GEO::index_t j = 0; j < nb_points; j++) {
        if(m[j] > 1e-30) {
            double s = 1.0 / m[j];
            for(GEO::index_t coord = 0; coord < dimension_; coord++) {
                P(j, coord) = s * mg[cur + coord];
            }
        }
        cur += dimension_;
    }
}


} // namespace





const char* lloyd_2d_doc = R"Qu8mg5v7(
Generate n samples in the unit square, [0, 1]^2 using Lloyd's algorithm
(https://en.wikipedia.org/wiki/Lloyd%27s_algorithm).

Parameters
----------
n : The number of 2d point samples to generate
num_lloyd : The number of Lloyd iterations to do (default 10)
num_newton : The number of Newton iterations to do when computing Voronoi diagrams (default 10)

Returns
-------
A n by 2 array of point samples in the unit square [0, 1]^2

)Qu8mg5v7";
npe_function(lloyd_2d)
npe_arg(n, int)
npe_default_arg(num_lloyd, int, 10)
npe_default_arg(num_newton, int, 10)
npe_doc(lloyd_2d_doc)
npe_begin_code()
{
  if (n <= 0) {
    throw pybind11::value_error("n must be a positive integer. Got n = " + std::to_string(n));
  }

  init_geogram_only_once();

  Eigen::MatrixXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       1.0, 1.0, 0.0,
       0.0, 1.0, 0.0;
  Eigen::MatrixXi F(2, 3);
  F << 0, 1, 3,
       1, 2, 3;

  Eigen::MatrixXd P;
  sample_tri_mesh_lloyd(V, F, n, num_lloyd, num_newton, P);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> P2d;
  P2d = P.block(0, 0, P.rows(), 2);
  return npe::move(P2d);
}
npe_end_code()



const char* lloyd_3d_doc = R"Qu8mg5v7(
Generate n samples in the unit cube, [0, 1]^3 using Lloyd's algorithm
(https://en.wikipedia.org/wiki/Lloyd%27s_algorithm).

Parameters
----------
n : The number of 3d point samples to generate
num_lloyd : The number of Lloyd iterations to do (default 10)
num_newton : The number of Newton iterations to do when computing Voronoi diagrams (default 10)

Returns
-------
A n by 3 array of point samples in the unit square [0, 1]^3

)Qu8mg5v7";
npe_function(lloyd_3d)
npe_arg(n, int)
npe_default_arg(num_lloyd, int, 10)
npe_default_arg(num_newton, int, 10)
npe_doc(lloyd_3d_doc)
npe_begin_code()
{
  if (n <= 0) {
    throw pybind11::value_error("n must be a positive integer. Got n = " + std::to_string(n));
  }

  init_geogram_only_once();

  Eigen::MatrixXd V(9, 3);
  V << 0. , 1. , 0. ,
       0. , 0. , 1. ,
       1. , 0. , 1. ,
       1. , 0. , 0. ,
       0. , 0. , 0. ,
       1. , 1. , 1. ,
       1. , 1. , 0. ,
       0. , 1. , 1. ,
       0.5, 0.5, 0.5;

  Eigen::MatrixXi T(12, 4);
  T << 8, 0, 7, 5,
       6, 2, 5, 8,
       3, 1, 2, 8,
       0, 8, 7, 4,
       7, 8, 5, 1,
       8, 0, 5, 6,
       3, 1, 8, 4,
       2, 1, 5, 8,
       8, 0, 6, 4,
       7, 8, 1, 4,
       6, 4, 3, 8,
       2, 6, 3, 8;
  Eigen::MatrixXi F(12, 3);
  F << 0, 5, 7,
       6, 2, 5,
       3, 1, 2,
       0, 7, 4,
       7, 5, 1,
       0, 6, 5,
       3, 4, 1,
       2, 1, 5,
       0, 4, 6,
       7, 1, 4,
       6, 4, 3,
       2, 6, 3;


  Eigen::MatrixXd P;
  sample_tet_mesh_lloyd(V, F, T, n, num_lloyd, num_newton, P);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> P3d;
  P3d = P;
  return npe::move(P3d);
}
npe_end_code()


const char* voronoi_centroids_unit_cube_doc = R"Qu8mg5v7(
Compute the centroids of a 3d voronoi diagram restricted to the unit cube [0, 1]^3

Parameters
----------
centers: A [n, 3] array where each row is a point in a voronoi diagram restricted to
         the unit cube [0, 1]^3

Returns
-------
A [n, 3] array of centroids where each row is the centroid of the corresponding
voronoi cell in the input centers.

)Qu8mg5v7";
npe_function(voronoi_centroids_unit_cube)
npe_arg(centers, dense_float, dense_double)
npe_doc(lloyd_3d_doc)
npe_begin_code()
{
  if (centers.maxCoeff() > 1.0 || centers.minCoeff() < 0.0) {
    throw pybind11::value_error("Centers of voronoi diagram must lie in the unit cube [0, 1]^3");
  }

  init_geogram_only_once();

  Eigen::MatrixXd V(9, 3);
  V << 0. , 1. , 0. ,
       0. , 0. , 1. ,
       1. , 0. , 1. ,
       1. , 0. , 0. ,
       0. , 0. , 0. ,
       1. , 1. , 1. ,
       1. , 1. , 0. ,
       0. , 1. , 1. ,
       0.5, 0.5, 0.5;

  Eigen::MatrixXi T(12, 4);
  T << 8, 0, 7, 5,
       6, 2, 5, 8,
       3, 1, 2, 8,
       0, 8, 7, 4,
       7, 8, 5, 1,
       8, 0, 5, 6,
       3, 1, 8, 4,
       2, 1, 5, 8,
       8, 0, 6, 4,
       7, 8, 1, 4,
       6, 4, 3, 8,
       2, 6, 3, 8;
  Eigen::MatrixXi F(12, 3);
  F << 0, 5, 7,
       6, 2, 5,
       3, 1, 2,
       0, 7, 4,
       7, 5, 1,
       0, 6, 5,
       3, 4, 1,
       2, 1, 5,
       0, 4, 6,
       7, 1, 4,
       6, 4, 3,
       2, 6, 3;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ctrs_copy = centers.template cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> P3d;
  voronoi_centroids(V, F, T, ctrs_copy, P3d);
  return npe::move(P3d);
}
npe_end_code()


const char* sample_mesh_lloyd_doc = R"Qu8mg5v7(
Generate n samples on a surface defined by a triangle mesh using
Lloyd's algorithm (https://en.wikipedia.org/wiki/Lloyd%27s_algorithm).

Parameters
----------
v : A #v by 3 array where each row is a vertex of the input mesh
f : A #f by 3 array of indices into `v` where each row is a triangle of the input mesh.
n : The number of surface samples to generate
num_lloyd : The number of Lloyd iterations to do (default 10)
num_newton : The number of Newton iterations to do when computing Voronoi diagrams (default 10)

Returns
-------
A n by 3 array of point samples on the input surface defined by (v, f)

)Qu8mg5v7";

npe_function(sample_mesh_lloyd)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_long, dense_longlong, dense_uint, dense_ulong, dense_ulonglong)
npe_arg(num_samples, int)
npe_default_arg(num_lloyd, int, 10)
npe_default_arg(num_newton, int, 10)
npe_default_arg(return_mesh, bool, false)
npe_doc(sample_mesh_lloyd_doc)
npe_begin_code()
{
  validate_mesh(v, f);
  if (num_samples <= 0) {
    throw pybind11::value_error("num_samples must be > 0");
  }

  init_geogram_only_once();

  Eigen::MatrixXd V = v.template cast<double>();
  Eigen::MatrixXi F = f.template cast<int>();

  GEO::Mesh M;
  vf_to_geogram_mesh(V, F, M);

  GEO::CentroidalVoronoiTesselation CVT(&M);
  bool was_quiet = GEO::Logger::instance()->is_quiet();
  GEO::Logger::instance()->set_quiet(true);
  CVT.compute_initial_sampling(num_samples);
  GEO::Logger::instance()->set_quiet(was_quiet);

  if (num_lloyd > 0) {
      CVT.Lloyd_iterations(num_lloyd);
  }

  if (num_newton > 0) {
      CVT.Newton_iterations(num_newton);
  }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> P;
  P.resize(num_samples, 3);
  std::copy_n(CVT.embedding(0), 3*num_samples, P.data());
  return npe::move(P);
}
npe_end_code()


