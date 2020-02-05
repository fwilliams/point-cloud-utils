#include <Eigen/Core>

#include <vector>
#include <mutex>

#include <geogram/basic/geometry.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <npe.h>
#include <pybind11/stl.h>

#include "common.h"


namespace {

/*
 * Convert a (V, F) pair to a geogram mesh
 */
template <typename DerivedV, typename DerivedF>
void vf_to_geogram_mesh(const Eigen::MatrixBase<DerivedV> &V, const Eigen::MatrixBase<DerivedF> &F, GEO::Mesh &M) {
    M.clear();
    // Setup vertices
    M.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < (int) M.vertices.nb(); ++i) {
        GEO::vec3 &p = M.vertices.point(i);
        p[0] = V(i, 0);
        p[1] = V(i, 1);
        p[2] = (V.cols() == 2 ? 0 : V(i, 2));
    }
    // Setup faces
    if (F.cols() == 3) {
        M.facets.create_triangles((int) F.rows());
    } else if (F.cols() == 4) {
        M.facets.create_quads((int) F.rows());
    } else {
        throw std::runtime_error("Mesh face type not supported");
    }
    for (int c = 0; c < (int) M.facets.nb(); ++c) {
        for (int lv = 0; lv < F.cols(); ++lv) {
            M.facets.set_vertex(c, lv, F(c, lv));
        }
    }
    M.facets.connect();
}

/*
 * Convert a (V, T) pair to a geogram mesh
 */
template <typename DerivedV, typename DerivedF, typename DerivedT>
void vft_to_geogram_tet_mesh(const Eigen::MatrixBase<DerivedV> &V, const Eigen::MatrixBase<DerivedF> &F, const Eigen::MatrixBase<DerivedT> &T, GEO::Mesh &M) {
    M.clear();

    // Setup vertices
    M.vertices.create_vertices((int) V.rows());
    for (int i = 0; i < (int) M.vertices.nb(); ++i) {
        GEO::vec3 &p = M.vertices.point(i);
        p[0] = V(i, 0);
        p[1] = V(i, 1);
        p[2] = (V.cols() == 2 ? 0 : V(i, 2));
    }

    // Setup faces
    if (F.cols() == 3) {
        M.facets.create_triangles((int) F.rows());
    } else {
        throw std::runtime_error("Mesh face type not supported");
    }
    for (int c = 0; c < (int) M.facets.nb(); ++c) {
        for (int lv = 0; lv < F.cols(); ++lv) {
            M.facets.set_vertex(c, lv, F(c, lv));
        }
    }

    // Setup tets
    if (T.cols() == 4) {
        M.cells.create_tets((int) T.rows());
    } else {
        throw std::runtime_error("Mesh cell type not supported");
    }
    for (int c = 0; c < (int) M.cells.nb(); ++c) {
        for (int lv = 0; lv < T.cols(); ++lv) {
            M.cells.set_vertex(c, lv, T(c, lv));
        }
    }

    M.cells.connect();
    M.facets.connect();
}

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

// Put geogram log messages in the trash where they belong
class GeoTrashCan: public GEO::LoggerClient {
protected:
    void div(const std::string& log) override { /*pybind11::print(log);*/ }
    void out(const std::string& log) override { /*pybind11::print(log);*/ }
    void warn(const std::string& log) override { /*pybind11::print(log);*/ }
    void err(const std::string& log) override { /*pybind11::print(log);*/ }
    void status(const std::string& log) override { /*pybind11::print(log);*/ }
};


// Geogram is stateful and needs to be initialized.
// These variables keep track of whether geogram is initialized in a thread-safe manner.
// I'm using p-threads so none of this will work on Windows.
bool geogram_is_initialized = false;
std::mutex geogram_init_mutex;

// Initialize geogram exactly once.
void init_geogram_only_once() {
  std::lock_guard<std::mutex> guard(geogram_init_mutex);

  if (!geogram_is_initialized) {
    GEO::initialize();

    GEO::Logger *geo_logger = GEO::Logger::instance();
    geo_logger->unregister_all_clients();
    geo_logger->register_client(new GeoTrashCan());
    geo_logger->set_pretty(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    geogram_is_initialized = true;
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

npe_end_code()

