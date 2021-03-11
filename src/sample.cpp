#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>

#include <fstream>
#include <iostream>
#include <functional>

#include "common.h"


namespace {

using namespace vcg;
class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};
class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3d, vertex::Normal3d, vertex::BitFlags> {};
class MyFace    : public Face<MyUsedTypes, face::FFAdj,  face::Normal3d, face::VertexRef, face::BitFlags> {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>> {};


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF, typename DerivedN>
static void vcg_mesh_from_vfn(const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedF>& F, const Eigen::MatrixBase<DerivedN>& N, MyMesh& m) {

  MyMesh::VertexIterator vit = Allocator<MyMesh>::AddVertices(m, V.rows());
  std::vector<MyMesh::VertexPointer> ivp(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    ivp[i] = &*vit;
    vit->P() = MyMesh::CoordType(V(i, 0), V(i, 1), V(i, 2));
    if (N.rows() > 0) {
      vit->N() = MyMesh::CoordType(N(i, 0), N(i, 1), N(i, 2));
    }
    vit++;
  }

  if (F.rows() > 0) {
    MyMesh::FaceIterator fit = Allocator<MyMesh>::AddFaces(m, F.rows());
    for (int i = 0; i < F.rows(); i++) {
      fit->V(0) = ivp[F(i, 0)];
      fit->V(1) = ivp[F(i, 1)];
      fit->V(2) = ivp[F(i, 2)];
      fit++;
    }
  }

  tri::UpdateBounding<MyMesh>::Box(m);
}


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF>
static void vcg_mesh_from_vf(const Eigen::MatrixBase<DerivedV>& V, const Eigen::MatrixBase<DerivedF>& F,MyMesh& m) {
    Eigen::MatrixXd N(0, 3);
    vcg_mesh_from_vfn(V, F, N, m);
}


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV>
static void vcg_mesh_from_v(const Eigen::MatrixBase<DerivedV>& V, MyMesh& m) {
    Eigen::MatrixXi F(0, 3);
    Eigen::MatrixXd N(0, 3);
    vcg_mesh_from_vfn(V, F, N, m);
}


/*
 * Use this to sample vertex indices when we are sampling from a point cloud
 */
template <class MeshType>
class EigenVertexIndexSampler
{
public:
    typedef typename MeshType::VertexType  VertexType;

    // The mesh we are sampling from
    MeshType& sampled_mesh;

    // Indices into the mesh vertex array, this is an eigen matrix of some type
    typedef Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, 1> IndexArray;
    IndexArray& indices;

    // Number of vertices
    int vcount = 0;

    EigenVertexIndexSampler(MeshType& in_mesh, IndexArray& out_inds) :
        sampled_mesh(in_mesh), indices(out_inds) {
    }

    void trim() {
        indices.conservativeResize(vcount, 1);
    }

    void reset() {
        vcount = 0;
    }

    void maybe_resize()  {
        // If we are about to overflow indexes, double its size
        if (indices.size() <= vcount) {
            const int n_rows = indices.size() == 0 ? 1024 : indices.size();
            indices.conservativeResize(2*n_rows, 1);
        }
    }

    void AddVert(const VertexType &p) {
        maybe_resize();
        std::ptrdiff_t p_offset = &p - &*sampled_mesh.vert.begin();
        indices(vcount, 0) = p_offset;
        vcount += 1;
    }
}; // end class EigenVertexIndexSampler


/*
 * Use this to sample vertex indices when we are sampling from a point cloud
 */
template <class MeshType, class BCType>
class EigenBarycentricSampler
{
public:
    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::CoordType   CoordType;

    // The mesh we are sampling from
    MeshType* sampled_mesh;
    MeshType* output_mesh;

    // Indices into the mesh vertex array, this is an eigen matrix of some type
    typedef Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, 1> IndexArray;
    IndexArray& indices;

    // Barycentric coordinates of each sampled vertex
    BCType& barycentric_coords;

    // Number of vertices
    int vcount = 0;

    bool perFaceNormal;  // default false; if true the sample normal is the face normal, otherwise it is interpolated

    EigenBarycentricSampler(MeshType* in_mesh, MeshType* output_mesh, BCType& out_bc, IndexArray& out_face_idxs) :
        sampled_mesh(in_mesh), output_mesh(output_mesh), barycentric_coords(out_bc), indices(out_face_idxs), perFaceNormal(false) {
    }

    void trim() {
        barycentric_coords.conservativeResize(vcount, 3);
        indices.conservativeResize(vcount, 1);
    }

    void reset() {
        vcount = 0;
    }

    void maybe_resize()  {
        // If we are about to overflow indexes, double its size
        if (barycentric_coords.rows() <= vcount) {
            const int n_rows = barycentric_coords.rows() == 0 ? 1024 : barycentric_coords.rows();
            barycentric_coords.conservativeResize(2 * n_rows, 3);
        }
        if (indices.size() <= vcount) {
            const int n_rows = indices.size() == 0 ? 1024 : indices.size();
            barycentric_coords.conservativeResize(2 * n_rows, 1);
        }
    }

    void AddFace(const FaceType &f, CoordType p) {
        maybe_resize();
        std::ptrdiff_t f_offset = &f - &*(sampled_mesh->face.begin());
        indices(vcount, 0) = f_offset;
        barycentric_coords(vcount, 0) = p[0];
        barycentric_coords(vcount, 1) = p[1];
        barycentric_coords(vcount, 2) = p[2];
        vcount += 1;

        if (output_mesh != nullptr) {
            tri::Allocator<MeshType>::AddVertices(*output_mesh, 1);
            output_mesh->vert.back().P() = f.cP(0) * p[0] + f.cP(1) * p[1] +f.cP(2) * p[2];
            if(perFaceNormal) {
                output_mesh->vert.back().N() = f.cN();
            } else {
                output_mesh->vert.back().N() = f.cV(0)->N() * p[0] + f.cV(1)->N() * p[1] + f.cV(2)->N() * p[2];
            }
            if(tri::HasPerVertexQuality(*output_mesh)) {
               output_mesh->vert.back().Q() = f.cV(0)->Q() * p[0] + f.cV(1)->Q() * p[1] + f.cV(2)->Q() * p[2];
            }
        }
    }
}; // end class EigenVertexIndexSampler


class AccumulatedPoint {
public:
    AccumulatedPoint()
            : num_of_points_(0),
              point_(0.0, 0.0, 0.0),
              normal_(0.0, 0.0, 0.0),
              color_(0.0, 0.0, 0.0),
              alpha_(0.0) {}

public:
    void AddPoint(Eigen::Vector3d point, Eigen::Vector3d normal, Eigen::Vector3d color, double alpha) {
        point_ += point;
        normal_ += normal;
        color_ += color;
        alpha_ += alpha;
        num_of_points_ += 1;
    }

    Eigen::Vector3d GetAveragePoint() const {
        return point_ / double(num_of_points_);
    }

    Eigen::Vector3d GetAverageNormal() const {
        // Call NormalizeNormals() afterwards if necessary
        return normal_ / double(num_of_points_);
    }

    Eigen::Vector3d GetAverageColor() const {
        return color_ / double(num_of_points_);
    }

    double GetAverageAlpha() const {
        return alpha_ / double(num_of_points_);
    }

    int GetNumPoints() const {
        return num_of_points_;
    }

public:
    int num_of_points_;
    Eigen::Vector3d point_;
    Eigen::Vector3d normal_;
    Eigen::Vector3d color_;
    double alpha_;
};


template <typename T>
struct hash_eigen {
    std::size_t operator()(T const& matrix) const {
        size_t seed = 0;
        for (int i = 0; i < (int)matrix.size(); i++) {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <typename DerivedV, typename DerivedN, typename DerivedC,
          typename DerivedOutV, typename DerivedOutN, typename DerivedOutC>
void downsample_point_cloud_to_voxels(const DerivedV& V,
                                      const DerivedN& N,
                                      const DerivedC& C,
                                      Eigen::Vector3d voxel_size3,
                                      Eigen::Vector3d voxel_min_bound,
                                      Eigen::Vector3d voxel_max_bound,
                                      DerivedOutV& outV,
                                      DerivedOutN& outN,
                                      DerivedOutC& outC,
                                      int min_pts_per_bin) {
    // The way this is templated is dangerous since you could pass in an Eigen expression template. However,
    // we always pass either an Eigen::Map or Eigen::Matrix so it should be okay.

    for (int i = 0; i < 3; i++) {
        if (voxel_size3[i] <= 0.0) {
            throw pybind11::value_error("Voxel size is negative");
        }
        if (voxel_size3[i] * std::numeric_limits<int>::max() < (voxel_max_bound - voxel_min_bound)[i]) {
            throw pybind11::value_error("Voxel size is too small");
        }
    }

    bool has_normals = (N.rows() != 0 && N.cols() != 0);
    bool has_colors = (C.rows() != 0 && C.cols() != 0);
    if (has_normals && N.rows() != V.rows()) {
        throw pybind11::value_error("Invalid number of normals (" + std::to_string(N.rows()) +
                                    "). Must match number of input vertices (" +
                                    std::to_string(V.rows()) + ") or b 0.");
    }
    if (has_normals && N.cols() != 3) {
        throw pybind11::value_error("Invalid dimension for normals. Got shape (" + std::to_string(N.rows()) +
                                    ", " + std::to_string(N.cols()) + "). Must have shape (N, 3)");
    }
    if (has_colors && C.rows() != V.rows()) {
        throw pybind11::value_error("Invalid number of colors (" + std::to_string(C.rows()) +
                                    "). Must match number of input vertices (" +
                                    std::to_string(V.rows()) + ") or b 0.");
    }
    if (has_colors && (C.cols() != 3 && C.cols() != 4)) {
        throw pybind11::value_error("Invalid dimension for colors. Got shape (" + std::to_string(C.rows()) +
                                    ", " + std::to_string(C.cols()) + "). Must have shape (N, 3) or (N, 4)");
    }
    bool color_has_alpha_channel = C.cols() == 4;

    std::unordered_map<Eigen::Vector3i, AccumulatedPoint, hash_eigen<Eigen::Vector3i>> voxelindex_to_accpoint;

    Eigen::Vector3d ref_coord;
    Eigen::Vector3d ref_point;
    Eigen::Vector3d ref_normal;
    Eigen::Vector3d ref_color;
    double ref_alpha;
    Eigen::Vector3i voxel_index;
    for (int i = 0; i < V.rows(); i++) {
        ref_point = Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2));
        ref_coord = (ref_point - voxel_min_bound).array() / voxel_size3.array();
        voxel_index << int(floor(ref_coord(0))),
                       int(floor(ref_coord(1))),
                       int(floor(ref_coord(2)));
        if (has_normals) {
            ref_normal = Eigen::Vector3d(N(i, 0), N(i, 1), N(i, 2));
        }
        if (has_colors) {
            ref_color = Eigen::Vector3d(C(i, 0), C(i, 1), C(i, 2));
            ref_alpha = color_has_alpha_channel ? C(i, 3) : 1.0;
        }
        voxelindex_to_accpoint[voxel_index].AddPoint(ref_point, ref_normal, ref_color, ref_alpha);
    }

    int num_output_points = voxelindex_to_accpoint.size();
    outV.resize(num_output_points, 3);
    if (has_normals) {
        outN.resize(num_output_points, 3);
    }
    if (has_colors) {
        int dim = color_has_alpha_channel ? 4 : 3;
        outC.resize(num_output_points, dim);
    }

    int count = 0;
    for (auto accpoint : voxelindex_to_accpoint) {
        if (accpoint.second.GetNumPoints() < min_pts_per_bin) {
            continue;
        }
        Eigen::RowVector3d accpoint_pos = accpoint.second.GetAveragePoint();
        for (int i = 0; i < 3; i++) {outV(count, i) = accpoint_pos[i];}

        if (has_normals) {
            Eigen::RowVector3d accpoint_nml = accpoint.second.GetAverageNormal();
            for (int i = 0; i < 3; i++) {outN(count, i) = accpoint_nml[i];}
        }
        if (has_colors) {
            Eigen::Vector3d accpoint_rgb = accpoint.second.GetAverageColor();
            for (int i = 0; i < 3; i++) {outC(count, i) = accpoint_rgb[i];}
            if (color_has_alpha_channel) {
                outC(count, 3) = accpoint.second.GetAverageAlpha();
            }
        }
        count += 1;
    }
    if (count < num_output_points) {
        outV.conservativeResize(count, 3);
        if (has_normals) {
            outN.conservativeResize(count, 3);
        }
        if (has_colors) {
            int dim = color_has_alpha_channel ? 4 : 3;
            outC.conservativeResize(count, dim);
        }
    }
}

} // namespace




const char* sample_mesh_poisson_disk_doc = R"Qu8mg5v7(
Downsample a point set (possibly on a mesh) so that samples are approximately evenly spaced.
This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

Parameters
----------
v : #v by 3 array of mesh vertex positions
f : #f by 3 array of mesh face indices
num_samples: desired number of Poisson Disk samples. Note that the actual number of returned samples
             will not be exactly this value (see sample_num_tolerance) to control the range of possible
             returned samples.
             Note: If this value <= 0, then the parameter radius is used to decide the number of samples
radius : desired separation between points, if num_samples <= 0, then this value is used to determine the
         sampling (-1.0, by default).
use_geodesic_distance : Use geodesic distance on the mesh downsampling. (True by default).
best_choice_sampling : When downsampling, always keep the sample that will remove the
                       fewest number of samples. (True by default).
random_seed : A random seed used to generate the samples.
              Passing in 0 will use the current time. (0 by default).
sample_num_tolerance: If you requested a target number of samples, by passsing num_samples > 0, then this function will return
                      between (1 - sample_num_tolerance) * num_samples and (1 + sample_num_tolerance) * num_samples.
                      Setting a very small value for this parameter will increase convergence time. (0.04 by default).
oversampling_factor: To generate Poisson disk samples, we first generate a very dense (uniform) random sampling of the mesh, then
                     prune these down to have the Poisson disk property. This parameter controls how many dense samples are generated.
                     i.e. we generate oversampling_factor * num_samples samples (if you passed in radius, we estimate num_samples from
                     the input points and radius). This parameter must be >= 1.0. (Default 40.0).
Returns
-------
A (m,) shaped array of face indices into f where m is the number of Poisson-disk samples
A (m, 3) shaped array of barycentric coordinates where m is the number of Poisson-disk samples

)Qu8mg5v7";
npe_function(sample_mesh_poisson_disk)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_arg(num_samples, int)
npe_default_arg(radius, double, 0.0)
npe_default_arg(use_geodesic_distance, bool, true)
npe_default_arg(best_choice_sampling, bool, true)
npe_default_arg(random_seed, unsigned int, 0)
npe_default_arg(sample_num_tolerance, float, 0.04)
npe_default_arg(oversampling_factor, float, 40.0)
npe_doc(sample_mesh_poisson_disk_doc)
npe_begin_code()
{
    validate_mesh(v, f);

    if (num_samples <= 0 && radius <= 0.0) {
        throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
    }
    if (sample_num_tolerance > 1.0 || sample_num_tolerance <= 0.0) {
        throw pybind11::value_error("sample_num_tolerance must be in (0, 1]");
    }

    typedef MyMesh MeshType;
    typedef EigenDenseLike<npe_Matrix_v> EigenRetBC;
    typedef EigenVertexIndexSampler<MeshType> PoissonDiskSampler;
    typedef EigenBarycentricSampler<MyMesh, EigenRetBC> MonteCarloSampler;

    MyMesh input_mesh;
    vcg_mesh_from_vf(v, f, input_mesh);
    typename tri::SurfaceSampling<MeshType, PoissonDiskSampler>::PoissonDiskParam pp;
    //    int t0 = clock();

    if(radius > 0 && num_samples <= 0) {
        num_samples = tri::SurfaceSampling<MeshType, PoissonDiskSampler>::ComputePoissonSampleNum(input_mesh, radius);
    }

    pp.pds.sampleNum = num_samples;
    pp.randomSeed = random_seed;
    pp.geodesicDistanceFlag = use_geodesic_distance;
    pp.bestSampleChoiceFlag = best_choice_sampling;

    // Dense barycentric coordinates and face indices
    const int num_dense_samples = std::max(10000, int(num_samples * oversampling_factor));
    EigenRetBC dense_bc(num_dense_samples, 3);
    typename MonteCarloSampler::IndexArray dense_fi(num_dense_samples);
    MeshType montecarlo_mesh;
    MonteCarloSampler mcSampler(&input_mesh, &montecarlo_mesh, dense_bc, dense_fi);

    // We overallocate a bit because we could end up with more samples
    typename PoissonDiskSampler::IndexArray dense_vi(int(num_samples * (1.0 + sample_num_tolerance)));
    PoissonDiskSampler pdSampler(montecarlo_mesh, dense_vi);

    if(random_seed > 0) {
        tri::SurfaceSampling<MeshType, MonteCarloSampler>::SamplingRandomGenerator().initialize(random_seed);
        tri::SurfaceSampling<MeshType, PoissonDiskSampler>::SamplingRandomGenerator().initialize(random_seed);
    }

    // Generate dense samples on the mesh
    tri::SurfaceSampling<MeshType, MonteCarloSampler>::Montecarlo(input_mesh, mcSampler, num_dense_samples);
    tri::UpdateBounding<MeshType>::Box(montecarlo_mesh);
    //    int t1=clock();
    //    pp.pds.montecarloTime = t1-t0;

    // TODO: Adaptive radius would be nice actually!
    const double radiusVariance = 1;
    if(radiusVariance !=1)
    {
      pp.adaptiveRadiusFlag = true;
      pp.radiusVariance = radiusVariance;
    }

    if (radius <= 0.0 && num_samples > 0) {
        tri::SurfaceSampling<MeshType, PoissonDiskSampler>::PoissonDiskPruningByNumber(pdSampler, montecarlo_mesh, num_samples, radius, pp, sample_num_tolerance);

    } else {
        tri::SurfaceSampling<MeshType, PoissonDiskSampler>::PoissonDiskPruning(pdSampler, montecarlo_mesh, radius, pp);
    }
    //    int t2=clock();
    //    pp.pds.totalTime = t2-t0;

    EigenRetBC ret_bc(pdSampler.vcount, 3);
    typename MonteCarloSampler::IndexArray ret_fi(pdSampler.vcount);
    for (int i = 0; i < pdSampler.vcount; i++) {
        const std::ptrdiff_t dense_idx = dense_vi[i];
        ret_bc.row(i) = dense_bc.row(dense_idx);
        ret_fi[i] = dense_fi[dense_idx];
    }
    return std::make_tuple(npe::move(ret_fi), npe::move(ret_bc));
}
npe_end_code()


const char* sample_mesh_random_doc = R"Qu8mg5v7(
Generate uniformly distributed random point samples on a mesh

Parameters
----------
v : #v by 3 array of mesh vertex positions
f : #f by 3 array of mesh face indices
num_samples : The number of samples to generate
random_seed : A random seed used to generate the samples.
              Passing in 0 will use the current time. (0 by default).
Returns
-------
A (num_samples,) shaped array of face indices into f where
A (num_samples, 3) shaped array of barycentric coordinates

)Qu8mg5v7";
npe_function(sample_mesh_random)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_arg(num_samples, int)
npe_default_arg(random_seed, unsigned int, 0)
npe_doc(sample_mesh_random_doc)
npe_begin_code()
{
    validate_mesh(v, f);

    typedef EigenDenseLike<npe_Matrix_v> EigenRetBC;
    typedef EigenBarycentricSampler<MyMesh, EigenRetBC> MonteCarloSampler;

    MyMesh m;
    vcg_mesh_from_vf(v, f, m);

    EigenRetBC ret_bc(num_samples, 3);
    typename MonteCarloSampler::IndexArray ret_fi(num_samples);

    MonteCarloSampler mrs(&m, nullptr, ret_bc, ret_fi);

    if(random_seed > 0) {
        tri::SurfaceSampling<MyMesh, MonteCarloSampler>::SamplingRandomGenerator().initialize(random_seed);
    }

    tri::SurfaceSampling<MyMesh, MonteCarloSampler>::Montecarlo(m, mrs, num_samples);

    mrs.trim();
    return std::make_tuple(npe::move(ret_fi), npe::move(ret_bc));
}
npe_end_code()


const char* downsample_point_cloud_poisson_disk_doc = R"Qu8mg5v7(
Downsample a point set so that samples are approximately evenly spaced.
This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

Parameters
----------
v : #v by 3 array of vertex positions
n : #v by 3 array of vertex normals
num_samples: desired number of Poisson Disk samples. Note that the actual number of returned samples
             will not be exactly this value (see sample_num_tolerance) to control the range of possible
             returned samples.
             Note: If this value <= 0, then the parameter radius is used to decide the number of samples
radius : desired separation between points, if num_samples <= 0, then this value is used to determine the
         sampling (-1.0, by default).
best_choice_sampling : When downsampling, always keep the sample that will remove the
                       fewest number of samples. (True by default).
random_seed : A random seed used to generate the samples.
              Passing in 0 will use the current time. (0 by default).
sample_num_tolerance: If you requested a target number of samples, by passsing num_samples > 0, then this function will return
                      between (1 - sample_num_tolerance) * num_samples and (1 + sample_num_tolerance) * num_samples.
                      Setting a very small value for this parameter will increase convergence time. (0.04 by default).

Returns
-------
A (m,) shaped array of indices into v where m is the number of Poisson-disk samples

)Qu8mg5v7";
npe_function(downsample_point_cloud_poisson_disk)
npe_arg(v, dense_float, dense_double)
npe_arg(num_samples, int)
npe_default_arg(radius, double, 0.0)
npe_default_arg(best_choice_sampling, bool, true)
npe_default_arg(random_seed, unsigned int, 0)
npe_default_arg(sample_num_tolerance, float, 0.04)
npe_doc(downsample_point_cloud_poisson_disk_doc)
npe_begin_code()
{
  MyMesh m;
  vcg_mesh_from_v(v, m);

  if (num_samples <= 0 && radius <= 0.0) {
      throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
  }
  if (sample_num_tolerance > 1.0 || sample_num_tolerance <= 0.0) {
      throw pybind11::value_error("sample_num_tolerance must be in (0, 1]");
  }

  typedef EigenVertexIndexSampler<MyMesh> PoissonDiskSampler;
  typedef PoissonDiskSampler::IndexArray EigenRetI;

  EigenRetI ret_i;
  PoissonDiskSampler mps(m, ret_i);

  typename tri::SurfaceSampling<MyMesh, PoissonDiskSampler>::PoissonDiskParam pp;
  typename tri::SurfaceSampling<MyMesh, PoissonDiskSampler>::PoissonDiskParam::Stat pds;
  pp.pds = pds;
  pp.bestSampleChoiceFlag = best_choice_sampling;
  pp.geodesicDistanceFlag = false;
  pp.randomSeed = random_seed;

  if(random_seed) {
      tri::SurfaceSampling<MyMesh, PoissonDiskSampler>::SamplingRandomGenerator().initialize(random_seed);
  }

  if (radius <= 0.0 && num_samples > 0) {
      num_samples = std::min(num_samples, (int)v.rows());
      tri::SurfaceSampling<MyMesh, PoissonDiskSampler>::PoissonDiskPruningByNumber(mps, m, num_samples, radius, pp, sample_num_tolerance);
  } else if (radius > 0.0 && num_samples <= 0) {
      tri::SurfaceSampling<MyMesh, PoissonDiskSampler>::PoissonDiskPruning(mps, m, radius, pp);
  }
  mps.trim();
  return npe::move(ret_i);
}
npe_end_code()


npe_function(__internal_downsample_point_cloud_voxel_grid)
npe_arg(v, dense_float, dense_double)
npe_arg(n, dense_float, dense_double)
npe_arg(c, dense_float, dense_double)
npe_arg(voxel_size_x, double)
npe_arg(voxel_size_y, double)
npe_arg(voxel_size_z, double)
npe_arg(voxel_min_x, double)
npe_arg(voxel_min_y, double)
npe_arg(voxel_min_z, double)
npe_arg(voxel_max_x, double)
npe_arg(voxel_max_y, double)
npe_arg(voxel_max_z, double)
npe_arg(min_points_per_voxel, int)
npe_begin_code()
{
    EigenDense<npe_Scalar_v> outV;
    EigenDense<npe_Scalar_n> outN;
    EigenDense<npe_Scalar_c> outC;
    Eigen::Vector3d voxel_size(voxel_size_x, voxel_size_y, voxel_size_z);
    Eigen::Vector3d voxel_min_bound(voxel_min_x, voxel_min_y, voxel_min_z);
    Eigen::Vector3d voxel_max_bound(voxel_max_x, voxel_max_y, voxel_max_z);

    downsample_point_cloud_to_voxels(v, n, c,
                                     voxel_size, voxel_min_bound, voxel_max_bound,
                                     outV, outN, outC, min_points_per_voxel);

    return std::make_tuple(npe::move(outV), npe::move(outN), npe::move(outC));
}


npe_end_code()
const char* estimate_normals_doc = R"Qu8mg5v7(
Estimate normals for a point cloud by locally fitting a plane to a small neighborhood of points

Parameters
----------
v : #v by 3 array of vertex positions (each row is a vertex)
k : Number of nearest neighbors to use in the estimate for the normal of a point. Default: 10.
smoothing_iterations : Number of smoothing iterations to apply to the estimated normals. Default: 0.

Returns
-------
A #v x 3 array of normals where each row i is the normal at vertex v[i]

)Qu8mg5v7";
npe_function(estimate_point_cloud_normals)
npe_arg(v, dense_float, dense_double)
npe_default_arg(k, int, 10)
npe_default_arg(smoothing_iterations, int, 0)
npe_doc(estimate_normals_doc)
npe_begin_code()
{
    MyMesh m;
    vcg_mesh_from_v(v, m);

    tri::PointCloudNormal<MyMesh>::Param p;
    p.fittingAdjNum = k;
    p.smoothingIterNum = smoothing_iterations;
    p.viewPoint = vcg::Point3d(0.0f, 0.0f, 0.0f);
    p.useViewPoint = false;

    tri::PointCloudNormal<MyMesh>::Compute(m, p, (vcg::CallBackPos*) nullptr);

    npe_Matrix_v ret(m.vn, 3);

    int vcount = 0;
    for (MyMesh::VertexIterator vit = m.vert.begin(); vit != m.vert.end(); vit++) {
        ret(vcount, 0) = vit->N()[0];
        ret(vcount, 1) = vit->N()[1];
        ret(vcount, 2) = vit->N()[2];
        vcount += 1;
    }

    return npe::move(ret);

}
npe_end_code()

