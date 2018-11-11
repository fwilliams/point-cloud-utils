#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

#include <fstream>
#include <iostream>


using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags> {};
class MyFace    : public Face<MyUsedTypes, face::FFAdj,  face::Normal3f, face::VertexRef, face::BitFlags> {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>> {};


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF, typename DerivedN>
static void vcg_mesh_from_vfn(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    const Eigen::MatrixBase<DerivedN>& N,
    MyMesh& m) {

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
 * Copy a VCG mesh into a #V x 3 matrix of vertices, V
 */
template <typename DerivedV, typename DerivedN>
static void vcg_mesh_to_vn(
    const MyMesh& m,
    Eigen::PlainObjectBase<DerivedV>& V,
    Eigen::PlainObjectBase<DerivedN>& N) {

  V.resize(m.vn, 3);
  N.resize(m.vn, 3);

  int vcount = 0;
  for (auto vit = m.vert.begin(); vit != m.vert.end(); vit++) {
    V(vcount, 0) = vit->P()[0];
    V(vcount, 1) = vit->P()[1];
    V(vcount, 2) = vit->P()[2];

    N(vcount, 0) = vit->N()[0];
    N(vcount, 1) = vit->N()[1];
    N(vcount, 2) = vit->N()[2];

    vcount += 1;
  }
}


/*
 * Copy a VCG mesh into a #V x 3 matrix of vertices, V
 */
template <typename DerivedV>
static void vcg_mesh_to_v(
    const MyMesh& m,
    Eigen::PlainObjectBase<DerivedV>& V) {

  V.resize(m.vn, 3);

  int vcount = 0;
  for (auto vit = m.vert.begin(); vit != m.vert.end(); vit++) {
    V(vcount, 0) = vit->P()[0];
    V(vcount, 1) = vit->P()[1];
    V(vcount, 2) = vit->P()[2];
    vcount += 1;
  }
}



const char* poisson_disk_sample_doc = R"Qu8mg5v7(
Downsample a point set (possibly on a mesh) so that samples are approximately evenly spaced.
This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

Parameters
----------
v : #v by 3 list of mesh vertex positions
f : #f by 3 list of mesh face indices
n : #v by 3 list of mesh vertex normals
radius : desired separation between points
use_geodesic_distance : Use geodesic distance on the mesh downsampling, False by default
best_choice_sampling : When downsampling, always keep the sample that will remove the
                       fewest number of samples, False by default

Returns
-------
A #pv x 3 matrix of points which are approximately evenly spaced and are a subset of the input v

)Qu8mg5v7";

npe_function(poisson_disk_sample)
npe_arg(v, dense_f32, dense_f64)
npe_arg(f, dense_i32, dense_i64)
npe_arg(n, dense_f32, dense_f64)
npe_arg(radius, double)
npe_default_arg(use_geodesic_distance, bool, false)
npe_default_arg(best_choice_sampling, bool, false)
npe_doc(poisson_disk_sample_doc)
npe_begin_code()

  if (v.rows() == 0 || v.cols() != 3) {
    std::stringstream ss;
    ss << "Invalid matrix size. Expected n by 3 (n > 0) but got shape =(" << v.rows() << ", " << v.cols() << ")";
    throw pybind11::value_error(ss.str());
  }

  MyMesh m;
  vcg_mesh_from_vfn(v, f, n, m);

  MyMesh subM;
  tri::MeshSampler<MyMesh> mps(subM);

  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam pp;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam::Stat pds;
  pp.pds = pds;
  pp.bestSampleChoiceFlag = best_choice_sampling;
  pp.geodesicDistanceFlag = use_geodesic_distance;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskPruning(mps, m, radius, pp);

  npe_Matrix_v ret_v;
  npe_Matrix_n ret_n;
  vcg_mesh_to_vn(subM, ret_v, ret_n);

  m.Clear();
  subM.Clear();
  mps.reset();
  return std::make_tuple(npe::move(ret_v), npe::move(ret_n));

npe_end_code()




const char* cluster_vertices_doc = R"Qu8mg5v7(
Divide the bounding box of a point cloud into cells and cluster vertices which lie in the samee cell

Parameters
----------
v : #v by 3 list of mesh vertex positions
cell_size : Dimension along one axis of the cells

Returns
-------
A #pv x 3 matrix of clustered points

)Qu8mg5v7";

npe_function(cluster_vertices)
npe_arg(v, dense_f32, dense_f64)
npe_arg(cell_size, double)
npe_doc(cluster_vertices_doc)
npe_begin_code()

  if (v.rows() == 0 || v.cols() != 3) {
    std::stringstream ss;
    ss << "Invalid matrix size. Expected n by 3 (n > 0) but got shape =(" << v.rows() << ", " << v.cols() << ")";
    throw pybind11::value_error(ss.str());
  }

  MyMesh m;
  npe_Matrix_v n(0, 0);
  Eigen::MatrixXi f(0, 0);
  vcg_mesh_from_vfn(v, f, n, m);

  MyMesh cluM;

  tri::Clustering<MyMesh, vcg::tri::AverageColorCell<MyMesh> > ClusteringGrid;
  ClusteringGrid.Init(m.bbox, 100000, cell_size);
  ClusteringGrid.AddPointSet(m);
  ClusteringGrid.ExtractMesh(cluM);

  npe_Matrix_v ret_v;
  vcg_mesh_to_v(cluM, ret_v);

  m.Clear();
  cluM.Clear();

  return npe::move(ret_v);

npe_end_code()


const char* random_sample_doc = R"Qu8mg5v7(
Generate uniformly distributed random point samples on a mesh

Parameters
----------
v : #v by 3 list of mesh vertex positions
f : #f by 3 list of mesh face indices
n : #v by 3 list of mesh vertex normals
num_samples : The number of samples to generate

Returns
-------
A #pv x 3 matrix of samples

)Qu8mg5v7";

npe_function(random_sample)
npe_arg(v, dense_f32, dense_f64)
npe_arg(f, dense_i32, dense_i64)
npe_arg(n, dense_f32, dense_f64)
npe_arg(num_samples, int)
npe_doc(random_sample_doc)
npe_begin_code()

  if (v.rows() == 0 || v.cols() != 3) {
    std::stringstream ss;
    ss << "Invalid matrix size. Expected n by 3 (n > 0) but got shape =(" << v.rows() << ", " << v.cols() << ")";
    throw pybind11::value_error(ss.str());
  }

  MyMesh m;
  vcg_mesh_from_vfn(v, f, n, m);

  MyMesh rndM;
  tri::MeshSampler<MyMesh> mrs(rndM);

  tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::VertexUniform(m, mrs, num_samples);

  npe_Matrix_v ret_v;
  npe_Matrix_n ret_n;
  vcg_mesh_to_vn(rndM, ret_v, ret_n);


  m.Clear();
  rndM.Clear();

  return std::make_tuple(npe::move(ret_v), npe::move(ret_n));

npe_end_code()

