#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags>{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};


/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V, and a #F x 3 matrix of face indices into a VCG mesh
 */
template <typename DerivedV, typename DerivedF>
static void vcg_mesh_from_vf(
    const Eigen::MatrixBase<DerivedV>& V,
    const Eigen::MatrixBase<DerivedF>& F,
    MyMesh& m) {

  MyMesh::VertexIterator vit = Allocator<MyMesh>::AddVertices(m, V.rows());
  std::vector<MyMesh::VertexPointer> ivp(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    ivp[i] = &*vit;
    vit->P() = MyMesh::CoordType(V(i, 0), V(i, 1), V(i, 2));
    vit++;
  }

  MyMesh::FaceIterator fit = Allocator<MyMesh>::AddFaces(m, F.rows());
  for (int i = 0; i < F.rows(); i++) {
    fit->V(0) = ivp[F(i, 0)];
    fit->V(1) = ivp[F(i, 1)];
    fit->V(2) = ivp[F(i, 2)];
    fit++;
  }

  tri::UpdateBounding<MyMesh>::Box(m);
}

/*
 * Copy a mesh stored as a #V x 3 matrix of vertices, V into a VCG mesh
 */
template <typename DerivedV>
static void vcg_mesh_from_v(
    const Eigen::MatrixBase<DerivedV>& V,
    MyMesh& m) {

  MyMesh::VertexIterator vit = Allocator<MyMesh>::AddVertices(m, V.rows());
  std::vector<MyMesh::VertexPointer> ivp(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    ivp[i] = &*vit;
    vit->P() = MyMesh::CoordType(V(i, 0), V(i, 1), V(i, 2));
    vit++;
  }

  tri::UpdateBounding<MyMesh>::Box(m);
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


npe_function(poisson_disk_sample)
npe_arg(v, dense_f32, dense_f64)
npe_arg(f, dense_i32, dense_i64)
npe_arg(radius, double)
npe_default_arg(approximate_geodesic_distance, bool, false)
npe_default_arg(best_choice_sampling, bool, false)
npe_begin_code()
  MyMesh m;
  vcg_mesh_from_vf(v, f, m);

  MyMesh subM;
  tri::MeshSampler<MyMesh> mps(subM);

  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam pp;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam::Stat pds;
  pp.pds = pds;
  pp.bestSampleChoiceFlag = best_choice_sampling;
  pp.geodesicDistanceFlag = approximate_geodesic_distance;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskPruning(mps, m, radius, pp);

  npe_Matrix_v ret;
  vcg_mesh_to_v(subM, ret);

  return npe::move(ret);
npe_end_code()


npe_function(cluster_vertices)
npe_arg(v, dense_f32, dense_f64)
npe_arg(cell_size, double)
npe_default_arg(approximate_geodesic_distance, bool, false)
npe_default_arg(best_choice_sampling, bool, false)
npe_begin_code()
  MyMesh m;
  vcg_mesh_from_v(v, m);

  MyMesh cluM;

  tri::Clustering<MyMesh, vcg::tri::AverageColorCell<MyMesh> > ClusteringGrid;
  ClusteringGrid.Init(m.bbox, 100000, cell_size);
  ClusteringGrid.AddPointSet(m);
  ClusteringGrid.ExtractMesh(cluM);

  npe_Matrix_v ret;
  vcg_mesh_to_v(cluM, ret);

  return npe::move(ret);
npe_end_code()


npe_function(random_sample)
npe_arg(v, dense_f32, dense_f64)
npe_arg(f, dense_i32, dense_i64)
npe_arg(num_samples, int)
npe_begin_code()
  MyMesh m;
  vcg_mesh_from_vf(v, f, m);

  MyMesh rndM;
  tri::MeshSampler<MyMesh> mrs(rndM);

  tri::SurfaceSampling<MyMesh, tri::MeshSampler<MyMesh>>::VertexUniform(m, mrs, num_samples);

  npe_Matrix_v ret;
  vcg_mesh_to_v(rndM, ret);

  return npe::move(ret);
npe_end_code()
