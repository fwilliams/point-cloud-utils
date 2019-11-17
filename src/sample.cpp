#include <npe.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

#include <fstream>
#include <iostream>

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

template <class MeshType, typename VType, typename FType, typename NType>
class EigenMeshSampler
{
public:
    typedef typename MeshType::VertexType  VertexType;
    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::CoordType   CoordType;

    VType& v;
    FType& f;
    NType& n;

    int vcount = 0;
    int ncount = 0;
    int fcount = 0;

    bool return_normals = false;

    EigenMeshSampler(VType& v, FType& f, NType& n, bool has_normals, bool per_face_normals=false) : v(v), n(n), f(f), return_normals(has_normals) {
        perFaceNormal = per_face_normals;
    }

    EigenMeshSampler(VType& v, FType& f) : v(v), n(n), f(f), return_normals(false) {
        perFaceNormal = false;
    }

    bool perFaceNormal;  // default false; if true the sample normal is the face normal, otherwise it is interpolated

    void trim() {
        v.conservativeResize(vcount, 3);
        n.conservativeResize(ncount, 3);
        f.conservativeResize(fcount, 3);
    }

    void reset() {
        vcount = 0;
        fcount = 0;
        ncount = 0;
    }

    void maybe_resize()  {
        if (v.rows() <= vcount) {
            const int n_rows = v.rows() == 0 ? 1024 : v.rows();
            v.conservativeResize(2*n_rows, 3);
        }

        if (return_normals) {
            if (n.rows() <= ncount) {
                const int n_rows = n.rows() == 0 ? 1024 : n.rows();
                n.conservativeResize(2*n_rows, 3);
            }
        }
    }

    void AddVert(const VertexType &p) {
        maybe_resize();

        v(vcount, 0) = p.cP()[0];
        v(vcount, 1) = p.cP()[1];
        v(vcount, 2) = p.cP()[2];
        vcount += 1;

        if (return_normals) {
            n(ncount, 0) = p.cN()[0];
            n(ncount, 1) = p.cN()[1];
            n(ncount, 2) = p.cN()[2];
            ncount += 1;
        }
    }

    void AddFace(const FaceType &f, CoordType p) {
        maybe_resize();

        auto vertex = f.cP(0)*p[0] + f.cP(1)*p[1] +f.cP(2)*p[2];
        v(vcount, 0) = vertex[0];
        v(vcount, 1) = vertex[1];
        v(vcount, 2) = vertex[2];
        vcount += 1;

        if (return_normals) {
            auto normal = perFaceNormal ? f.cN() : f.cV(0)->N()*p[0] + f.cV(1)->N()*p[1] + f.cV(2)->N()*p[2];
            n(ncount, 0) = normal[0];
            n(ncount, 1) = normal[0];
            n(ncount, 2) = normal[0];
            ncount += 1;
        }
    }
}; // end class BaseSampler



} // namespace








const char* sample_mesh_poisson_disk_doc = R"Qu8mg5v7(
Downsample a point set (possibly on a mesh) so that samples are approximately evenly spaced.
This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

Parameters
----------
v : #v by 3 list of mesh vertex positions
f : #f by 3 list of mesh face indices
n : #v by 3 list of mesh vertex normals
num_samples: desired number of Poisson Disk samples. If this value <= 0, then the parameter radius is used to decide the number of samples
radius : desired separation between points, if num_samples <= 0, then this value is used to determine the sampling (-1.0, by default)
use_geodesic_distance : Use geodesic distance on the mesh downsampling, False by default
best_choice_sampling : When downsampling, always keep the sample that will remove the
                       fewest number of samples, False by default
random_seed : A random seed used to generate the samples, 0 by default will use the current time

Returns
-------
A #pv x 3 matrix of points which are approximately evenly spaced
A #pv x 3 matrix of normals if normals were passed in (else an empty array)

)Qu8mg5v7";

npe_function(sample_mesh_poisson_disk)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_arg(n, npe_matches(v))
npe_arg(num_samples, int)
npe_default_arg(radius, double, 0.0)
npe_default_arg(use_geodesic_distance, bool, false)
npe_default_arg(best_choice_sampling, bool, false)
npe_default_arg(random_seed, unsigned int, 0)
npe_doc(sample_mesh_poisson_disk_doc)
npe_begin_code()

    validate_mesh(v, f, n);

    if (num_samples <= 0 && radius <= 0.0) {
        throw pybind11::value_error("Cannot have both num_samples <= 0 and radius <= 0");
    }

    double radiusVariance = 1;
    double PruningByNumberTolerance = 0.04;

    typedef MyMesh MeshType;
    typedef EigenDenseLike<npe_Matrix_v> EigenRetV;
    typedef EigenDenseLike<npe_Matrix_n> EigenRetN;
    typedef EigenDenseLike<npe_Matrix_f> EigenRetF;
    typedef EigenMeshSampler<MeshType, EigenRetV, EigenRetF, EigenRetN> BaseSampler;
    typedef tri::MeshSampler<MeshType> MontecarloSampler;

    MyMesh m;
    vcg_mesh_from_vfn(v, f, n, m);
    typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;
    int t0=clock();

    if(radius > 0 && num_samples <= 0) {
        num_samples = tri::SurfaceSampling<MeshType,BaseSampler>::ComputePoissonSampleNum(m, radius);
    }

    pp.pds.sampleNum = num_samples;
    pp.randomSeed = random_seed;
    pp.geodesicDistanceFlag = use_geodesic_distance;
    pp.bestSampleChoiceFlag = best_choice_sampling;

    MeshType MontecarloMesh;

    // First step build the sampling
    MontecarloSampler mcSampler(MontecarloMesh);
    EigenRetV ret_v(int(num_samples*1.2), 3);
    EigenRetN ret_n(int(num_samples*1.2), 3);
    EigenRetF ret_f(0, 3);
    const bool return_normals = (n.rows() != 0);
    BaseSampler pdSampler(ret_v, ret_f, ret_n, return_normals);

    if(random_seed) {
        tri::SurfaceSampling<MeshType,MontecarloSampler>::SamplingRandomGenerator().initialize(random_seed);
    }
    tri::SurfaceSampling<MeshType,MontecarloSampler>::Montecarlo(m, mcSampler, std::max(10000, num_samples*40));
    tri::UpdateBounding<MeshType>::Box(MontecarloMesh);
    //    int t1=clock();
    //    pp.pds.montecarloTime = t1-t0;

    if(radiusVariance !=1)
    {
      pp.adaptiveRadiusFlag=true;
      pp.radiusVariance=radiusVariance;
    }
    if(num_samples == 0) {
        tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskPruning(pdSampler, MontecarloMesh, radius, pp);
    } else {
        tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskPruningByNumber(pdSampler, MontecarloMesh, num_samples, radius,pp,PruningByNumberTolerance);
    }
    //    int t2=clock();
    //    pp.pds.totalTime = t2-t0;


    pdSampler.trim();
    return std::make_tuple(npe::move(ret_v), npe::move(ret_n));

npe_end_code()



const char* sample_mesh_random_doc = R"Qu8mg5v7(
Generate uniformly distributed random point samples on a mesh

Parameters
----------
v : #v by 3 list of mesh vertex positions
f : #f by 3 list of mesh face indices
n : #v by 3 list of mesh vertex normals (or 0 by 3 if no normals available)
num_samples : The number of samples to generate

Returns
-------
A #pv x 3 array of samples
A #pv x 3 array of normals if they were passed in, otherwise a (0, 3) empty array

)Qu8mg5v7";

npe_function(sample_mesh_random)
npe_arg(v, dense_float, dense_double)
npe_arg(f, dense_int, dense_longlong, dense_uint, dense_ulonglong)
npe_arg(n, npe_matches(v))
npe_arg(num_samples, int)
npe_doc(sample_mesh_random_doc)
npe_begin_code()

    validate_mesh(v, f, n);

    typedef EigenDenseLike<npe_Matrix_v> EigenRetV;
    typedef EigenDenseLike<npe_Matrix_n> EigenRetN;
    typedef EigenDenseLike<npe_Matrix_f> EigenRetF;
    typedef EigenMeshSampler<MyMesh, EigenRetV, EigenRetF, EigenRetN> MonteCarloSampler;

    MyMesh m;
    vcg_mesh_from_vfn(v, f, n, m);

    EigenRetV ret_v(num_samples, 3);
    EigenRetN ret_n(num_samples, 3);
    EigenRetF ret_f(0, 3);
    const bool return_normals = (n.rows() != 0);
    MonteCarloSampler mrs(ret_v, ret_f, ret_n, return_normals);

    tri::SurfaceSampling<MyMesh, MonteCarloSampler>::Montecarlo(m, mrs, num_samples);

//    m.Clear();

    mrs.trim();
    return std::make_tuple(npe::move(ret_v), npe::move(ret_n));

npe_end_code()


//const char* prune_point_cloud_poisson_disk_doc = R"Qu8mg5v7(
//Downsample a point set so that samples are approximately evenly spaced.
//This function uses the method in "Parallel Poisson Disk Sampling with Spectrum Analysis on Surface"
//(http://graphics.cs.umass.edu/pubs/sa_2010.pdf)

//Parameters
//----------
//v : #v by 3 list of mesh vertex positions
//n : #v by 3 list of mesh vertex normals
//radius : desired separation between points
//best_choice_sampling : When downsampling, always keep the sample that will remove the
//                       fewest number of samples, False by default

//Returns
//-------
//A #pv x 3 matrix of points which are approximately evenly spaced and are a subset of the input v

//)Qu8mg5v7";

//npe_function(sample_point_cloud_poisson_disk)
//npe_arg(v, dense_float, dense_double)
//npe_arg(n, npe_matches(v))
//npe_arg(radius, double)
//npe_default_arg(best_choice_sampling, bool, false)
//npe_doc(sample_mesh_poisson_disk_doc)
//npe_begin_code()

//  Eigen::MatrixXi f(0, 3);

//  MyMesh m;
//  vcg_mesh_from_vfn(v, f, n, m);

//  MyMesh subM;
//  tri::MeshSampler<MyMesh> mps(subM);

//  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam pp;
//  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam::Stat pds;
//  pp.pds = pds;
//  pp.bestSampleChoiceFlag = best_choice_sampling;
//  pp.geodesicDistanceFlag = false;
//  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskPruning(mps, m, radius, pp);

//  EigenDenseLike<npe_Matrix_v> ret_v;
//  EigenDenseLike<npe_Matrix_n> ret_n;
//  vcg_mesh_to_vn(subM, ret_v, ret_n, (n.rows() == 0) /*skip_normals*/);

//  m.Clear();
//  subM.Clear();
//  mps.reset();
//  return std::make_tuple(npe::move(ret_v), npe::move(ret_n));

//npe_end_code()





