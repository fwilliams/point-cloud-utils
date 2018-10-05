
// stuff to define the mesh
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/component_ep.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/normal.h>

// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

#include <cstdlib>

#include <sys/timeb.h>
#include <iostream>
#include <string>


class BaseVertex;
class BaseEdge;
class BaseFace;

struct BaseUsedTypes: public vcg::UsedTypes<vcg::Use<BaseVertex>::AsVertexType,vcg::Use<BaseEdge>::AsEdgeType,vcg::Use<BaseFace>::AsFaceType>{};

class BaseVertex  : public vcg::Vertex< BaseUsedTypes,
	vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  > {};

class BaseEdge : public vcg::Edge< BaseUsedTypes> {};

class BaseFace    : public vcg::Face< BaseUsedTypes,
	vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Mark, vcg::face::EmptyEdgePlane > {};

class BaseMesh    : public vcg::tri::TriMesh<std::vector<BaseVertex>, std::vector<BaseFace> > {};


class RTVertex;
class RTEdge;
class RTFace;

struct RTUsedTypes: public vcg::UsedTypes<vcg::Use<RTVertex>::AsVertexType,vcg::Use<RTEdge>::AsEdgeType,vcg::Use<RTFace>::AsFaceType>{};

class RTVertex  : public vcg::Vertex< RTUsedTypes,
	vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  > {};

class RTEdge : public vcg::Edge< RTUsedTypes> {};

class RTFace    : public vcg::Face< RTUsedTypes,
	vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::EdgePlane, vcg::face::Mark, vcg::face::BitFlags > {};

class RTMesh    : public vcg::tri::TriMesh<std::vector<RTVertex>, std::vector<RTFace> > {};


using namespace vcg;

void Usage()
{
	printf( "\nUsage:  trimesh_closest mesh.ply samplenum sampledistance(as fraction of bboxdiag)");
	exit(-1);
}

// Testing of closest point on a mesh functionalities
// Two main options
// - using or not precomputed edges and planes
// - using the simple wrapper or the basic functions of the grid.
// - using the fn as size of the grid or the edge lenght as cell side

template <class MeshType, bool useEdge,bool useWrap, bool useFaceNumForGrid>
bool UnitTest_Closest(const char *filename1, int sampleNum, float dispPerc, std::vector<int> resultVec)
{
  MeshType mr;
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::FaceType FaceType;
  typedef GridStaticPtr<FaceType, ScalarType> TriMeshGrid;

  int startOpen=clock();
  int err=vcg::tri::io::Importer<MeshType>::Open(mr,filename1);
  tri::UpdateBounding<MeshType>::Box(mr);
//  tri::UpdateNormals<MeshType>::PerFaceNormalized(mr);
  tri::UpdateNormal<MeshType>::PerFace(mr);
  float dispAbs = mr.bbox.Diag()*dispPerc;
  if(err)
  {
      std::cerr << "Unable to open mesh " << filename1 << " : " << vcg::tri::io::Importer<MeshType>::ErrorMsg(err) << std::endl;
      exit(-1);
  }
  int endOpen = clock();
  printf("Loading %6.3f - ",float(endOpen-startOpen)/CLOCKS_PER_SEC);

  int startSampling = clock();

  std::vector<Point3f> MontecarloSamples;
  // First step build the sampling
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  BaseSampler mcSampler(MontecarloSamples);
  tri::SurfaceSampling<MeshType,BaseSampler>::SamplingRandomGenerator().initialize(123);
  tri::SurfaceSampling<MeshType,BaseSampler>::Montecarlo(mr, mcSampler, sampleNum);
  math::MarsenneTwisterRNG rnd;
  rnd.initialize(123);
  for(size_t i=0;i<MontecarloSamples.size();++i)
  {
    Point3f pp(rnd.generate01(),rnd.generate01(),rnd.generate01());
    pp = (pp+Point3f(-0.5f,-0.5f,-0.5f))*2.0f;
    pp*=rnd.generate01()*dispAbs;
    MontecarloSamples[i]+=pp;
  }
  int endSampling = clock();

  printf("Sampling  %6.3f - ",float(endSampling-startSampling)/CLOCKS_PER_SEC);

  int startGridInit = clock();
  TriMeshGrid TRGrid;
  if(useFaceNumForGrid)
  {
    TRGrid.Set(mr.face.begin(),mr.face.end(),mr.FN()*2);
  }
  else
  {
    float avgEdge = tri::Stat<MeshType>::ComputeEdgeLengthAverage(mr);
    TRGrid.SetWithRadius(mr.face.begin(),mr.face.end(),avgEdge*2);
  }

  if(useEdge)
     tri::UpdateComponentEP<MeshType>::Set(mr);

  int endGridInit = clock();
  printf("Grid Init %6.3f - ",float(endGridInit-startGridInit)/CLOCKS_PER_SEC);

  const ScalarType maxDist=std::max(dispAbs*10.0f,mr.bbox.Diag()/1000.f);
  CoordType closest;
  ScalarType dist;
  int startGridQuery = clock();
  double avgDist=0;
  resultVec.resize(MontecarloSamples.size());
  if(useEdge && useWrap)
    for(size_t i=0;i<MontecarloSamples.size();++i)
    {
      resultVec[i]=tri::Index(mr,tri::GetClosestFaceEP(mr,TRGrid,MontecarloSamples[i], maxDist,dist,closest));
      if(resultVec[i]) avgDist += double(dist);
    }
  if(!useEdge && useWrap)
    for(size_t i=0;i<MontecarloSamples.size();++i)
    {
      resultVec[i]=tri::Index(mr,tri::GetClosestFaceBase(mr,TRGrid,MontecarloSamples[i], maxDist,dist,closest));
      if(resultVec[i]) avgDist += double(dist);
    }
  if(useEdge && !useWrap)
  {
    typedef tri::FaceTmark<MeshType> MarkerFace;
    MarkerFace mf;
    mf.SetMesh(&mr);
    face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
    for(size_t i=0;i<MontecarloSamples.size();++i)
    {
      resultVec[i]=tri::Index(mr,TRGrid.GetClosest(PDistFunct,mf,MontecarloSamples[i],maxDist,dist,closest));
      if(resultVec[i]) avgDist += double(dist);
    }
  }
  if(!useEdge && !useWrap)
  {
    typedef tri::FaceTmark<MeshType> MarkerFace;
    MarkerFace mf;
    mf.SetMesh(&mr);
    face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
    for(size_t i=0;i<MontecarloSamples.size();++i)
    {
      resultVec[i]=tri::Index(mr,TRGrid.GetClosest(PDistFunct,mf,MontecarloSamples[i],maxDist,dist,closest));
      if(resultVec[i]) avgDist += double(dist);
    }
  }

  int endGridQuery = clock();
  printf("Grid Size %3i %3i %3i - ",TRGrid.siz[0],TRGrid.siz[1],TRGrid.siz[2]);
  printf("Avg dist %6.9lf - ",avgDist / float(MontecarloSamples.size()));
  printf("Grid Query %6.3f \n", float(endGridQuery-startGridQuery)/CLOCKS_PER_SEC);
  return true;
}

int main(int argc ,char**argv)
{
  if(argc<3) Usage();
  float dispPerc = atof(argv[3]);
  int sampleNum = atoi(argv[2]);
  std::vector<int> resultVecRT11;
  std::vector<int> resultVecRT01;
  std::vector<int> resultVecRT00;
  std::vector<int> resultVecRT10;
  std::vector<int> resultVecBS01;
  std::vector<int> resultVecBS00;
  UnitTest_Closest<RTMesh, true,  true,  true>     (argv[1],sampleNum,dispPerc,resultVecRT11);
  UnitTest_Closest<RTMesh, true,  true,  false>    (argv[1],sampleNum,dispPerc,resultVecRT11);
  UnitTest_Closest<RTMesh, true,  false, true>    (argv[1],sampleNum,dispPerc,resultVecRT01);
  UnitTest_Closest<RTMesh, true,  false, false>    (argv[1],sampleNum,dispPerc,resultVecRT00);
  UnitTest_Closest<RTMesh, false, true,  true>   (argv[1],sampleNum,dispPerc,resultVecRT10);
  UnitTest_Closest<RTMesh, false, true,  false>   (argv[1],sampleNum,dispPerc,resultVecRT10);
  UnitTest_Closest<RTMesh, false, false, true>   (argv[1],sampleNum,dispPerc,resultVecRT10);
  UnitTest_Closest<RTMesh, false, false, false>   (argv[1],sampleNum,dispPerc,resultVecRT10);

  UnitTest_Closest<BaseMesh,false, true,  true>  (argv[1],sampleNum,dispPerc,resultVecBS01);
  UnitTest_Closest<BaseMesh,false, true,  false> (argv[1],sampleNum,dispPerc,resultVecBS01);
  UnitTest_Closest<BaseMesh,false, false, true>(argv[1],sampleNum,dispPerc,resultVecBS01);
  UnitTest_Closest<BaseMesh,false, false, false>(argv[1],sampleNum,dispPerc,resultVecBS01);

  for(size_t i=0;i<resultVecRT11.size();++i)
  {
    if(resultVecRT11[i]!=resultVecRT01[i]) printf("%lu is diff",i);
    if(resultVecRT11[i]!=resultVecRT00[i]) printf("%lu is diff",i);
    if(resultVecRT11[i]!=resultVecRT10[i]) printf("%lu is diff",i);
    if(resultVecRT11[i]!=resultVecBS00[i]) printf("%lu is diff",i);
    if(resultVecRT11[i]!=resultVecBS01[i]) printf("%lu is diff",i);
  }
  return 0;
}
