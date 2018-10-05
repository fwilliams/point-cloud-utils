#include <vector>

using namespace std;

// VCG headers for triangular mesh processing
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/closest.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
                                        Use<MyEdge>			::AsEdgeType,
																				Use<MyFace>			::AsFaceType>{};


class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Mark,vertex::Color4b, vertex::Qualityf>{};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyFace    : public Face  <MyUsedTypes, face::VertexRef,face::BitFlags,face::Mark, face::Normal3f> {};

class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};

// Uncomment only one of the two following lines to test different data structures
typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
//typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv)
{
  if (argc<2)
	{
		printf("\n");
    printf("    Compute an approximation of the shape diameter function\n");
    printf("    Usage: trimesh_intersection <filename> [angle samplenum]\n\n");
    printf("       <filename>        Mesh model for which to compute the sdf (PLY format).\n");
    printf("       angle             the wideness (degree) of the cone of ray that must be shot from each vertex (default 45)\n");
    printf("       samplenum         the oversampling factor (0 -> one ray, 1, 9 ray, 2-> 25 rays (default 2)\n");

		return 0;
	}

	MyMesh m;
 int t0=clock();
	// open a mesh
	int err = tri::io::Importer<MyMesh>::Open(m,argv[1]);
  if(err) {
		printf("Error in reading %s: '%s'\n",argv[1],tri::io::Importer<MyMesh>::ErrorMsg(err));
		exit(-1);  
	}
 // the other parameters
  float widenessRad = math::ToRad(20.0);

  if(argc>2) {
    widenessRad = math::ToRad(atof(argv[2]));
    printf("Setting wideness to %f degree\n",atof(argv[2]));
  }
  int n_samples=2;
  if(argc>3) n_samples = atoi(argv[3]);
  int samplePerVert = (n_samples*2+ 1)*(n_samples*2+ 1);
  printf("Using oversampling to %i  (%i sample per vertex)\n",n_samples,samplePerVert);


  // some cleaning to get rid of bad stuff
	int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
	int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);

	if (dup > 0 || unref > 0)
		printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);

  // updating
  tri::UpdateBounding<MyMesh>::Box(m);
  tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
	// Create a static grid (for fast indexing) and fill it
	TriMeshGrid static_grid;
	static_grid.Set(m.face.begin(), m.face.end());

  typedef MyMesh::ScalarType ScalarType;
  int t1=clock();
  float t;
  MyMesh::FaceType *rf;
  MyMesh::VertexIterator vi;
  float maxDist=m.bbox.Diag();
  float offset= maxDist / 10000.0;
  int totRay=0;

  ScalarType deltaRad=widenessRad/(ScalarType)(n_samples*2);
  if(n_samples==0) deltaRad=0;

  tri::UpdateQuality<MyMesh>::VertexConstant(m,0);
  for(vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    vcg::Ray3f ray;
    ray.SetOrigin((*vi).cP()-((*vi).cN()*offset));
    Point3f dir0 = -(*vi).cN();
    int cnt=0;
    ScalarType theta_init,phi_init,ro;
    dir0.ToPolarRad(ro,theta_init,phi_init);
    for (int x=-n_samples;x<=n_samples;x++)
      for (int y=-n_samples;y<=n_samples;y++)
      {
        ScalarType theta=theta_init+x*deltaRad;
        ScalarType phi=phi_init+y*deltaRad;

        if (theta<0) theta=2.0*M_PI+theta;

        Point3f dir;
        dir.FromPolarRad(ro,theta,phi);
        dir.Normalize();
        ray.SetDirection(dir);

        rf = tri::DoRay<MyMesh,TriMeshGrid>(m,static_grid,ray,maxDist,t);
        if(rf)
        {
          (*vi).Q()+=t;
          cnt++;
        }
      }
    if(cnt>0){
      (*vi).Q()/=cnt;
      totRay+=cnt;
    }
  }
  int t2 = clock();
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);
  tri::io::ExporterPLY<MyMesh>::Save(m,"SDF.ply",tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);

  printf("Initializated in %i msec\n",t1-t0);
  printf("Completed in %i msec\n",t2-t1);
  printf("Shoot %i rays and found %i intersections\n",m.VN()*samplePerVert,totRay);

return 0;
}

