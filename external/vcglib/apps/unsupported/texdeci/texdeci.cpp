#include <vector>
#include <limits>

#include <stdio.h>
#include <stdlib.h>

// stuff to define the mesh
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/complex/complex.h>

#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/clean.h>

// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// update
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/smooth.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric_tex.h>

using namespace vcg;
using namespace tri;

// The class prototypes.
class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes,
  vertex::VFAdj,
  vertex::Coord3f,
  vertex::Normal3f,
  vertex::Mark,
  vertex::BitFlags  >{
  };

class MyEdge : public Edge< MyUsedTypes> {};

typedef BasicVertexPair<MyVertex> VertexPair;

class MyFace    : public Face< MyUsedTypes,
  face::VFAdj,
  face::VertexRef,
  face::BitFlags,
  face::WedgeTexCoord2f> {};

// the main mesh class
class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};

class MyTriEdgeCollapseQTex: public TriEdgeCollapseQuadricTex< MyMesh, VertexPair, MyTriEdgeCollapseQTex, QuadricTexHelper<MyMesh> > {
            public:
            typedef  TriEdgeCollapseQuadricTex< MyMesh,  VertexPair, MyTriEdgeCollapseQTex, QuadricTexHelper<MyMesh> > TECQ;
            inline MyTriEdgeCollapseQTex(  const VertexPair &p, int i,BaseParameterClass *pp) :TECQ(p,i,pp){}
};


void TexDecimation(MyMesh &m, bool CleaningFlag,int TargetFaceNum)
{
  tri::TriEdgeCollapseQuadricTexParameter pp;

  pp.SetDefaultParams();
  if(CleaningFlag){
        int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
        int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
        printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",dup,unref);
    }

    printf("reducing it to %i\n",TargetFaceNum);
    int t1=clock();

    tri::UpdateBounding<MyMesh>::Box(m);
    math::Quadric<double> QZero;
    QZero.SetZero();
    QuadricTexHelper<MyMesh>::QuadricTemp TD3(m.vert,QZero);
    QuadricTexHelper<MyMesh>::TDp3()=&TD3;

    std::vector<std::pair<vcg::TexCoord2<float>,Quadric5<double> > > qv;

    QuadricTexHelper<MyMesh>::Quadric5Temp TD(m.vert,qv);
    QuadricTexHelper<MyMesh>::TDp()=&TD;


    vcg::LocalOptimization<MyMesh> DeciSession(m, &pp);
//    cb(1,"Initializing simplification");
    DeciSession.Init<MyTriEdgeCollapseQTex>();

    DeciSession.SetTargetSimplices(TargetFaceNum);
    DeciSession.SetTimeBudget(0.1f);
  //	int startFn=m.fn;

    int faceToDel=m.fn-TargetFaceNum;
    int t2=clock();

    while( DeciSession.DoOptimization() && m.fn>TargetFaceNum )
    {
      printf("Simplifing heap size %i ops %i\n",int(DeciSession.h.size()),DeciSession.nPerfmormedOps);
    };

    DeciSession.Finalize<MyTriEdgeCollapseQTex>();

    int t3=clock();
    printf("mesh  %d %d Error %g \n",m.vn,m.fn,DeciSession.currMetric);
    printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
}

// mesh to simplify

int main(int argc, char**argv){

  int meshNum=argc-1;

//std::vector<MyMesh> meshVec(meshNum);

MyMesh meshVec[10];
int tt0=clock();
char buf[255];
int i;

for(i=0;i<meshNum;++i)
{
  int err=vcg::tri::io::Importer<MyMesh>::Open(meshVec[i],argv[i+1]);
  if(err)
  {
    printf("Unable to open mesh %s : '%s'\n",argv[i+1], vcg::tri::io::Importer<MyMesh>::ErrorMsg(err));
    exit(-1);
  }
  printf("mesh loaded %d %d \n",meshVec[i].vn,meshVec[i].fn);

 int t1=clock();
 tri::Smooth<MyMesh>::VertexCoordLaplacian(meshVec[i],5*i);

 TexDecimation(meshVec[i],true,meshVec[i].fn/2);
 int t2=clock();
 printf("%i %5.3f sec\n",i,float(t2-t1)/CLOCKS_PER_SEC);
 sprintf(buf,"out%i.ply",i);
 tri::io::ExporterPLY<MyMesh>::Save(meshVec[i],buf,false);
}

int tt1=clock();
printf("---Total %5.3f sec\n",float(tt1-tt0)/CLOCKS_PER_SEC);

for(int i=0;i<meshNum;++i)
{
 char buf[255];
 sprintf(buf,"out%i.ply",i);
 tri::io::ExporterPLY<MyMesh>::Save(meshVec[i],buf,tri::io::Mask::IOM_WEDGTEXCOORD,false);
}

//  TriEdgeCollapseQuadricParameter qparams;
//  qparams.QualityThr  =.3;
//  float TargetError=std::numeric_limits<float>::max();
//  bool CleaningFlag =false;
//     // parse command line.
//    for(int i=4; i < argc;)
//    {
//      if(argv[i][0]=='-')
//        switch(argv[i][1])
//      {
//        case 'H' : qparams.SafeHeapUpdate=true; printf("Using Safe heap option\n"); break;
//        case 'Q' : if(argv[i][2]=='y') { qparams.QualityCheck	= true;  printf("Using Quality Checking\n");	}
//                                  else { qparams.QualityCheck	= false; printf("NOT Using Quality Checking\n");	}                break;
//        case 'N' : if(argv[i][2]=='y') { qparams.NormalCheck	= true;  printf("Using Normal Deviation Checking\n");	}
//                                  else { qparams.NormalCheck	= false; printf("NOT Using Normal Deviation Checking\n");	}        break;
//        case 'O' : if(argv[i][2]=='y') { qparams.OptimalPlacement	= true;  printf("Using OptimalPlacement\n");	}
//                                  else { qparams.OptimalPlacement	= false; printf("NOT Using OptimalPlacement\n");	}        break;
//        case 'S' : if(argv[i][2]=='y') { qparams.ScaleIndependent	= true;  printf("Using ScaleIndependent\n");	}
//                                  else { qparams.ScaleIndependent	= false; printf("NOT Using ScaleIndependent\n");	}        break;
//        case 'B' : if(argv[i][2]=='y') { qparams.PreserveBoundary	= true;  printf("Preserving Boundary\n");	}
//                                  else { qparams.PreserveBoundary	= false; printf("NOT Preserving Boundary\n");	}        break;
//        case 'T' : if(argv[i][2]=='y') { qparams.PreserveTopology	= true;  printf("Preserving Topology\n");	}
//                                  else { qparams.PreserveTopology	= false; printf("NOT Preserving Topology\n");	}        break;
//        case 'q' :	qparams.QualityThr	= atof(argv[i]+2);	           printf("Setting Quality Thr to %f\n",atof(argv[i]+2)); 	 break;
//        case 'n' :	qparams.NormalThrRad = math::ToRad(atof(argv[i]+2));  printf("Setting Normal Thr to %f deg\n",atof(argv[i]+2)); break;
//        case 'b' :	qparams.BoundaryWeight  = atof(argv[i]+2);			printf("Setting Boundary Weight to %f\n",atof(argv[i]+2)); break;
//        case 'e' :	TargetError = float(atof(argv[i]+2));			printf("Setting TargetError to %g\n",atof(argv[i]+2)); break;
//        case 'P' :	CleaningFlag=true;  printf("Cleaning mesh before simplification\n"); break;

//        default  :  printf("Unknown option '%s'\n", argv[i]);
//          exit(0);
//      }
//      i++;
//    }



//  if(CleaningFlag){
//      int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
//      int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
//      printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",dup,unref);
//  }


//  printf("reducing it to %i\n",FinalSize);
	
//  vcg::tri::UpdateBounding<MyMesh>::Box(mesh);
  
//  // decimator initialization
//  vcg::LocalOptimization<MyMesh> DeciSession(mesh,&qparams);
	
//  int t1=clock();
//  DeciSession.Init<MyTriEdgeCollapse>();
//  int t2=clock();
//  printf("Initial Heap Size %i\n",int(DeciSession.h.size()));

//  DeciSession.SetTargetSimplices(FinalSize);
//  DeciSession.SetTimeBudget(0.5f);
//  if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

//  while(DeciSession.DoOptimization() && mesh.fn>FinalSize && DeciSession.currMetric < TargetError)
//    printf("Current Mesh size %7i heap sz %9i err %9g \r",mesh.fn, int(DeciSession.h.size()),DeciSession.currMetric);

//  int t3=clock();
//  printf("mesh  %d %d Error %g \n",mesh.vn,mesh.fn,DeciSession.currMetric);
//  printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
	
//  vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh,argv[2]);
	return 0;

}
