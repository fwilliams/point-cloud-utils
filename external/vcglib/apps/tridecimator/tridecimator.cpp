// stuff to define the mesh
#include <vcg/complex/complex.h>

// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace vcg;
using namespace tri;

/**********************************************************
Mesh Classes for Quadric Edge collapse based simplification

For edge collpases we need verteses with:
- V->F adjacency
- per vertex incremental mark
- per vertex Normal


Moreover for using a quadric based collapse the vertex class
must have also a Quadric member Q();
Otherwise the user have to provide an helper function object
to recover the quadric.

******************************************************/
// The class prototypes.
class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType,Use<MyEdge>::AsEdgeType,Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes,
    vertex::VFAdj,
    vertex::Coord3f,
    vertex::Normal3f,
    vertex::Mark,
    vertex::Qualityf,
    vertex::BitFlags  >{
public:
  vcg::math::Quadric<double> &Qd() {return q;}
private:
  math::Quadric<double> q;
  };

class MyEdge : public Edge< MyUsedTypes> {};

typedef BasicVertexPair<MyVertex> VertexPair;

class MyFace    : public Face< MyUsedTypes,
  face::VFAdj,
  face::VertexRef,
  face::BitFlags > {};

// the main mesh class
class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};


class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > {
            public:
            typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh,  VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > TECQ;
            typedef  MyMesh::VertexType::EdgeType EdgeType;
            inline MyTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
};

void Usage()
{
    printf(
          "---------------------------------\n"
          "        TriDecimator 1.0 \n"
          "     http://vcg.isti.cnr.it\n"
          "   release date: " __DATE__ 
          "\n---------------------------------\n\n"
          "Copyright 2003-2016 Visual Computing Lab I.S.T.I. C.N.R.\n"
          "\nUsage:  "\
          "tridecimator fileIn fileOut face_num [opt]\n"\
          "Where opt can be:\n"\
          "     -e# QuadricError threshold  (range [0,inf) default inf)\n"
          "     -b# Boundary Weight (default .5)\n"
          "     -q# Quality threshold (range [0.0, 0.866],  default .3 )\n"
          "     -n# Normal threshold  (degree range [0,180] default 90)\n"
          "     -w# Quality weight factor  (10)\n"
          "     -E# Minimal admitted quadric value (default 1e-15, must be >0)\n"
          "     -Q[y|n]  Use or not Face Quality Threshold (default yes)\n"
          "     -P[y|n]  Add or not QualityQuadric (default no)\n"
          "     -N[y|n]  Use or not Face Normal Threshold (default no)\n"
          "     -A[y|n]  Use or not Area Weighted Quadrics (default yes)\n"
          "     -O[y|n]  Use or not vertex optimal placement (default yes)\n"
          "     -S[y|n]  Use or not Scale Independent quadric measure(default yes) \n"
          "     -B[y|n]  Preserve or not mesh boundary (default no)\n"
          "     -T[y|n]  Preserve or not Topology (default no)\n"
          "     -W[y|n]  Use or not per vertex Quality to weight the quadric error (default no)\n"
          "     -C       Before simplification, remove duplicate & unreferenced vertices\n"
          );
  exit(-1);
}

int main(int argc ,char**argv)
{
  if(argc<4) Usage();

  MyMesh mesh;
  
  int FinalSize=atoi(argv[3]);
  int err=vcg::tri::io::Importer<MyMesh>::Open(mesh,argv[1]);
  if(err)
  {
    printf("Unable to open mesh %s : '%s'\n",argv[1],vcg::tri::io::Importer<MyMesh>::ErrorMsg(err));
    exit(-1);
  }
  printf("mesh loaded %d %d \n",mesh.vn,mesh.fn);

  TriEdgeCollapseQuadricParameter qparams;
  qparams.QualityThr  =.3;
  float TargetError=std::numeric_limits<float>::max();
  bool CleaningFlag =false;
     // parse command line.
    for(int i=4; i < argc;)
    {
      if(argv[i][0]=='-')
        switch(argv[i][1])
      {
        case 'Q' : if(argv[i][2]=='y') { qparams.QualityCheck	= true;  printf("Using Quality Checking\n");	}
                                  else { qparams.QualityCheck	= false; printf("NOT Using Quality Checking\n");	}                break;
        case 'N' : if(argv[i][2]=='y') { qparams.NormalCheck	= true;  printf("Using Normal Deviation Checking\n");	}
                                  else { qparams.NormalCheck	= false; printf("NOT Using Normal Deviation Checking\n");	}        break;
        case 'O' : if(argv[i][2]=='y') { qparams.OptimalPlacement	= true;  printf("Using OptimalPlacement\n");	}
                                  else { qparams.OptimalPlacement	= false; printf("NOT Using OptimalPlacement\n");	}        break;
        case 'S' : if(argv[i][2]=='y') { qparams.ScaleIndependent	= true;  printf("Using ScaleIndependent\n");	}
                                  else { qparams.ScaleIndependent	= false; printf("NOT Using ScaleIndependent\n");	}        break;
        case 'B' : if(argv[i][2]=='y') { qparams.PreserveBoundary	= true;  printf("Preserving Boundary\n");	}
                                  else { qparams.PreserveBoundary	= false; printf("NOT Preserving Boundary\n");	}        break;
        case 'T' : if(argv[i][2]=='y') { qparams.PreserveTopology	= true;  printf("Preserving Topology\n");	}
                                  else { qparams.PreserveTopology	= false; printf("NOT Preserving Topology\n");	}        break;
        case 'P' : if(argv[i][2]=='y') { qparams.QualityQuadric 	= true;  printf("Adding Quality Quadrics\n");	}
                                  else { qparams.QualityQuadric 	= false; printf("NOT Adding Quality Quadrics\n");	}        break;
        case 'W' : if(argv[i][2]=='y') { qparams.QualityWeight 	= true;  printf("Using per vertex Quality as Weight\n");	}
                                  else { qparams.QualityWeight 	= false; printf("NOT Using per vertex Quality as Weight\n");	}        break;
        case 'w' :	qparams.QualityWeightFactor	= atof(argv[i]+2);	           printf("Setting Quality Weight factor to %f\n",atof(argv[i]+2)); 	 break;
        case 'q' :	qparams.QualityThr	= atof(argv[i]+2);	           printf("Setting Quality Thr to %f\n",atof(argv[i]+2)); 	 break;
        case 'n' :	qparams.NormalThrRad = math::ToRad(atof(argv[i]+2));  printf("Setting Normal Thr to %f deg\n",atof(argv[i]+2)); break;
        case 'b' :	qparams.BoundaryWeight  = atof(argv[i]+2);			printf("Setting Boundary Weight to %f\n",atof(argv[i]+2)); break;
        case 'e' :	TargetError = float(atof(argv[i]+2));			printf("Setting TargetError to %g\n",atof(argv[i]+2)); break;
        case 'C' :	CleaningFlag=true;  printf("Cleaning mesh before simplification\n"); break;

        default  :  printf("Unknown option '%s'\n", argv[i]);
          exit(0);
      }
      i++;
    }



  if(CleaningFlag){
      int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
      int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
      printf("Removed %i duplicate and %i unreferenced vertices from mesh \n",dup,unref);
  }


  printf("reducing it to %i\n",FinalSize);

  vcg::tri::UpdateBounding<MyMesh>::Box(mesh);

  // decimator initialization
  vcg::LocalOptimization<MyMesh> DeciSession(mesh,&qparams);

  int t1=clock();
  DeciSession.Init<MyTriEdgeCollapse>();
  int t2=clock();
  printf("Initial Heap Size %i\n",int(DeciSession.h.size()));

  DeciSession.SetTargetSimplices(FinalSize);
  DeciSession.SetTimeBudget(0.5f);
  if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

  while(DeciSession.DoOptimization() && mesh.fn>FinalSize && DeciSession.currMetric < TargetError)
    printf("Current Mesh size %7i heap sz %9i err %9g \r",mesh.fn, int(DeciSession.h.size()),DeciSession.currMetric);

  int t3=clock();
  printf("mesh  %d %d Error %g \n",mesh.vn,mesh.fn,DeciSession.currMetric);
  printf("\nCompleted in (%5.3f+%5.3f) sec\n",float(t2-t1)/CLOCKS_PER_SEC,float(t3-t2)/CLOCKS_PER_SEC);

  vcg::tri::io::ExporterPLY<MyMesh>::Save(mesh,argv[2]);
    return 0;

}
