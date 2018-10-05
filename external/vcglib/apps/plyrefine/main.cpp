	
// mesh definition 
#include <vcg/simplex/vertex/with/afvn.h>
#include <vcg/simplex/face/with/af.h>
#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/refine.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

// std
#include <vector>

using namespace vcg;
using namespace std;

struct MyFace;
struct MyTetra;
struct MyEdge;
struct MyVertex: public VertexAFVNf<MyEdge,MyFace,MyTetra>{};
struct MyFace: public FaceAF<MyVertex,MyEdge,MyFace>{};
struct MyMesh: public tri::TriMesh< vector<MyVertex>, vector<MyFace> >{};




#define FLAT	0
#define ARC		1
#define BUTTERFLY 2
//#define BUTTERFLY2 3

//#define PLANE 4
//#define SPHERE 5

#define LENGTH 6
#define ONLY_SEL 7


int  main(int argc, char **argv){
		if(argc<4)
	{
		printf(
		"\n                  PlyRefine ("__DATE__")\n"
			"						Visual Computing Group I.S.T.I. C.N.R.\n"
			"Usage: PlyRefine filein.ply fileout.ply [command list]\n"
			"Commands: \n"
			" Refinement rules:\n"
			"     -m# midpoint flat \n"
			"     -a# midpoint arc\n"
			"     -b# butterfly\n"
			//"     -p# clip with plane   \n"
			//"     -s# clip with sphere   \n"
			" Selective Refinement\n"
			"     -L# refine only if the the edge is longer than #(default 0.0)\n"
			"     -S(0|1)  refine only selected faces\n"
			);
		exit(0);
	}

	typedef	pair<int,float> OP_TYPE;
	vector<OP_TYPE > operations;
	bool only_selected=false;
	int i=3; int n_steps;float length=0;
	while(i<argc)
		{
			if(argv[i][0]!='-')
				{printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			operations.push_back(OP_TYPE()); OP_TYPE & op = operations.back();
			switch(argv[i][1])
			{
				
				case 'm' :	
										n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
										op.first = FLAT;op.second = n_steps; break;
				case 'a' :	
										n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
										op.first = ARC;op.second = n_steps; break;
				case 'b' :	
										n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
										op.first = BUTTERFLY;op.second = n_steps; break;
				//case 'v' :	
				//						n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
				//						op.first = BUTTERFLY2;op.second = n_steps; break;
				//case 'p' :	
				//						n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
				//						op.first = PLANE;op.second = n_steps; break;
				//case 's' :	
				//						n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
				//						op.first = SPHERE;op.second = n_steps; break;
				case 'L' :	
										n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
										op.first = LENGTH;op.second = n_steps; break;
				case 'S' :	
										n_steps = atof(argv[i]+2); n_steps=max(1,n_steps);
										op.first = ONLY_SEL;op.second = n_steps; break;
				default : {printf("Error unable to parse option '%s'\n",argv[i]); exit(0);}
			}
			++i;
		}


	MyMesh m;
		if(vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
			{
			printf("Error reading file  %s\n",argv[1]);
			exit(0);
		}

		

	vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
	vcg::tri::UpdateTopology<MyMesh>::FaceBorderFlags(m);
	vcg::tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);

	int h;
	for(i=0;i < operations.size();++i){
					
				switch(operations[i].first){
					case FLAT:
										for(h=0;h<operations[i].second;++h){
											Refine(m,MidPoint<MyMesh>(),length,only_selected);
											}
									 break;
					case ARC:
								for(h=0;h<operations[i].second;++h){
											Refine(m,MidPointArc<MyMesh>(),length,only_selected);}
								 break;
					case BUTTERFLY:
								for(h=0;h<operations[i].second;++h){
											Refine(m,MidPointButterfly<MyMesh>(),length,only_selected); }
										break;
					//case BUTTERFLY2:
					//			for(h=0;h<operations[i].second;++h){
					//						m.ComputeNormalizedVertexNormal();
					//						Refine(m,MidPointButterfly2<MyMesh>(),length,only_selected);
					//				}
					//						break;
	/*				case PLANE:
								for(h=0;h<operations[i].second;++h){
										m.ComputeNormalizedVertexNormal();
										Refine(m,MidPointPlane<MyMesh>(),length,only_selected);
									}
									 break;
					case SPHERE:
								for(h=0;h<operations[i].second;++h){
											m.ComputeNormalizedVertexNormal();
											Refine(m,MidPointSphere<MyMesh>(),length,only_selected);
									}
							break;
	*/				case LENGTH:
										length = operations[i].second; break;
					case ONLY_SEL:
										only_selected = (bool)operations[i].second; break;
			}
		}
	
	
	//m.ComputeNormalizedVertexNormal();
	//Refine(m,MidPointArc<MyMesh>(),0);
	vcg::tri::io::PlyInfo pi;
	vcg::tri::io::ExporterPLY<MyMesh>::Save(m,argv[2],pi.mask);
	return 0;
	}
