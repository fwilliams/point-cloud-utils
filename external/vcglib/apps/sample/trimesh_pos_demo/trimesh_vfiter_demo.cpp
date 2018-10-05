#include <vector>

#include <vcg/simplex/vertex/base.h>   
#include <vcg/simplex/vertex/component.h>   
#include <vcg/simplex/face/base.h>   
#include <vcg/simplex/face/component.h>   

#include <vcg/complex/complex.h>   
#include<vcg/complex/algorithms/create/platonic.h>

#include<vcg/complex/algorithms/update/topology.h>

#include <vcg/simplex/face/pos.h> 

class MyEdge;
class MyFace;

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::Coord3d, vcg::vert::Normal3f,vcg::vert::VFAdj>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace, vcg::face::VertexRef,vcg::face::VFAdj>{};
class MyMesh: public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

void OneRingNeighborhoodVF( MyVertex * v)
{  
	vcg::face::VFIterator<MyFace> vfi(v); //initialize the iterator tohe first face
	for(;!vfi.End();++vfi)
	{
		MyFace* f = vfi.F();
		// ...do something with face f
	}
}

int main()
{
	MyMesh m;
	vcg::tri::Tetrahedron(m);
	vcg::tri::UpdateTopology<MyVCGMesh>::VertexFace(m);
	OneRingNeighborhoodVF(&(*m.vert.begin()));
	return 0;
}