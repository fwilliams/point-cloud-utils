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

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::Coord3d, vcg::vert::Normal3f>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace, vcg::face::VertexRef,vcg::face::FFAdj>{};

class MyMesh: public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

void OneRingNeighborhood( MyFace * f)
{
	MyVertex * v = f->V(0);  
	MyFace* start = f;
	vcg::face::Pos<MyFace> p(f,0,v);// constructor that takes face, edge and vertex
	do
	{
		p.FlipF();
		p.FlipE();
	}while(p.f!=start);
}

#include <vcg/simplex/face/jumping_pos.h> // include the definition of jumping pos

void OneRingNeighborhoodJP( MyFace * f)
{
	MyVertex * v = f->V(0);  
	MyFace* start = f;
	vcg::face::JumpingPos<MyFace> p(f,0,v);// constructor that takes face, edge and vertex
	do
	{
		p.NextFE();
	}while(p.f!=start);
}

int main()
{
	MyMesh m;
	vcg::tri::Tetrahedron(m);
	vcg::tri::UpdateTopology<MyVCGMesh>::FaceFace(m);
	OneRingNeighborhood(&(*m.face.begin()));
	OneRingNeighborhoodJP(&(*m.face.begin()));
	return 0;
}