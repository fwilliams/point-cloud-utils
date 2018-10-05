#include <vector>

#include <vcg/simplex/vertex/base.h>   
#include <vcg/simplex/vertex/component.h>   
#include <vcg/simplex/face/base.h>   
#include <vcg/simplex/face/component.h>   

#include <vcg/complex/complex.h>   
#include<vcg/container/simple_temporary_data.h>

#include<vcg/complex/algorithms/create/platonic.h>

class MyEdge;
class MyFace;

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::Coord3d, vcg::vert::Normal3f>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace, vcg::face::VertexRef>{};

class MyMesh: public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};


int main()
{
	MyMesh m;
	vcg::tri::Tetrahedron(m);
	vcg::SimpleTempData<MyMesh::VertContainer, short> MyTempData(m.vert);

	MyTempData.Start();         // enable the user defined attribute (memory is allocated)

	MyMesh::VertexIterator vi;  // declare the iterator over the vertices
	for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
	{
		MyTempData[*vi]    = 10;    // assign the value for the 'short' attribute
		MyTempData[vi]     = 10;    // you can pass the element or an iterator to it
	}

	MyTempData.Stop();          // disable the user defined attribute (memory is freed)

	return 0;
}