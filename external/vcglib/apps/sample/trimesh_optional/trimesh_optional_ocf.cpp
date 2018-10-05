#include <vector>
#include <vcg/simplex/vertex/base.h>   
#include <vcg/simplex/vertex/component_occ.h>  
#include <vcg/simplex/face/base.h>   
#include <vcg/simplex/face/component.h>   

#include <vcg/complex/complex.h>   

#include<vcg/complex/algorithms/create/platonic.h>

class MyEdge;
class MyFace;

class MyVertex: public vcg::VertexSimp2<MyVertex,MyEdge,MyFace, vcg::vert::InfoOcf,vcg::vert::Coord3d, vcg::vert::Normal3fOcc>{};
class MyFace: public vcg::FaceSimp2<MyVertex,MyEdge,MyFace,vcg::face::VertexRef>{};
class MyMesh: public vcg::tri::TriMesh< vcg::vert::vector_occ<MyVertex>, std::vector<MyFace> >{};

int main()
{
 MyMesh m;
 vcg::tri::Tetrahedron(m);
 MyMesh::VertexIterator vi = m.vert.begin();

 (*vi).N() = vcg::Point3f(1.0,1.0,1.0); // ERROR   
 m.vert.EnableAttribute<vcg::vert::Normal3fOcc::NormalType>(); // this allocate the memory to store the normal
 (*vi).N() = vcg::Point3f(1.0,1.0,1.0); // OK
 m.vert.DisableAttribute<vcg::vert::Normal3fOcc::NormalType>(); // this deallocate the memory to store the normal

 (*vi).N() = vcg::Point3f(1.0,1.0,1.0); // ERROR  (again)! 
 return 0;
}
