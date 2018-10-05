#include <math.h>
#include <stdio.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/cylinder_clipping.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType >{};

class MyVertex : public Vertex<MyUsedTypes, vertex::Normal3f, vertex::Coord3f, vertex::BitFlags, vertex::Color4b >{};
class MyFace   : public Face<MyUsedTypes, face::Mark, face::Normal3f, face::FFAdj, face::BitFlags, face::VertexRef, face::Color4b > {};
class MyEdge   : public Edge<MyUsedTypes, edge::BitFlags>{};
class MyMesh   : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge> > {};

int main()
{
    MyMesh m;
    tri::Hexahedron(m);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::UpdateTopology<MyMesh>::FaceFace(m);    //for FFAdj
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    tri::UpdateFlags<MyMesh>::Clear(m);

    Point3f origin(0.8f,-0.4,0);
    Point3f end= origin+Point3f(0,1,0);
    float radius = 0.5f;
    MyMesh cm;
    tri::OrientedCylinder(cm,origin,end,radius,64,4);
    tri::io::ExporterPLY<MyMesh>::Save(cm,"cyl.ply");

    tri::CylinderClipping<MyMesh>::Apply(m,origin,end,radius);
    tri::CylinderClipping<MyMesh>::Apply(m,origin,end,radius/2.0f);

    tri::io::ExporterPLY<MyMesh>::Save(m,"cube.ply");
    return 0;
}
