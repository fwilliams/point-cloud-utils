#include <simplex\vertex\with\a.h>
#include <simplex\tetrahedron\with\avq.h>
#include <simplex\face\with\fa.h>
#include <complex\tetramesh\base.h>
#include <complex\tetramesh\update\topology.h>
#include <complex\tetramesh\update\allocate.h>
#include<vector>
#include<complex\tetramesh\update\triconvert.h>
#include<complex\trimesh\base.h>
using namespace vcg;
using namespace tetra;
using namespace tri;
#include<apps\test\tetramesh\test\myclasses.h>

int main( int argc, char **argv )
{
 MyMesh tm=MyMesh();
 UpdateTopology<vector<MyVertex>,vector<MyTetrahedron> > ut= UpdateTopology<vector<MyVertex>,vector<MyTetrahedron> >();
 Allocator<MyMesh> All= Allocator<MyMesh>();
 vector<MyVertex **> local_var=vector<MyVertex **>();
 char* filename="sphere.ts";
 tm.LoadTs(filename,1);
 ut.TTTopology(tm.vert,tm.tetra);
 ut.TestTTTopology(tm.vert,tm.tetra);
 ut.VTTopology(tm.vert,tm.tetra);
 All.AddVertices(tm,10,local_var);
 All.AddVertices(tm,10);
 ut.TestTTTopology(tm.vert,tm.tetra);
 MyTriMesh mesh;
 TriConverter <MyMesh,MyTriMesh>tric=TriConverter<MyMesh,MyTriMesh>();
 tric.Convert(tm,mesh);
}