/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file trimesh_allocate.cpp
\ingroup code_sample

\brief the minimal example about creating and deleting elements

Attributes are a simple mechanism to associate user-defined 'attributes' to the simplicies and to the mesh.
\ref attributes for more Details
*/

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>		::AsVertexType,
                                            vcg::Use<MyFace>			::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f,vcg::vertex::Normal3f>{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

int main()
{
  MyMesh m;
  MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(m,3);
  MyMesh::FaceIterator fi = vcg::tri::Allocator<MyMesh>::AddFaces(m,1);

  MyMesh::VertexPointer ivp[4];
  ivp[0]=&*vi; vi->P()=MyMesh::CoordType ( 0.0, 0.0, 0.0); ++vi;
  ivp[1]=&*vi; vi->P()=MyMesh::CoordType ( 1.0, 0.0, 0.0); ++vi;
  ivp[2]=&*vi; vi->P()=MyMesh::CoordType ( 0.0, 1.0, 0.0); ++vi;

  fi->V(0)=ivp[0];
  fi->V(1)=ivp[1];
  fi->V(2)=ivp[2];

  // Alternative, more compact, method for adding a single vertex
  ivp[3]= &*vcg::tri::Allocator<MyMesh>::AddVertex(m,MyMesh::CoordType ( 1.0, 1.0, 0.0));

  // Alternative, more compact, method for adding a single face (once you have the vertex pointers)
  vcg::tri::Allocator<MyMesh>::AddFace(m, ivp[1],ivp[0],ivp[3]);

  // a potentially dangerous pointer to a mesh element
  MyMesh::FacePointer fp = &m.face[0];
  vcg::tri::Allocator<MyMesh>::PointerUpdater<MyMesh::FacePointer> pu;

  // now the fp pointer could be no more valid due to eventual re-allocation of the m.face vector.
  vcg::tri::Allocator<MyMesh>::AddVertices(m,3);
  vcg::tri::Allocator<MyMesh>::AddFaces(m,1,pu);

  // check if an update of the pointer is needed and do it.
  if(pu.NeedUpdate()) pu.Update(fp);

  // Now fill the mesh with an Icosahedron and then delete some faces
  vcg::tri::Icosahedron(m);
  vcg::tri::Allocator<MyMesh>::DeleteFace(m,m.face[1]);
  vcg::tri::Allocator<MyMesh>::DeleteFace(m,m.face[3]);

  // If you loop in a mesh with deleted elements you have to skip them!
  for(fi = m.face.begin(); fi!=m.face.end(); ++fi )
  {
     if(!fi->IsD()) //    <---- Check added
       {
        MyMesh::CoordType b = vcg::Barycenter(*fi);
       }
  }

  // WRONG WAY of iterating: FN() != m.face.size() if there are deleted elements
  for(int i=0;i<m.FN();++i)
  {
     if(!fi->IsD())
       {
        MyMesh::CoordType b = vcg::Barycenter(*fi);
       }
  }

  // To remove the elements marked as deleted use
  vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
  vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);

  // To clean all the containers from deleted elements
  vcg::tri::Allocator<MyMesh>::CompactEveryVector(m);

  // finally lets copy this mesh onto another one.
  MyMesh m2;
  vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(m2,m);

}
