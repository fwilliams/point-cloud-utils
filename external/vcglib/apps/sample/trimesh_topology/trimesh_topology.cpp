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

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>

#include<vcg/complex/algorithms/update/topology.h>

#include<vcg/simplex/face/pos.h>

using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,face::FFAdj, face::Mark, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace > >{};


int main(int ,char ** )
{
  MyMesh m;

  //generate a mesh
  vcg::tri::Icosahedron(m);

  //update the face-face topology
  vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

  // Now for each face the FFp() FFi() members are correctly initialized 

  if(face::IsBorder(m.face[0],0)) printf("Edge 0 of face 0 is a border\n");
  else printf("Edge 0 of face 0 is NOT a border\n"); // always this path!

  vcg::face::FFDetach<MyFace>(m.face[0],0);  // Detach the face [0] from the mesh
  vcg::face::FFDetach<MyFace>(m.face[0],1);
  vcg::face::FFDetach<MyFace>(m.face[0],2);

  if(face::IsBorder(m.face[0],0)) printf("Edge 0 of face 0 is a border\n"); // always this path!
  else printf("Edge 0 of face 0 is NOT a border\n");

  tri::Allocator<MyMesh>::DeleteFace(m,m.face[0]);

  // declare an iterator on the mesh
  vcg::face::Pos<MyMesh::FaceType> he, hei;

  UnMarkAll(m);

  // Now a simple search and trace of all the borders of the mesh
  int BorderEdgeNum=0;
  int HoleNum=0;
  for(MyMesh::FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
  {
    for(int j=0;j<3;j++)
    {
      if ( face::IsBorder(*fi,j) && !tri::IsMarked(m,&*fi))
      {
        tri::Mark(m,&*fi);
        hei.Set(&*fi,j,fi->V(j));
        he=hei;
        do
        {
          BorderEdgeNum++;
          he.NextB(); // next pos along a border
          tri::Mark(m,he.f);
        }
        while (he.f!=hei.f);
        HoleNum++;
      }
    }
  }

  printf("Mesh has %i holes and %i border edges\n",HoleNum,BorderEdgeNum);
  return 0;
}

