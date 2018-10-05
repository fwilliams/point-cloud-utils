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
/*! \file trimesh_kdtree.cpp
\ingroup code_sample

\brief An example about using a kdtree to spatially index the vertexes of a mesh

KdTree are one of the Spatial indexing data structures available.
They are tailored for storing point-based structures and performing k-neighbours queries.
In this simple example we simply compute the average distance of a vertex from its neighbours.
\ref spatial_indexing for more Details
*/

#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/color.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,  face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2) argv[1]="../../meshes/torus_irregular.ply";

  MyMesh m;
  if(tri::io::Importer<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  VertexConstDataWrapper<MyMesh> ww(m);

  KdTree<float> tree(ww);
  KdTree<float>::PriorityQueue queue;

  for (int j = 0; j < m.VN(); j++) {
      tree.doQueryK(m.vert[j].cP(), 3, queue);
      int neighbours = queue.getNofElements();
      float avgDist=0;
      for (int i = 0; i < neighbours; i++) {
          int neightId = queue.getIndex(i);
          avgDist += Distance(m.vert[j].cP(),m.vert[neightId].cP());
      }
      m.vert[j].Q() = avgDist/=neighbours;
  }
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);
  tri::io::ExporterPLY<MyMesh>::Save(m,"out.ply",tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  return 0;
}
