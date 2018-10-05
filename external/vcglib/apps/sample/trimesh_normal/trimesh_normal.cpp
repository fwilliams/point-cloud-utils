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
/*! \file trimesh_normal.cpp
\ingroup code_sample

\brief An example of all the methods for computing normals over a mesh.

*/
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import_off.h>

#include<vcg/complex/algorithms/update/normal.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.obj>\n");
    return -1;
  }

  MyMesh m;

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  Matrix44f m44 = Matrix44f::Identity();

  tri::UpdateNormal<MyMesh>::PerVertexClear(m);
  tri::UpdateNormal<MyMesh>::PerVertex(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::PerVertexNelsonMaxWeighted(m);
  tri::UpdateNormal<MyMesh>::PerFace(m);
  tri::UpdateNormal<MyMesh>::PerVertexFromCurrentFaceNormal(m);
  tri::UpdateNormal<MyMesh>::PerFaceFromCurrentVertexNormal(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  tri::UpdateNormal<MyMesh>::NormalizePerFace(m);
  tri::UpdateNormal<MyMesh>::NormalizePerFaceByArea(m);
  tri::UpdateNormal<MyMesh>::PerBitQuadFaceNormalized(m);
  tri::UpdateNormal<MyMesh>::PerVertexMatrix(m,m44);
  tri::UpdateNormal<MyMesh>::PerFaceMatrix(m,m44);
  tri::UpdateNormal<MyMesh>::PerWedgeCrease(m, math::ToRad(45.0f));

  return 0;
}
