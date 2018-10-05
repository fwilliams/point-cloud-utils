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
/*! \file trimesh_inertia.cpp
\ingroup code_sample

\brief An example of computing the inertia properties of meshes

Two meshes are created a rectangular box and a torus and their mass properties are computed and shown.
The result should match the closed formula for these objects (with a reasonable approximation)

*/

#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/implicit_smooth.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                            vcg::Use<MyEdge>     ::AsEdgeType,
                                            vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::FFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};


int main( int argc, char **argv )
{
  MyMesh m;
  if(argc < 2 )
  {
    printf("Usage: trimesh_implicit_smooth mesh.ply\n");
    return -1;
  }

  int ret= tri::io::ImporterPLY<MyMesh>::Open(m,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  tri::UpdateTopology<MyMesh>::FaceFace(m);

  ImplicitSmoother<MyMesh>::Parameter par;
  par.lambda = 0.5f;
  ImplicitSmoother<MyMesh>::Compute(m,par);

  tri::io::ExporterPLY<MyMesh>::Save(m,"smooth.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);

  return 0;
}
