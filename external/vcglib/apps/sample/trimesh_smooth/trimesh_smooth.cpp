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
/*! \file trimesh_smooth.cpp
\ingroup code_sample

\brief the minimal example of using the lib

This file contain a minimal example of the library

*/

#include<vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>

#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/smooth.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>


using namespace vcg;
using namespace std;


class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh<vector<MyVertex>, vector<MyFace> > {};

int main(int argc,char ** argv)
{
  if(argc<4)
  {
    printf("Usage: trimesh_smooth <filename> <steps> <sigma> <fitstep>\n");
    return 0;
  }

  int Step= atoi(argv[2]);

  MyMesh m;

  //open a mesh
  int err = tri::io::Importer<MyMesh>::Open(m,argv[1]);
  if(err) { // all the importers return 0 in case of success
    printf("Error in reading %s: '%s'\n",argv[1], tri::io::Importer<MyMesh>::ErrorMsg(err));
    exit(-1);
  }

  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);

  tri::UpdateTopology<MyMesh>::VertexFace(m);

  for(int i=0;i<Step;++i)
  {
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
    tri::Smooth<MyMesh>::VertexCoordPasoDoble(m,atoi(argv[3]),atof(argv[4]),atoi(argv[5]));
  }

  //LaplacianSmooth(m,atoi(argv[2]));
  tri::io::ExporterPLY<MyMesh>::Save(m,"out.ply");

  return 0;
}

