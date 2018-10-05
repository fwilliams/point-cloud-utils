/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2009                                           \/)\/    *
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

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/vertex/component.h>

#include <vcg/complex/used_types.h>

#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component.h>
#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/parametrization/poisson_solver.h>
#include<vcg/complex/algorithms/update/texture.h>
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::TexCoord2f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::VertexRef, face::BitFlags, face::FFAdj , face::WedgeTexCoord2f> {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};


int main( int argc, char **argv )
{
  MyMesh m;
  if(argc < 2 ) return -1;

  printf("Reading %s\n",argv[1]);
  int ret= tri::io::ImporterPLY<MyMesh>::Open(m,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  printf("Mesh has %i vn %i fn\n",m.VN(),m.FN());
  tri::PoissonSolver<MyMesh> PS(m);

  if(!PS.IsFeaseable())
  {
    printf("mesh is not homeomorphic to a disk\n");
    return -1;
  } else
    printf("OK - mesh is homeomorphic to a disk\n");

  PS.Init();
  PS.SetBorderAsFixed();
  PS.FixDefaultVertices();
  PS.SolvePoisson(true);

  tri::UpdateTexture<MyMesh>::WedgeTexFromVertexTex(m);
  tri::io::ExporterPLY<MyMesh>::Save(m,"pippo.ply",tri::io::Mask::IOM_WEDGTEXCOORD);

  return 0;
}
