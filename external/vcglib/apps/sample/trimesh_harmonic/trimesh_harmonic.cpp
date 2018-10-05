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
#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/harmonic.h>


using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::Normal3f, face::Mark, face::BitFlags, face::VFAdj, face::FFAdj > {};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace>   > {};




int main( int argc, char **argv )
{
  MyMesh m;
  if(argc < 2 )
  {
    printf("Usage: trimesh_harmonic mesh.ply\n");
    return -1;
  }

  int ret= tri::io::ImporterPLY<MyMesh>::Open(m,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  tri::UpdateTopology<MyMesh>::FaceFace(m);


  srand(time(0));
  int ind0=rand()%m.vn;
  int ind1=rand()%m.vn;
  printf("Computing harmonic field from vertex %i to vertex %i\n",ind0,ind1);

  // Get the two vertices with value set
  tri::Harmonic<MyMesh, double>::ConstraintVec constraints;
  constraints.push_back(tri::Harmonic<MyMesh, double>::Constraint(&(m.vert[ind0]), 1.0f));
  constraints.push_back(tri::Harmonic<MyMesh, double>::Constraint(&(m.vert[ind1]), 2.0f));

  MyMesh::PerVertexAttributeHandle<double> handle = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(m, "harmonic");
  bool ok = tri::Harmonic<MyMesh, double>::ComputeScalarField(m, constraints, handle);
  assert(ok);

  tri::UpdateQuality<MyMesh>::VertexFromAttributeHandle(m,handle);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m,1,2);
  tri::io::ExporterPLY<MyMesh>::Save(m,"harmonic.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);

  return 0;
}
