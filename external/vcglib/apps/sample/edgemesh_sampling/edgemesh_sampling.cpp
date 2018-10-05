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

#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>
#include<wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VEAdj, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes,edge::VertexRef, edge::VEAdj, edge::EEAdj, edge::BitFlags >{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};


int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage edgemesh_sampling <meshfilename.off> radius\n");
    return -1;
  }

  MyMesh m,e;

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  printf("Input mesh has %i vert %i faces\n",m.vn,m.fn);


  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m);

  tri::UpdateFlags<MyMesh>::FaceFauxSignedCrease(m,math::ToRad(-40.f),math::ToRad(40.f  ),false);
  assert(tri::Clean<MyMesh>::IsFaceFauxConsistent(m));
  tri::BuildFromNonFaux(m,e);
  printf("Out mesh has %i vert %i edges\n",e.vn,e.en);

  tri::UpdateTopology<MyMesh>::VertexEdge(e);
  tri::Clean<MyMesh>::SelectNonManifoldVertexOnEdgeMesh(e);
  printf("Selected vertices %i\n",tri::UpdateSelection<MyMesh>::VertexCount(e));
  tri::Clean<MyMesh>::SplitSelectedVertexOnEdgeMesh(e);
  printf("Out mesh has %i vert %i edges\n",e.vn,e.en);
  tri::Clean<MyMesh>::SelectCreaseVertexOnEdgeMesh(e,math::ToRad(30.f));
  printf("Selected vertices %i\n",tri::UpdateSelection<MyMesh>::VertexCount(e));
  tri::Clean<MyMesh>::SplitSelectedVertexOnEdgeMesh(e);
  printf("Out mesh has %i vert %i edges\n",e.vn,e.en);

  tri::io::ExporterPLY<MyMesh>::Save(e,"out.ply", tri::io::Mask::IOM_EDGEINDEX);

  std::vector<Point3f> sampleVec;
  tri::TrivialSampler<MyMesh> ps(sampleVec);
  tri::SurfaceSampling<MyMesh>::EdgeMeshUniform(e,ps,m.bbox.Diag()/90.0f);
  MyMesh sampleMesh;
  tri::BuildMeshFromCoordVector(sampleMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(sampleMesh,"sampleMesh.ply");
  return 0;
}
