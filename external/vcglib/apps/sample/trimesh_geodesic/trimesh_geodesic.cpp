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

#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_ply.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/geodesic.h>
#include<vcg/complex/algorithms/update/color.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Mark, vertex::VFAdj, vertex::Color4b, vertex::Qualityf, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::VFAdj,  face::VertexRef, face::Normal3f, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_geodesic <meshfilename.obj>\n");
//    return -1;
  }

  MyMesh m;

//  if(tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/disk_irregular_1k.ply")!=0)
    if(tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/disk_irregular_650k.ply")!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }

  Point3f c=m.bbox.Center();
  MyVertex*closest=&*m.vert.begin();
  float minDist = Distance(closest->P(),c);
  for(MyMesh::VertexIterator vi=m.vert.begin();vi!=m.vert.end(); ++vi)
  {
    if(Distance(vi->P(),c)<minDist)
    {
      minDist = Distance(vi->P(),c);
      closest = &*vi;
    }
  }
  vector<MyVertex*> seedVec;
  seedVec.push_back(closest);
  tri::EuclideanDistance<MyMesh> ed;
  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  tri::Allocator<MyMesh>::CompactEveryVector(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  tri::Geodesic<MyMesh>::Compute(m,seedVec,ed);
  pair<float,float> minmax = tri::Stat<MyMesh>::ComputePerVertexQualityMinMax(m);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);
  printf("min %f max %f\n",minmax.first,minmax.second);
  tri::io::ExporterPLY<MyMesh>::Save(m,"base.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);
  int t0=clock();
  tri::Geodesic<MyMesh>::PerVertexDijsktraCompute(m,seedVec,ed);
  int t1=clock();
  printf("Geodesic dijkstra %6.3f\n",float(t1-t0)/CLOCKS_PER_SEC);
  tri::UpdateColor<MyMesh>::PerVertexQualityRamp(m);
  tri::io::ExporterPLY<MyMesh>::Save(m,"base_d.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY);

  return 0;
}
