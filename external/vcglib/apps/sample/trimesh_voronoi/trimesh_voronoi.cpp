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
#include <stdio.h>
#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<wrap/io_trimesh/import_ply.h>
#include<wrap/io_trimesh/export_off.h>
#include<wrap/io_trimesh/export_ply.h>
#include<wrap/io_trimesh/export_dxf.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/voronoi_processing.h>


using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::Normal3f, face::Mark, face::BitFlags, face::VFAdj, face::FFAdj > {};
class MyEdge    : public Edge< MyUsedTypes, edge::VertexRef, edge::BitFlags>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyEdge>, vector<MyFace>   > {};

int main( int argc, char **argv )
{
  MyMesh baseMesh,voronoiMesh, voronoiPoly, delaunayMesh;
  if(argc < 6 )
  {
    printf("Usage: trimesh_voronoi mesh [sampleNum] voronoiRelaxIter delaunayRefinementStep delaunayRelaxStep \n");
     return -1;
  }
  int sampleNum = atoi(argv[2]);
  int iterNum   = atoi(argv[3]);

  int refineStep   = atoi(argv[4]);
  int relaxStep   = atoi(argv[5]);

  int t0=clock();
  int ret= tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }
  tri::VoronoiProcessingParameter vpp;

  tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  int t1=clock();
  printf("Read         %30s (%7i vn %7i fn) in %6.3f \n",argv[1],baseMesh.vn,baseMesh.fn,float(t1-t0)/CLOCKS_PER_SEC);


  vector<Point3f> pointVec;
  float radius;
  tri::PoissonSampling<MyMesh>(baseMesh,pointVec,sampleNum,radius);
  vector<MyVertex *> seedVec;
  tri::VoronoiProcessing<MyMesh>::PreprocessForVoronoi(baseMesh,radius,vpp);
  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,pointVec,seedVec);

  int t2=clock();
  printf("Preprocessed %30s (%7i vn %7i fn) Computed   %i seed (asked %i) (radius %f) in %6.3f\n",
         argv[1],baseMesh.vn,baseMesh.fn, seedVec.size(), sampleNum, radius, float(t2-t1)/CLOCKS_PER_SEC);

  tri::EuclideanDistance<MyMesh> df;
  vpp.geodesicRelaxFlag=false;
  int actualIter = tri::VoronoiProcessing<MyMesh>::VoronoiRelaxing(baseMesh, seedVec, iterNum, df, vpp);

  int t3=clock();
  printf("relaxed %lu seeds for %i(up to %i) iterations in %f secs\n",
         seedVec.size(), actualIter, iterNum,float(t3-t2)/CLOCKS_PER_SEC);

  tri::io::ExporterPLY<MyMesh>::Save(baseMesh,"baseMesh.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY );
  if(tri::VoronoiProcessing<MyMesh>::CheckVoronoiTopology(baseMesh,seedVec))
  {
     tri::VoronoiProcessing<MyMesh>::ConvertVoronoiDiagramToMesh(baseMesh,voronoiMesh,voronoiPoly,seedVec, vpp);
  }
  else
  {
    printf("WARNING some voronoi region are not disk like; the resulting delaunay triangulation is not manifold.\n");
    refineStep=1;
  }

  tri::VoronoiProcessing<MyMesh>::ConvertDelaunayTriangulationToMesh(baseMesh,delaunayMesh,seedVec,true);
  tri::io::ExporterPLY<MyMesh>::Save(delaunayMesh,"delaunayBaseMesh.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTFLAGS,false );
  tri::VoronoiProcessing<MyMesh>::RelaxRefineTriangulationSpring(baseMesh,delaunayMesh,refineStep,relaxStep);

  int t4=clock();
  printf("Refined %i times and relaxed %i to a %i v %i f mesh in %f secs\n",
         refineStep, relaxStep, delaunayMesh.vn,delaunayMesh.fn,float(t4-t3)/CLOCKS_PER_SEC);

  tri::io::ExporterPLY<MyMesh>::Save(baseMesh,"baseMesh.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTQUALITY );
  tri::io::ExporterPLY<MyMesh>::Save(voronoiMesh,"voronoiMesh.ply",tri::io::Mask::IOM_VERTCOLOR );
  tri::io::ExporterPLY<MyMesh>::Save(delaunayMesh,"delaunayMesh.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTFLAGS,false );
  tri::io::ExporterPLY<MyMesh>::Save(voronoiPoly,"voronoiPoly.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_EDGEINDEX,false);
  return 0;
}
