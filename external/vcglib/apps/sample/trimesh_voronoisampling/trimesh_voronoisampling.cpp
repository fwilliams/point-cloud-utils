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
class MyFace    : public Face< MyUsedTypes,   face::VertexRef, face::Normal3f, face::BitFlags, face::Mark, face::VFAdj, face::FFAdj > {};
class MyEdge    : public Edge< MyUsedTypes, edge::VertexRef, edge::BitFlags>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyEdge>, vector<MyFace>   > {};

class EmEdge;
class EmFace;
class EmVertex;
struct EmUsedTypes : public UsedTypes<	Use<EmVertex>   ::AsVertexType,
                                        Use<EmEdge>     ::AsEdgeType,
                                        Use<EmFace>     ::AsFaceType>{};

class EmVertex  : public Vertex<EmUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::VFAdj , vertex::Qualityf, vertex::Color4b, vertex::BitFlags  >{};
class EmFace    : public Face< EmUsedTypes,   face::VertexRef, face::BitFlags, face::VFAdj > {};
class EmEdge    : public Edge< EmUsedTypes, edge::VertexRef> {};
class EmMesh    : public tri::TriMesh< vector<EmVertex>, vector<EmEdge>, vector<EmFace>   > {};


int main( int argc, char **argv )
{
  MyMesh baseMesh;
  if(argc < 4 )
  {
    printf("Usage trimesh_voronoisampling mesh sampleNum iterNum\n");
     return -1;
  }
  int sampleNum = atoi(argv[2]);
  int iterNum   = atoi(argv[3]);

  bool fixCornerFlag=true;
  bool uniformEdgeSamplingFlag = true;

  printf("Reading %s and sampling %i points with %i iteration\n",argv[1],sampleNum,iterNum);
  int ret= tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
  if(ret!=0)
  {
    printf("Unable to open %s for '%s'\n",argv[1],tri::io::ImporterPLY<MyMesh>::ErrorMsg(ret));
    return -1;
  }

  tri::VoronoiProcessingParameter vpp;

  tri::Clean<MyMesh>::RemoveUnreferencedVertex(baseMesh);
  tri::Allocator<MyMesh>::CompactEveryVector(baseMesh);
  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);

  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskParam pp;
  float radius = tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::ComputePoissonDiskRadius(baseMesh,sampleNum);
  tri::VoronoiProcessing<MyMesh>::PreprocessForVoronoi(baseMesh,radius,vpp);

  tri::UpdateFlags<MyMesh>::FaceBorderFromVF(baseMesh);
  tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(baseMesh);


  // -- Build a sampling with just corners (Poisson filtered)
  MyMesh poissonCornerMesh;
  std::vector<Point3f> sampleVec;
  tri::TrivialSampler<MyMesh> mps(sampleVec);
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::VertexBorderCorner(baseMesh,mps,math::ToRad(150.f));
  tri::BuildMeshFromCoordVector(poissonCornerMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(poissonCornerMesh,"cornerMesh.ply");
  sampleVec.clear();
  MyMesh borderMesh,poissonBorderMesh;


  if(uniformEdgeSamplingFlag)
  {

  }
  else
  {
    if(fixCornerFlag)
    {
      tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, poissonCornerMesh, radius, pp);
      tri::BuildMeshFromCoordVector(poissonCornerMesh,sampleVec);
      tri::io::ExporterPLY<MyMesh>::Save(poissonCornerMesh,"poissonCornerMesh.ply");
      // Now save the corner as Fixed Seeds for later...
      std::vector<MyVertex *> fixedSeedVec;
      tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,sampleVec,fixedSeedVec);
      tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::FixVertexVector(baseMesh,fixedSeedVec);
      vpp.preserveFixedSeed=true;
    }

  // -- Build a sampling with points on the border
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::VertexBorder(baseMesh,mps);
  tri::BuildMeshFromCoordVector(borderMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(borderMesh,"borderMesh.ply");

  // -- and then prune the border sampling with poisson strategy using the precomputed corner vertexes.
  pp.preGenMesh = &poissonCornerMesh;
  pp.preGenFlag=true;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, borderMesh, radius*0.8f, pp);
  tri::BuildMeshFromCoordVector(poissonBorderMesh,sampleVec);
  }

  tri::io::ExporterPLY<MyMesh>::Save(poissonBorderMesh,"PoissonEdgeMesh.ply");

  // -- Build the montercarlo sampling of the surface
  MyMesh MontecarloSurfaceMesh;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::Montecarlo(baseMesh,mps,50000);
  tri::BuildMeshFromCoordVector(MontecarloSurfaceMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(MontecarloSurfaceMesh,"MontecarloSurfaceMesh.ply");

  // -- Prune the montecarlo sampling with poisson strategy using the precomputed vertexes on the border.
  pp.preGenMesh = &poissonBorderMesh;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, MontecarloSurfaceMesh, radius, pp);
  MyMesh PoissonMesh;
  tri::BuildMeshFromCoordVector(PoissonMesh,sampleVec);
  tri::io::ExporterPLY<MyMesh>::Save(PoissonMesh,"PoissonMesh.ply");

  std::vector<MyVertex *> seedVec;
  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,sampleVec,seedVec);

  // Select all the vertexes on the border to define a constrained domain.
  // In our case we select the border vertexes to make sure that the seeds on the border
  // relax themselves remaining on the border
  for(size_t i=0;i<baseMesh.vert.size();++i){
    if(baseMesh.vert[i].IsB())
        baseMesh.vert[i].SetS();
  }

//  vpp.deleteUnreachedRegionFlag=true;
  vpp.deleteUnreachedRegionFlag=false;
  vpp.triangulateRegion = false;
  vpp.geodesicRelaxFlag=false;
  vpp.constrainSelectedSeed=true;

  tri::EuclideanDistance<MyMesh> dd;
  int t0=clock();
  // And now, at last, the relaxing procedure!
  int actualIter = tri::VoronoiProcessing<MyMesh, tri::EuclideanDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum, dd, vpp);
  int t1=clock();

  MyMesh voroMesh, voroPoly, delaMesh;
  // Get the result in some pleasant form converting it to a real voronoi diagram.
  if(tri::VoronoiProcessing<MyMesh>::CheckVoronoiTopology(baseMesh,seedVec))
     tri::VoronoiProcessing<MyMesh>::ConvertVoronoiDiagramToMesh(baseMesh,voroMesh,voroPoly,seedVec, vpp);
  else
    printf("WARNING some voronoi region are not disk like; the resulting delaunay triangulation is not manifold.\n");

  tri::io::ExporterPLY<MyMesh>::Save(baseMesh,"base.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );
  tri::io::ExporterPLY<MyMesh>::Save(voroMesh,"voroMesh.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_FLAGS );
  tri::io::ExporterPLY<MyMesh>::Save(voroPoly,"voroPoly.ply",tri::io::Mask::IOM_VERTCOLOR| tri::io::Mask::IOM_EDGEINDEX ,false);

  tri::VoronoiProcessing<MyMesh>::ConvertDelaunayTriangulationToMesh(baseMesh,delaMesh, seedVec);
  tri::io::ExporterPLY<MyMesh>::Save(delaMesh,"delaMesh.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );
  tri::VoronoiProcessing<MyMesh>::RelaxRefineTriangulationSpring(baseMesh,delaMesh,2,10);
  tri::io::ExporterPLY<MyMesh>::Save(delaMesh,"delaMeshRef.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );


//  tri::io::ImporterPLY<MyMesh>::Open(baseMesh,argv[1]);
//  tri::UpdateTopology<MyMesh>::VertexFace(baseMesh);
//  tri::PoissonSampling<MyMesh>(baseMesh,pointVec,sampleNum,radius,radiusVariance);
//  tri::VoronoiProcessing<MyMesh>::SeedToVertexConversion(baseMesh,pointVec,seedVec);
//  tri::IsotropicDistance<MyMesh> id(baseMesh,radiusVariance);
//  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::VoronoiRelaxing(baseMesh, seedVec, iterNum,id,vpp);
//  tri::VoronoiProcessing<MyMesh, tri::IsotropicDistance<MyMesh> >::ConvertVoronoiDiagramToMesh(baseMesh,outMesh,polyMesh,seedVec, id, vpp);

//  tri::io::ExporterPLY<MyMesh>::Save(outMesh,"outW.ply",tri::io::Mask::IOM_VERTCOLOR );
//  tri::io::ExporterPLY<MyMesh>::Save(polyMesh,"polyW.ply",tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_EDGEINDEX,false);
//  tri::io::ExporterDXF<MyMesh>::Save(polyMesh,"outW.dxf");
  printf("Completed! %i (%i) iterations in %f sec for %lu seeds \n",actualIter, iterNum,float(t1-t0)/CLOCKS_PER_SEC,seedVec.size());
  return 0;
}
