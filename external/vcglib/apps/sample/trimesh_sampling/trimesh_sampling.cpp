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

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.obj> radius\n");
    return -1;
  }

  MyMesh m;

  if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::SamplingRandomGenerator().initialize(time(0));

  //----------------------------------------------------------------------
  // Basic Sample of a mesh surface
  // Build a point cloud with points with a plain poisson disk distribution

  int t0=clock();
  vector<Point3f> pointVec;
  float rad=0;
  if(argc>2) rad=atof(argv[2]);
  int sampleNum=rad?0:1000;
  tri::PoissonSampling<MyMesh>(m,pointVec,sampleNum,rad);
  int t1=clock();
  MyMesh BasicPoissonMesh;
  tri::BuildMeshFromCoordVector(BasicPoissonMesh,pointVec);

  tri::io::ExporterOFF<MyMesh>::Save(BasicPoissonMesh,"BasicPoissonMesh.off");
  printf("Computed a basic poisson disk distribution of %i vertices radius is %6.3f in %5.2f sec\n",BasicPoissonMesh.VN(),rad,float(t1-t0)/CLOCKS_PER_SEC);

  //----------------------------------------------------------------------
  // Advanced Sample
  // Make a feature dependent Poisson Disk sampling
  MyMesh MontecarloSurfaceMesh;
  MyMesh MontecarloEdgeMesh;
  MyMesh PoissonEdgeMesh;
  MyMesh PoissonMesh;

  std::vector<Point3f> sampleVec;
  tri::TrivialSampler<MyMesh> mps(sampleVec);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateNormal<MyMesh>::PerFace(m);
  tri::UpdateFlags<MyMesh>::FaceFauxCrease(m,math::ToRad(40.0f));
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::EdgeMontecarlo(m,mps,10000,false);
  tri::BuildMeshFromCoordVector(MontecarloEdgeMesh,sampleVec);
  tri::io::ExporterOFF<MyMesh>::Save(MontecarloEdgeMesh,"MontecarloEdgeMesh.off");

  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::VertexCrease(m, mps);
  tri::BuildMeshFromCoordVector(PoissonEdgeMesh,sampleVec);
  tri::io::ExporterOFF<MyMesh>::Save(PoissonEdgeMesh,"VertexCreaseMesh.off");

  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskParam pp;
  pp.preGenMesh = &PoissonEdgeMesh;
  pp.preGenFlag=true;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, MontecarloEdgeMesh, rad, pp);
  tri::BuildMeshFromCoordVector(PoissonEdgeMesh,sampleVec);
  tri::io::ExporterOFF<MyMesh>::Save(PoissonEdgeMesh,"PoissonEdgeMesh.off");

  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::Montecarlo(m,mps,50000);
  tri::BuildMeshFromCoordVector(MontecarloSurfaceMesh,sampleVec);
  tri::io::ExporterOFF<MyMesh>::Save(MontecarloSurfaceMesh,"MontecarloSurfaceMesh.off");

  pp.preGenMesh = &PoissonEdgeMesh;
  pp.preGenFlag=true;
  sampleVec.clear();
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::PoissonDiskPruning(mps, MontecarloSurfaceMesh, rad, pp);
  tri::BuildMeshFromCoordVector(PoissonMesh,sampleVec);
  tri::io::ExporterOFF<MyMesh>::Save(PoissonMesh,"PoissonMesh.off");
  printf("Computed a feature aware poisson disk distribution of %i vertices radius is %6.3f\n",PoissonMesh.VN(),rad);

  return 0;
}
