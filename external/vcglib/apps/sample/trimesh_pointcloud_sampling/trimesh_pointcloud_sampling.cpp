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

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/clustering.h>

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
  if(argc<3)
  {
    printf("Usage trimesh_base <meshfilename> radius (as perc) (\n");
    return -1;
  }

  MyMesh m;
  MyMesh subM;
  MyMesh cluM;
  MyMesh rndM;

  tri::MeshSampler<MyMesh> mps(subM);
  tri::MeshSampler<MyMesh> mrs(rndM);

  if(tri::io::Importer<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> >::SamplingRandomGenerator().initialize(time(0));
  float perc = atof(argv[2]);

  float radius = m.bbox.Diag() * perc;
  printf("Subsampling a PointCloud of %i vert with %f radius\n",m.VN(),radius);
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam pp;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskParam::Stat pds; pp.pds=&pds;
  pp.bestSampleChoiceFlag=false;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::PoissonDiskPruning(mps, m, radius, pp);
  tri::io::ExporterPLY<MyMesh>::Save(subM,"PoissonMesh.ply");
  printf("Sampled %i vertices in %5.2f\n",subM.VN(), float(pds.pruneTime+pds.gridTime)/float(CLOCKS_PER_SEC));

  int t0=clock();
  tri::Clustering<MyMesh, vcg::tri::AverageColorCell<MyMesh> > ClusteringGrid;
  ClusteringGrid.Init(m.bbox,100000,radius*sqrt(2.0f));
  ClusteringGrid.AddPointSet(m);
  ClusteringGrid.ExtractMesh(cluM);
  int t1=clock();
  tri::io::ExporterPLY<MyMesh>::Save(cluM,"ClusterMesh.ply");
  printf("Sampled %i vertices in %5.2f\n",cluM.VN(), float(t1-t0)/float(CLOCKS_PER_SEC));

  int t2=clock();
  int sampleNum = (cluM.VN()+subM.VN())/2;
  tri::SurfaceSampling<MyMesh,tri::MeshSampler<MyMesh> >::VertexUniform(m, mrs,sampleNum);
  int t3=clock();
  tri::io::ExporterPLY<MyMesh>::Save(rndM,"RandomMesh.ply");
  printf("Sampled %i vertices in %5.2f\n",rndM.VN(), float(t3-t2)/float(CLOCKS_PER_SEC));


  return 0;
}
