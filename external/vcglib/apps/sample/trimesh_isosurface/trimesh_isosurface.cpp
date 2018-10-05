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
#include <vcg/complex/complex.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace std;
using namespace vcg;

typedef float ScalarType;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
                                                                                Use<MyFace>			::AsFaceType>{};

class MyVertex     : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags>{};
class MyFace       : public Face< MyUsedTypes, face::VertexRef, face::BitFlags> {};

class MyMesh		: public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};



typedef SimpleVolume<SimpleVoxel<float> > MyVolume;

int main(int /*argc*/ , char **/*argv*/)
{
    MyVolume	volume;

  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
  typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
  MyWalker walker;


  // Simple initialization of the volume with some cool perlin noise
  vcg::Box3f bb(vcg::Point3f(-1,-1,-1),vcg::Point3f(1,1,1));
  volume.Init(Point3i(64,64,64),bb);
  for(int i=0;i<64;i++)
    for(int j=0;j<64;j++)
      for(int k=0;k<64;k++)
        volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);


	// MARCHING CUBES
	MyMesh		mc_mesh;
	printf("[MARCHING CUBES] Building mesh...");
	MyMarchingCubes					mc(mc_mesh, walker);
	walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, 20*20);
	vcg::tri::io::ExporterPLY<MyMesh>::Save( mc_mesh, "marching_cubes.ply");

	printf("OK!\n");
};
