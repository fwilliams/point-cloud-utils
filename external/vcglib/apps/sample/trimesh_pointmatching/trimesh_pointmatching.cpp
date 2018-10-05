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
/*! \file trimesh_fitting.cpp
\ingroup code_sample

\brief A small example about sampling and fitting

Given a mesh (an icosahedron) for each face we get a few random samples over it, and then we recover:
- the plane fitting them (that coincide with the face plane and exactly approximate all the sample points)
- the plane fitting the perturbed version of this set
- the plane fitting the perturbed version of the set but using a weighted fitting scheme.
*/

#include<vcg/complex/complex.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<wrap/io_trimesh/import_off.h>
#include<vcg/space/point_matching.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                        vcg::Use<MyEdge>     ::AsEdgeType,
                                        vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

float EvalPlane(vcg::Plane3f &pl, std::vector<vcg::Point3f> posVec)
{
  float off=0;
  for(size_t i=0;i<posVec.size();++i)
    off += fabs(vcg::SignedDistancePlanePoint(pl,posVec[i]));

  off/=float(posVec.size());
  return off;
}


int main( )
{
  MyMesh m;
  vcg::tri::Icosahedron(m);
  vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m);
  vcg::tri::UpdateBounding<MyMesh>::Box(m);

  // As a simple test we get a few random points on a mesh,
  // we rot and trans them
  // and we fit them

  std::vector<vcg::Point3f> ExactVec;
  std::vector<vcg::Point3f> PerturbVec;
  tri::MontecarloSampling(m,ExactVec,10);
  PerturbVec=ExactVec;

  Matrix44f RotM;
  Matrix44f TraM;
  Point3f dir;
  vcg::math::MarsenneTwisterRNG rnd;

  vcg::math::GeneratePointInUnitBallUniform<float>(rnd);
  RotM.SetRotateDeg(rand()%360,dir);
  TraM.SetTranslate(1,2,3);
  Matrix44f RigidM = RotM*TraM;

  for(size_t i=0;i<ExactVec.size();++i)
    PerturbVec[i]=RigidM*ExactVec[i];

  Quaternionf q;
  Point3f tr;
  Matrix44f res;
  ComputeRigidMatchMatrix(PerturbVec,ExactVec,res);

  res.print();
  RigidM.print();

  return 0;
}
