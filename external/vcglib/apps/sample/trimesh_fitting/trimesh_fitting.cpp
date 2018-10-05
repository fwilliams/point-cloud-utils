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
#include<vcg/space/fitting3.h>


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

  // As a simple test just run over all the faces of a mesh
  // get a few random points over it, perturb them and fit a plane on them.
  vcg::Plane3f ple,plf,plw;
  int cnt=0;
  float scaleFac = m.bbox.Diag()/10.0f;
  printf("ScaleFac %f\n\n",scaleFac);
  vcg::math::MarsenneTwisterRNG rnd;
  for(int i=0;i<m.FN();++i)
  {
    std::vector<vcg::Point3f> ExactVec;
    std::vector<vcg::Point3f> PerturbVec;
    std::vector<float> WeightVec;
    vcg::Plane3f pl;
    pl.Init(vcg::Barycenter(m.face[i]),m.face[i].N());
    for(int j=0;j<200;++j)
    {
      vcg::Point3f p = vcg::tri::SurfaceSampling<MyMesh>::RandomPointInTriangle(m.face[i]);
      ExactVec.push_back(p);
      vcg::Point3f off = vcg::math::GeneratePointInUnitBallUniform<float>(rnd);
      p+=off*scaleFac;
      float w =  std::max(0.0f, 1.0f-fabs(vcg::SignedDistancePlanePoint(pl,p))/scaleFac);
      PerturbVec.push_back(p);
      WeightVec.push_back(w*w); // as weight we use the square of  (1-distance)
    }

    vcg::FitPlaneToPointSet(ExactVec,ple);
    float err=EvalPlane(ple,ExactVec);

    vcg::FitPlaneToPointSet(PerturbVec,plf);
    float err0=EvalPlane(plf,ExactVec);

    vcg::WeightedFitPlaneToPointSet(PerturbVec,WeightVec,plw);
    float err1=EvalPlane(plw,ExactVec);
    printf("Exact %5.3f Fit to Perturbed %5.3f Weighted fit to perturbed %5.3f\n",err,err0,err1);
    if(err0>err1) cnt++;
  }

  printf("\nWeighted Fitting was better %i on %i\n",cnt,m.FN());
  return 0;
}
