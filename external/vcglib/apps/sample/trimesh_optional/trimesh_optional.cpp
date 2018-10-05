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

/*! \file trimesh_optional.cpp
\ingroup code_sample

\brief the minimal example of using the \ref optional_component "Optional Components".

This file shows how to dynamically allocate component that you do not need most of the time.

*/


#include<vcg/complex/complex.h>

#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <wrap/io_trimesh/export_off.h>

class MyFace;
class MyVertex;

struct MyUsedTypes:    public	vcg::UsedTypes<
    vcg::Use<MyVertex>::AsVertexType,
    vcg::Use<MyFace  >::AsFaceType>{};

class MyVertex     : public vcg::Vertex<	MyUsedTypes,
    vcg::vertex::Coord3f,  vcg::vertex::Qualityf,
    vcg::vertex::Color4b,  vcg::vertex::BitFlags,
    vcg::vertex::Normal3f, vcg::vertex::VFAdj >{};

class MyFace       : public vcg::Face< MyUsedTypes,
    vcg::face::FFAdj,   vcg::face::VFAdj,
    vcg::face::Color4b, vcg::face::VertexRef,
    vcg::face::BitFlags,   vcg::face::Normal3f > {};
class MyMesh       : public vcg::tri::TriMesh<     std::vector<MyVertex   >,           std::vector<MyFace   > > {};


class MyVertexOcf;
class MyFaceOcf;

struct MyUsedTypesOcf: public	vcg::UsedTypes<
    vcg::Use<MyVertexOcf>::AsVertexType,
    vcg::Use<MyFaceOcf>::AsFaceType>{};

class MyVertexOcf  : public vcg::Vertex< MyUsedTypesOcf,
    vcg::vertex::InfoOcf,       //   <--- Note the use of the 'special' InfoOcf component
    vcg::vertex::Coord3f,  vcg::vertex::QualityfOcf,
    vcg::vertex::Color4b,  vcg::vertex::BitFlags,
    vcg::vertex::Normal3f, vcg::vertex::VFAdjOcf >{};

class MyFaceOcf    : public vcg::Face< MyUsedTypesOcf,
    vcg::face::InfoOcf,         //   <--- Note the use of the 'special' InfoOcf component
    vcg::face::FFAdjOcf,   vcg::face::VFAdjOcf,
    vcg::face::Color4bOcf, vcg::face::VertexRef,
    vcg::face::BitFlags,   vcg::face::Normal3fOcf > {};

// the mesh class must make use of the 'vector_ocf' containers instead of the classical std::vector
class MyMeshOcf : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<MyVertexOcf>, vcg::face::vector_ocf<MyFaceOcf> > {};


using namespace vcg;
using namespace std;

int main(int , char **)
{
  MyMesh cm;
  MyMeshOcf cmof;

  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmof);

  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.VN(),cm.FN());

  assert(tri::HasFFAdjacency(cmof) == false);
  cmof.face.EnableFFAdjacency();
  assert(tri::HasFFAdjacency(cmof) == true);

  assert(tri::HasVFAdjacency(cmof) == false);
  cmof.vert.EnableVFAdjacency();
  cmof.face.EnableVFAdjacency();
  assert(tri::HasVFAdjacency(cmof) == true);

  tri::UpdateTopology<MyMesh   >::FaceFace(cm);
  tri::UpdateTopology<MyMeshOcf>::FaceFace(cmof);

  tri::UpdateFlags<MyMesh   >::FaceBorderFromFF(cm);
  tri::UpdateFlags<MyMeshOcf>::FaceBorderFromFF(cmof);

  tri::UpdateNormal<MyMesh   >::PerVertexPerFace(cm);
  cmof.face.EnableNormal();  // remove this line and you will throw an exception for a missing 'normal' component
  tri::UpdateNormal<MyMeshOcf>::PerVertexPerFace(cmof);

  cmof.vert.EnableQuality();

  tri::UpdateColor<MyMeshOcf>::PerVertexConstant(cmof,Color4b::LightGray);
  printf("cmof IsColorEnabled() %s\n",cmof.face.back().IsColorEnabled()?"Yes":"No");
  cmof.face.EnableColor();
  tri::UpdateColor<MyMeshOcf>::PerFaceConstant(cmof,Color4b::LightGray);
  cmof.vert[0].C()=Color4b::Red;

  printf("cmof IsColorEnabled() %s\n",cmof.face.back().IsColorEnabled()?"Yes":"No");
  printf("cm IsColorEnabled() %s\n",cm.face.back().IsColorEnabled()?"Yes":"No");
  printf("cm IsColorEnabled() %s\n",cm.face.back().HasColor()?"Yes":"No");

  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0,t2=0;
  while(float(t1-t0)/CLOCKS_PER_SEC < 0.1f)
  {
    t0=clock();
    tri::RefineOddEven<MyMesh> (cm, tri::OddPointLoop<MyMesh>(cm), tri::EvenPointLoop<MyMesh>(), 0);
    t1=clock();
    tri::RefineOddEven<MyMeshOcf> (cmof, tri::OddPointLoop<MyMeshOcf>(cmof), tri::EvenPointLoop<MyMeshOcf>(), 0);
    t2=clock();
  }
  printf("Last Iteration: Refined a tetra up to a mesh of %i faces in:\n"
         "Standard Component %5.2f sec\n"
         "OCF      Component %5.2f sec\n",cm.FN(),float(t1-t0)/CLOCKS_PER_SEC,float(t2-t1)/CLOCKS_PER_SEC);
  tri::io::ExporterOFF<MyMeshOcf>::Save(cmof,"test.off",tri::io::Mask::IOM_VERTCOLOR);

  return 0;
}
