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
/*! \file trimesh_color.cpp
\ingroup code_sample

\brief a brief overview of the various color oriented functionalities

*/
#include<vcg/complex/complex.h>

#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clustering.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,    vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};


int  main()
{
 MyMesh m;

 vcg::tri::io::ImporterPLY<MyMesh>::Open(m,"../../meshes/torus_irregular.ply");

 vcg::tri::UpdateColor<MyMesh>::PerVertexConstant(m, vcg::Color4b::LightGray);
 vcg::tri::UpdateColor<MyMesh>::PerFaceConstant(m, vcg::Color4b::LightGray);

 vcg::tri::UpdateColor<MyMesh>::PerVertexPerlinNoise(m,vcg::Point3f(0.5,0.75,1.0));
 vcg::tri::UpdateColor<MyMesh>::PerFaceFromVertex(m);


 vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"out.ply",vcg::tri::io::Mask::IOM_FACECOLOR+vcg::tri::io::Mask::IOM_VERTCOLOR);
  return 0;
}
