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
/*! \file trimesh_inertia.cpp
\ingroup code_sample

\brief An example of computing the inertia properties of meshes

Two meshes are created a rectangular box and a torus and their mass properties are computed and shown.
The result should match the closed formula for these objects (with a reasonable approximation)

*/

#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import_off.h>

#include<vcg/complex/algorithms/inertia.h>
#include<vcg/complex/algorithms/create/platonic.h>

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                            vcg::Use<MyEdge>     ::AsEdgeType,
                                            vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::FFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

int main( int argc, char **argv )
{
  MyMesh boxMesh,torusMesh;
  vcg::Matrix33f IT;
  vcg::Point3f ITv;

  vcg::tri::Hexahedron(boxMesh);
  vcg::Matrix44f ScaleM,TransM;
  ScaleM.SetScale(1.0f, 2.0f, 5.0f);
  TransM.SetTranslate(2.0f,3.0f,4.0f);
  vcg::tri::UpdatePosition<MyMesh>::Matrix(boxMesh,ScaleM);
  vcg::tri::UpdatePosition<MyMesh>::Matrix(boxMesh,TransM);
  vcg::tri::Inertia<MyMesh> Ib(boxMesh);
  vcg::Point3f cc = Ib.CenterOfMass();
  Ib.InertiaTensorEigen(IT,ITv);

  printf("Box of size 2,4,10, centered in (2,3,4)\n");
  printf("Volume %f \n",Ib.Mass());
  printf("CenterOfMass %f %f %f\n",cc[0],cc[1],cc[2]);
  printf("InertiaTensor Values  %6.3f %6.3f %6.3f\n",ITv[0],ITv[1],ITv[2]);
  printf("InertiaTensor Matrix\n");

  printf(" %6.3f %6.3f %6.3f\n",IT[0][0],IT[0][1],IT[0][2]);
  printf(" %6.3f %6.3f %6.3f\n",IT[1][0],IT[1][1],IT[1][2]);
  printf(" %6.3f %6.3f %6.3f\n",IT[2][0],IT[2][1],IT[2][2]);

  // Now we have a box with sides (h,w,d) 2,4,10, centered in (2,3,4)
  // Volume is 80
  // inertia tensor should be:
  // I_h = 1/12 m *(w^2+d^2)  = 1/12 * 80 * (16+100) = 773.33
  // I_w = 1/12 m *(h^2+d^2)  = 1/12 * 80 * (4+100)  = 693.33
  // I_d = 1/12 m *(h^2+w^2)  = 1/12 * 80 * (4+16)   = 133.33


  vcg::tri::Torus(torusMesh,2,1,1024,512);
  vcg::tri::Inertia<MyMesh> It(torusMesh);
  cc = It.CenterOfMass();
  It.InertiaTensorEigen(IT,ITv);

  printf("\nTorus of radius 2,1\n");
  printf("Mass %f \n",It.Mass());
  printf("CenterOfMass %f %f %f\n",cc[0],cc[1],cc[2]);
  printf("InertiaTensor Values  %6.3f %6.3f %6.3f\n",ITv[0],ITv[1],ITv[2]);
  printf("InertiaTensor Matrix\n");

  printf(" %6.3f %6.3f %6.3f\n",IT[0][0],IT[0][1],IT[0][2]);
  printf(" %6.3f %6.3f %6.3f\n",IT[1][0],IT[1][1],IT[1][2]);
  printf(" %6.3f %6.3f %6.3f\n",IT[2][0],IT[2][1],IT[2][2]);

  /*
     Now we have a torus with c = 2, a = 1
     c = radius of the ring
     a = radius of the section

    Volume is:
    V= 2 PI^2 * a^2 * c = ~39.478

    Inertia tensor should be:

    | ( 5/8 a^2 + 1/2 c^2 ) M             0                         0     |
    |            0             ( 5/8 a^2 + 1/2 c^2 ) M              0     | =
    |            0                        0             (3/4 a^2 + c^2) M |

    | ( 5/8+2 ) M        0           0     |   | 103.630    0        0     |
  = |      0         ( 5/8+2 ) M     0     | = |    0    103.630     0     |
    |      0             0       (3/4+2) M |   |    0       0      187.52  |

  */

  return 0;
}
