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
#include<vcg/complex/append.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/clip.h>
#include<vcg/complex/algorithms/update/bounding.h>

// input output
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>

// std
#include <vector>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
                                                                                Use<MyFace>			::AsFaceType>{};

class MyVertex  : public Vertex <MyUsedTypes, vertex::Coord3f, vertex::BitFlags  >{};
class MyFace    : public Face   < MyUsedTypes, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> > {};



int main(int argc,char **argv )
{
 if(argc<2)
    {
        printf(	"Usage: trimesh_join [opt] filename.ply [filename.ply | *] \n"
                "where opt can be:\n"
                "   -b xmin ymin zmin xmax ymax zmax  : \n"
                "      Returns only mesh composed by faces inside specified bbox\n"
                "   -t Just scan all the input files computing the total bbox\n"
            );
        exit(0);
    }

  MyMesh ml,mr;
  Box3f ClipBB,TotBB;
  bool ClipFlag=false,MergeFlag=true;
    int i=1;
  // Parsing option loop
  while(argv[i][0]=='-')
  {
    switch(argv[i][1])
    {
      case 'b': {
              if(argc<i+7) {
                printf("Error in parsing bbox option");
                  exit(0);
              }
              ClipBB.min=Point3f::Construct(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
              ClipBB.max=Point3f::Construct(atof(argv[i+4]),atof(argv[i+5]),atof(argv[i+6]));
              i+=6;
              printf("Clipping incoming meshes with box:\n (%7.4f %7.4f %7.4f) - (%7.4f %7.4f %7.4f)\n",
                    ClipBB.min[0],ClipBB.min[1],ClipBB.min[2],
                    ClipBB.max[0],ClipBB.max[1],ClipBB.max[2]);
              ClipFlag=true;
          } break;
    case 't': MergeFlag=false; break;
    default : printf("Unknown option '%s'\n",argv[i]);
    }
  ++i;
  }
    while(i<argc)
        {
          if(vcg::tri::io::ImporterPLY<MyMesh>::Open(mr,argv[i])!=0)
          {
        printf("Error reading file  %s\n",argv[1]);
              exit(0);
          }
      printf("Input mesh %3i           vn:%9i fn:%9i\n",i, mr.VN(), mr.FN());
      if(ClipFlag)
      {
        tri::GenericVertexInterpolator<MyMesh> interp(mr);
        tri::TriMeshClipper<MyMesh>::Box(ClipBB,interp,mr);
        printf("              clipped to vn:%9i fn:%9i\n", mr.VN(), mr.FN());
      }
      tri::UpdateBounding<MyMesh>::Box(mr);
      TotBB.Add(mr.bbox);

      if(MergeFlag) tri::Append<MyMesh,MyMesh>::Mesh(ml,mr); // append mesh mr to ml
        ++i;
        }

  printf("Output mesh vn:%i fn:%i\n",ml.VN(),ml.FN());

  tri::io::ExporterPLY<MyMesh>::Save(ml,"joined.ply");
  int dv=tri::Clean<MyMesh>::RemoveDuplicateVertex(ml);
  printf("Removed %i duplicated vertices\n",dv);
  tri::io::ExporterPLY<MyMesh>::Save(ml,"joined_unif.ply");
  printf("Final BBox of mesh :\n (%7.4f %7.4f %7.4f) - (%7.4f %7.4f %7.4f)\n",
                    TotBB.min[0],TotBB.min[1],TotBB.min[2],
                    TotBB.max[0],TotBB.max[1],TotBB.max[2]);

}

