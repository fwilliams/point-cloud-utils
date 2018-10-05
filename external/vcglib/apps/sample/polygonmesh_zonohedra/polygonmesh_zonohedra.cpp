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

#include <stdio.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/zonohedron.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <wrap/io_trimesh/export_off.h>



class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes: public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType,vcg::Use<MyEdge>::AsEdgeType,vcg::Use<MyFace>::AsFaceType>{};

class MyVertex : public vcg::Vertex< MyUsedTypes,vcg::vertex::Coord3f,vcg::vertex::BitFlags >{};

class MyEdge : public vcg::Edge< MyUsedTypes > {};

class MyFace : public vcg::Face< MyUsedTypes,
    vcg::face::FFAdj,
    vcg::face::VertexRef,
  vcg::face::Normal3f,
    vcg::face::BitFlags > {};

// the main mesh class
class MyMesh : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};



// example 1: build a cube as a Zonohedron
void example1(){

    vcg::tri::Zonohedron<float> z;
    z.addVector( 0,0,1 );
    z.addVector( 0,1,0 );
    z.addVector( 1,0,0 );

    MyMesh m;
    z.createMesh(m); // this will be a cube

    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m); // needed by exporter

    int savemask = vcg::tri::io::Mask::IOM_BITPOLYGONAL;
    vcg::tri::io::ExporterOFF<MyMesh>::Save(m,"cube.off",savemask);
}


// example2: reads input file, builds zonohedra as described there
void example2(){

    FILE* f = fopen("input.txt","rt");
    if (!f) return;

    while (1) {

        // read mesh name
        char meshFilename[1024], fullMeshFilename[1024];
        if (fscanf(f,"%s",meshFilename)!=1) break;
        sprintf(fullMeshFilename,"%s.off",meshFilename);

        // build input vector
        vcg::tri::Zonohedron<float> z;
        while (1) {
            float a,b,c;
            if (fscanf(f,"%f %f %f",&a, &b, &c)!=3) break;
            z.addVector(a,b,c);
        }

        printf("Building %s from %d vectors...\n",fullMeshFilename, z.vectors().size() );

        MyMesh m;
        z.createMesh(m);

        vcg::tri::UpdateTopology<MyMesh>::FaceFace(m); // needed by exporter

        // normally, faces with more than 4sides are split into parallelograms
        // this merges them (optional, try removing it!)
        vcg::tri::PolygonSupport<MyMesh,MyMesh>::MergeFlatFaces(m);

        int savemask = vcg::tri::io::Mask::IOM_BITPOLYGONAL;
        vcg::tri::io::ExporterOFF<MyMesh>::Save(m,fullMeshFilename,savemask);
    }

}

int main(int argc, char *argv[]){
    example1();
    example2();
    return 0;
}
