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
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/miq/MIQ.h>
#include <wrap/miq/quadrangulator.h>


using namespace vcg;

class CFace;
class CVertex;

struct MyUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType, Use<CFace>::AsFaceType >{};

class CVertex : public Vertex< MyUsedTypes,
        vertex::Coord3d,     vertex::Normal3d,
        vertex::BitFlags,    vertex::VFAdj,
        vertex::TexCoord2d,  vertex::Qualityd>{};

class CFace   : public Face<  MyUsedTypes, face::VertexRef,
        face::VFAdj,            face::FFAdj,      face::Normal3d,
        face::WedgeTexCoord2d,  face::BitFlags,   face::CurvatureDird,
        face::Qualityd,         face::Color4b,    face::Mark>{};

class CMesh   : public tri::TriMesh< std::vector<CVertex>, std::vector<CFace> >{};


class MyPolyFace;
class MyPolyVertex;
struct PolyUsedTypes: public UsedTypes<Use<MyPolyVertex>	::AsVertexType,
                                            Use<MyPolyFace>	::AsFaceType
                                            >{};

class MyPolyVertex:public Vertex<	PolyUsedTypes,
                                        vertex::Coord3f,
                                        vertex::Normal3f,
                                        vertex::BitFlags>{} ;

class MyPolyFace:public Face<   PolyUsedTypes,
    face::PolyInfo, // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
    face::PFVAdj,	 // Pointer to the vertices (just like FVAdj )
    face::BitFlags, // bit flags
    face::Normal3f // normal
> {};

class MyPolyMesh: public
    tri::TriMesh<     std::vector<MyPolyVertex>,
                           std::vector<MyPolyFace > 	  >{};


using namespace std;

int main(int argc, const char * argv[])
{
    const char* configfile;

    if (argc != 2)
    {
        cout << "Not enough parameters;\nUsage:\n\n./quadrangulator configfile" << endl;
        exit(EXIT_FAILURE);
    }

    configfile = argv[1];
    printf("configuration file %s",configfile);
    fflush(stdout);

    //  ********* Read parameters from config File *********
    FILE *fp= fopen(configfile,"r");
    if (fp==NULL)
    {
        printf("Cannot open config file %s\n",configfile);
        return -1;
    }

    char  meshFilename[200];      fscanf(fp,"mesh=%s\n", meshFilename);
    char  fieldFilename[200];     fscanf(fp,"field=%s\n", fieldFilename);
    int   scalegradient;          fscanf(fp,"scalegradient=%d\n", &scalegradient);
    float GradientSize;           fscanf(fp,"gradient=%f\n", &GradientSize);
    int   DirectRounding;         fscanf(fp,"directround=%d\n", &DirectRounding);
    char  outFilename[200];       fscanf(fp,"out=%s\n", outFilename);

    fclose(fp);

    // other parameters not in the config file...
    float Stiffness=4;
    int iterationNum=10;
    int localIterationNum=5;

    // ********* Start the processing *********
    // 1) Load, Prepare and Validate the mesh and field

    CMesh trimesh;
    if (tri::io::Importer<CMesh>::Open(trimesh, meshFilename)!=0)
    {
        printf("error loading mesh file %s\n",meshFilename);
        exit(0);
    }

    tri::UpdateBounding<CMesh>::Box(trimesh);
    tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(trimesh);
    tri::UpdateTopology<CMesh>::FaceFace(trimesh);
    tri::UpdateTopology<CMesh>::VertexFace(trimesh);

    if (std::string(fieldFilename).find(".ffield")==-1)
    {
        printf("error loading field file %s\n",fieldFilename);
        exit(0);
    }
    bool field_loaded=tri::io::ImporterFIELD<CMesh>::LoadFFIELD(trimesh,fieldFilename);
    if (!field_loaded) return false;

    if (scalegradient!=0) // If the gradient has not been specified try an educated guess
      GradientSize*=1.0/trimesh.bbox.Diag();

    if (! MIQ_parametrization<CMesh>::IsValid(trimesh) )
    {
        printf("Mesh not valid for quadrangulation: \n"
               "It must be a 2-manifold single connected component \n");
        exit(0);
    }

    // 2) Do the actual Quadrangulation and save the result

    MyPolyMesh polymesh;
    MIQ_parametrization<CMesh>::InitSeamsSing(trimesh,true,true,true);
    MIQ_parametrization<CMesh>::Parametrize(trimesh,MIQ_parametrization<CMesh>::ITERATIVE,Stiffness,GradientSize,(bool)DirectRounding,iterationNum,localIterationNum,true);

    Quadrangulator<CMesh,MyPolyMesh> Quad;
    Quad.TestIsProper(trimesh);
    Quad.Quadrangulate(trimesh,polymesh);

    tri::io::Exporter<MyPolyMesh>::Save(polymesh,outFilename);

}

