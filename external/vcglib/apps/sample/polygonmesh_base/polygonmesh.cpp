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

/*include the algorithms for updating: */
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/create/platonic.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

/* include the support for polygon meshes (function to convert from/to trimesh)*/
//#include <vcg/complex/algorithms/polygon_support.h>

/* include the support for half edges */
#include <vcg/complex/algorithms/update/halfedge_indexed.h>


using namespace vcg;
using namespace std;

// forward declarations
class TFace;
class TVertex;

struct TUsedTypes: public vcg::UsedTypes< vcg::Use<TVertex>::AsVertexType, vcg::Use<TFace>::AsFaceType >{};



/* Definition of a mesh of triangles
*/
class TVertex : public Vertex< TUsedTypes,
    vertex::BitFlags,
    vertex::Coord3f,
    vertex::Normal3f,
    vertex::Mark >{};

class TFace   : public Face<   TUsedTypes,
    face::VertexRef,	// three pointers to vertices
    face::Normal3f,		// normal
    face::BitFlags,		// flags
    face::FFAdj			// three pointers to adjacent faces
> {};

/* the mesh is a container of vertices and a container of faces */
class TMesh   : public vcg::tri::TriMesh< vector<TVertex>, vector<TFace> > {};


/* Definition of a mesh of polygons that also supports half-edges
*/
class PFace;
class PVertex;
class PHEdge;
class PEdge;

struct PUsedTypes: public vcg::UsedTypes<vcg::Use<PVertex>  ::AsVertexType,
                                            vcg::Use<PEdge>	::AsEdgeType,
                                            vcg::Use<PHEdge>::AsHEdgeType,
                                            vcg::Use<PFace>	::AsFaceType
                                            >{};

//class DummyEdge: public vcg::Edge<PolyUsedTypes>{};
class PVertex:public vcg::Vertex<	PUsedTypes,
                                        vcg::vertex::Coord3f,
                                        vcg::vertex::Normal3f,
                                        vcg::vertex::Mark,
                                        vcg::vertex::BitFlags,
                                        vcg::vertex::VHAdj>{} ;

class PEdge : public Edge<PUsedTypes>{};
class PHEdge : public HEdge< PUsedTypes, hedge::BitFlags,
    //hedge::HFAdj,		// pointer to the face
    //hedge::HOppAdj,	// pointer to the opposite edge
    //hedge::HVAdj,		// pointer to the vertex
    //hedge::HNextAdj,	// pointer to the next halfedge
    hedge::HEdgeData		// the previous 4 components (just more handy, you can comment this and uncomment the previous four lines)
    //,hedge::HPrevAdj	// pointer to the previous halfedge
>{};

class PFace:public vcg::Face<
     PUsedTypes
    ,vcg::face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
                         // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
                         // will be allocated dynamically")
    ,vcg::face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
    ,vcg::face::PFVAdj
    ,vcg::face::PFFAdj	 // Pointer to edge-adjacent face (just like FFAdj )
    ,vcg::face::PFHAdj	 // Pointer its half -edges  ( you may need this if you use half edges)
    ,vcg::face::BitFlags // bit flags
    ,vcg::face::Normal3f // normal
> {};

class PMesh: public
    vcg::tri::TriMesh<
    std::vector<PVertex>,	// the vector of vertices
    std::vector<PFace >, 						// the vector of faces
    std::vector<PHEdge>		,						// the vector of edges
    std::vector<PEdge> 								// the vector of edges
    >{};

PMesh pm;
TMesh tm0;

int	main(int argc, char *argv[]) {

    int loadmask;

if(true){
    /*
    first way:
    1) read a polygon mesh that will be automatically converted in a triangle mesh tagging
      the internal edges (i.e. the edges that have been added for triangulating the polygons)
    2) make some cleaning
    3) import the tagged triangle mesh in a polygon mesh
    */
//	vcg::tri::io::ImporterOBJ<CMesh>::Open(mesh,argv[1],loadmask);
//    vcg::tri::io::ImporterOFF<TMesh>::Open(tm0,argv[1],loadmask);
    vcg::tri::Hexahedron(tm0);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(tm0);
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(tm0);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(tm0);
    vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(tm0);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(tm0);
    assert(vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(tm0)==0);
    assert(vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(tm0)==0);

    // create a polygon meshe from a trimesh with tagged faces
    vcg::tri::PolygonSupport<TMesh,PMesh>::ImportFromTriMesh(pm,tm0);
}
else
{
    /* second way:
    Load into a polygon mesh straight away.
    */
    vcg::tri::io::ImporterOBJ<PMesh>::Open(pm,argv[1],loadmask);
    vcg::tri::UpdateTopology<PMesh>::FaceFace(pm);
    vcg::tri::Clean<PMesh>::RemoveNonManifoldFace(pm);
    vcg::tri::UpdateTopology<PMesh>::FaceFace(pm);
    assert(vcg::tri::Clean<PMesh>::CountNonManifoldEdgeFF(pm));
}


    // compute the half edges because I'm a half-edge programmer
    vcg::tri::UpdateHalfEdges<PMesh>::FromIndexed(pm);

    // .... my half edge based code ......

    // check for consistency
    assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(pm));

    int size =  pm.face.size();

    // add a face to each face with more than 3 vertices ( just one pass)

    for(int i = 0; i < size; ++i)
        if(!(pm.face[i].IsD()))
        if(pm.face[i].VN()>3){
            PMesh::HEdgePointer ef =  pm.face[i].FHp();
            PMesh::HEdgePointer ef1 = ef -> HNp();
            ef1 = ef1->HNp();
            vcg::tri::UpdateHalfEdges<PMesh>::AddHEdge(pm, ef, ef1 );
        }
    assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(pm));
    size =  pm.face.size();

    // remove an edge for each face

    for(int i = 0; i < size; ++i)
        if(!(pm.face[i].IsD() ))
        {
            PMesh::HEdgePointer ef =  pm.face[i].FHp();
            if( ef->HOp()->HFp() !=NULL){
                vcg::tri::UpdateHalfEdges<PMesh>::RemoveHEdge(pm,ef);
            }
        }

    // check for consistency
    assert(vcg::tri::UpdateHalfEdges<PMesh>::CheckConsistency(pm));

    // recompute indexed data structure from the half edge data structure
//    vcg::tri::UpdateIndexed<PMesh>::FromHalfEdges(pm );

    // create a triangle mesh from a polygon mesh
    TMesh tm1;
    vcg::tri::PolygonSupport<TMesh,PMesh>::ImportFromPolyMesh(tm1,pm);

    vcg::tri::io::PlyInfo pi;
    vcg::tri::io::ExporterPLY<TMesh>::Save(tm1,"converted_tri.ply",false,pi);
    vcg::tri::io::ExporterPLY<PMesh>::Save(pm,"converted_poly.ply",false,pi);
}



