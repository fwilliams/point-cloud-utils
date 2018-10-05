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

#include <sstream>

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/clean.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_off.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/algorithms/local_optimization/quad_diag_collapse.h>
#include <vcg/complex/algorithms/update/fitmaps.h>

using namespace vcg;
using namespace std;

// forward declarations
class CFace;
class CVertex;
class CHEdge;
class CEdge;
class MyPolyVertex;

struct CUsedTypes: public vcg::UsedTypes< vcg::Use<CVertex>::AsVertexType, vcg::Use<CFace>::AsFaceType >{};

// Mesh of triangles
class CVertex : public Vertex<
        CUsedTypes,
    vertex::BitFlags,
    vertex::Coord3f,
    vertex::Normal3f,
        vertex::VFAdj,
        vertex::Mark,
        vcg::vertex::Curvaturef,
        vcg::vertex::CurvatureDirf,
        vertex::Color4b,
        vertex::Qualityf
        >{};

class CFace   : public Face<
        CUsedTypes,
        face::VertexRef,
        face::Normal3f,
        face::BitFlags,
        face::FFAdj,
        face::VFAdj,
        face::Mark,
        face::EdgePlane
        > {};

class CMesh : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};




// Poly mesh
class MyPolyFace;
class MyPolyVertex;

struct PolyUsedTypes: public vcg::UsedTypes<
    vcg::Use<MyPolyVertex>  ::AsVertexType,
    vcg::Use<CEdge>         ::AsEdgeType,
    vcg::Use<CHEdge>        ::AsHEdgeType,
    vcg::Use<MyPolyFace>    ::AsFaceType
    >{};

class MyPolyVertex:public Vertex<
    PolyUsedTypes,
    vertex::Coord3f,
    vertex::Normal3f,
    vertex::Mark,
    vertex::BitFlags,
    vertex::VHAdj,
    vertex::VFAdj
    >{};

class CEdge : public Edge<PolyUsedTypes>{};

class CHEdge : public HEdge<
    PolyUsedTypes,
    hedge::BitFlags,
    hedge::HFAdj,
    hedge::HOppAdj,
    hedge::HNextAdj,
    hedge::HVAdj,
    hedge::HPrevAdj,
    hedge::Mark
    >{};

class MyPolyFace:public Face<
    PolyUsedTypes,
    face::PolyInfo,
    face::PFVAdj,
    face::PFFAdj,
    face::PFHAdj,
    face::BitFlags,
    face::Normal3f,
    face::Mark
    > {};

class MyPolyMesh: public tri::TriMesh<
    std::vector<MyPolyVertex>,
    std::vector<MyPolyFace>,
    std::vector<CHEdge>,
    std::vector<CEdge>
    >{};



/*!
  * \brief Collapse operation for adaptive simplification using fitmaps
  *
  */
class MyCollapseAdaptive: public vcg::tri::QuadDiagonalCollapse< MyPolyMesh, MyCollapseAdaptive, CMesh , vcg::tri::VertReg<MyPolyMesh> ,vcg::tri::FitmapsCollapse<MyPolyMesh, CMesh> , vcg::tri::FitmapsCollapse<MyPolyMesh, CMesh> >
{
public:

    typedef vcg::tri::QuadDiagonalCollapse< MyPolyMesh, MyCollapseAdaptive, CMesh , vcg::tri::VertReg<MyPolyMesh>, vcg::tri::FitmapsCollapse<MyPolyMesh, CMesh> , vcg::tri::FitmapsCollapse<MyPolyMesh, CMesh> > constructor;

    MyCollapseAdaptive(HEdgePointer he, int mark):constructor(he,mark){}
};

/*!
  * \brief Collapse for uniform simplification
  *
  */
class MyCollapse: public vcg::tri::QuadDiagonalCollapseBase< MyPolyMesh, MyCollapse, CMesh , vcg::tri::VertReg<MyPolyMesh> >
{
public:

    typedef vcg::tri::QuadDiagonalCollapseBase< MyPolyMesh, MyCollapse, CMesh , vcg::tri::VertReg<MyPolyMesh> > constructor;

    MyCollapse(HEdgePointer he, int mark):constructor(he,mark){}
};


typedef CMesh::FaceType TriFaceType;
typedef vcg::GridStaticPtr<CMesh::FaceType, TriFaceType::ScalarType> GRID;

typedef CMesh::PerVertexAttributeHandle<float> Fitmap_attr;

/*! Initializes the grid for smoothing and fitmaps
  *
  * \param m Reference mesh
  *
  */
void initGrid(CMesh & m)
{

  GRID* grid = new GRID();

  vcg::tri::UpdateBounding<CMesh>::Box(m);
  vcg::tri::UpdateEdges<CMesh>::Set(m);

  grid->Set(m.face.begin(), m.face.end());

//  grid->ShowStats(stdout);
  MyCollapse::grid() = grid;
  MyCollapseAdaptive::grid() = grid;

}

/*! Initializes the heap of operations on a mesh
  *
  * \param m Mesh
  * \param loc
  * \param Adaptive Specifies if simplificaton will be adaptive
  *
  */
void init_heap(MyPolyMesh &m, vcg::LocalOptimization<MyPolyMesh> &loc, bool adaptive)
{
    if(adaptive)
        MyCollapseAdaptive::Init(m, loc.h);
    else
        MyCollapse::Init(m,loc.h);

    std::make_heap(loc.h.begin(),loc.h.end());

    if(!loc.h.empty())
        loc.currMetric=loc.h.front().pri;
}

/*! Read fitmaps values from a file and loads them into a mesh
  *
  * \param m Mesh
  * \param fn Name of the file to read
  *
  */
bool read_fitmaps(CMesh &m, const char *fn)
{
    ifstream fitmaps;
    fitmaps.open(fn);

    if(! fitmaps.is_open())
        return false;

    Fitmap_attr S_Fit = tri::Allocator<CMesh>::GetPerVertexAttribute<float>(m,"S-Fitmap");
    Fitmap_attr M_Fit = tri::Allocator<CMesh>::GetPerVertexAttribute<float>(m,"M-Fitmap");

    int index;
    float S_fit, M_fit;
    do
    {
        fitmaps >> index >> S_fit >> M_fit;
        S_Fit[m.vert[index]] = S_fit;
        M_Fit[m.vert[index]] = M_fit;

    }while(fitmaps.good());


    bool eof = fitmaps.eof();

    fitmaps.close();
    return eof;

}

/*! Writes fitmaps values into a file
  *
  * \param m Mesh
  * \param fn Name of the file to write
  *
  */
bool store_fitmaps(CMesh &m, const char *fn)
{
    ofstream fitmaps;
    fitmaps.open(fn);

    if(! fitmaps.is_open())
        return false;

    Fitmap_attr S_Fit = tri::Allocator<CMesh>::GetPerVertexAttribute<float>(m,"S-Fitmap");
    Fitmap_attr M_Fit = tri::Allocator<CMesh>::GetPerVertexAttribute<float>(m,"M-Fitmap");

    for(unsigned int i =0; i< m.vert.size(); i++)
    {
        if( !(m.vert[i].IsD()) )
        {
            fitmaps << i << " " << S_Fit[m.vert[i]] << " " << M_Fit[m.vert[i]] << endl;

            if(!fitmaps.good())
            {
                fitmaps.close();
                return false;
            }
        }
    }

    fitmaps.close();
    return true;
}

/*! Load fitmaps for a mesh, computing them or reading values from a file
  *
  * \param m Mesh
  * \param fn Name of the mesh file
  *
  */
void load_fitmaps(CMesh &m, char* fn)
{

    Fitmap_attr S_Fit = tri::Allocator<CMesh>::AddPerVertexAttribute<float>  (m, string("S-Fitmap"));
    Fitmap_attr M_Fit = tri::Allocator<CMesh>::AddPerVertexAttribute<float>  (m, string("M-Fitmap"));

    string filename(fn);

    int found = filename.find_last_of("/");

    string name = filename.substr(found+1);

    string suffix = ".fmp";

    if( !read_fitmaps( m, (name + suffix).c_str()) )
    {
        tri::Fitmaps<CMesh>::computeSFitmap(m);

        for(CMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        S_Fit[vi] = vi->Q();

        tri::Fitmaps<CMesh>::computeMFitmap(m, 5);

        for(CMesh::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        M_Fit[vi] = vi->Q();

        store_fitmaps(m, ( name + suffix).c_str());
    }

}

/*! Load a mesh, from a file
  *
  * \param m Mesh that will be filled with data from a file
  * \param fn Name of the mesh file
  * \param loadFitmaps Specifies if fitmaps have to be loaded
  *
  */
void loadMesh(CMesh & m, char* fn, bool loadFitmaps = false)
{

    int ret = vcg::tri::io::Importer<CMesh>::Open(m,fn);

    if(ret != 0)
    {
        cerr << "Error reading file " << fn << endl;
        exit(1);
    }

    tri::Clean<CMesh>::RemoveDegenerateFace(m);
    tri::Clean<CMesh>::RemoveDuplicateFace(m);
    tri::Clean<CMesh>::RemoveDuplicateVertex(m);
    tri::Clean<CMesh>::RemoveUnreferencedVertex(m);

    tri::UpdateTopology<CMesh>::FaceFace(m);

    tri::Clean<CMesh>::RemoveNonManifoldFace(m);
    tri::UpdateTopology<CMesh>::FaceFace(m);

    tri::Clean<CMesh>::RemoveNonManifoldVertex(m);
    tri::UpdateTopology<CMesh>::FaceFace(m);

    // update bounding box
    vcg::tri::UpdateBounding<CMesh>::Box (m);

    // update Normals
    vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizedPerFace(m);
    vcg::tri::UpdateNormals<CMesh>::PerFaceNormalized(m);

    if(loadFitmaps)
        load_fitmaps(m,fn);

}

int main(int argc, char *argv[]) {


    // HE mesh
    MyPolyMesh pm;

    // Tri meshes
    CMesh mesh,refMesh;

    char* meshfile = NULL;
    char* trimeshfile = NULL;
    char* outfile = "output.off";
    int faces;
    bool adaptive = false;

    if(argc < 2)
    {
        cerr << "Usage: " << argv[0] << " -meshfile filename [-trimeshfile filename] -faces num_faces [-adaptive] [-outfile filename]" << endl;
    }

    for(int i=1; i< argc; i++)
    {
        string arg = string(argv[i]);

        if ( arg == "-meshfile")
            meshfile = argv[++i];

        else if (arg == "-trimeshfile")
            trimeshfile = argv[++i];

        else if (arg == "-faces")
        {
            stringstream myStream(argv[++i], stringstream::in | stringstream::out);
            myStream >> faces;
        }

        else if (arg == "-outfile")
            outfile = argv[++i];

        else if (arg == "-adaptive")
            adaptive = true;
    }


    if( !meshfile)
    {
        cerr << "Missing mesh filename" << endl;
        exit(1);
    }

    if(faces < 0)
    {
        cerr << "Missing faces number" << endl;
        exit(1);
    }


    // Load the mesh to simplify
    loadMesh(mesh, meshfile);

    // Load the reference mesh
    if(trimeshfile)
        loadMesh(refMesh, trimeshfile, adaptive);
    else
        loadMesh(refMesh, meshfile, adaptive);

    initGrid(refMesh);

    MyCollapse::refMesh() = &refMesh;
    MyCollapseAdaptive::refMesh() = &refMesh;


    vcg::tri::PolygonSupport<CMesh,MyPolyMesh>::ImportFromTriMesh(pm,mesh);
    vcg::tri::UpdateHalfEdges<MyPolyMesh>::FromIndexed(pm);

    // After loading check mesh consistency
    assert(vcg::tri::UpdateHalfEdges<MyPolyMesh>::CheckConsistency(pm));

    HalfedgeQuadClean<MyPolyMesh>::remove_singlets(pm);
    HalfedgeQuadClean<MyPolyMesh>::remove_doublets(pm);

    vcg::LocalOptimization<MyPolyMesh> loc(pm);
    init_heap(pm, loc, adaptive);

    loc.HeapSimplexRatio = 9;
    loc.SetTargetSimplices(faces);

    // Perform simplification
    loc.DoOptimization();


    assert(vcg::tri::UpdateHalfEdges<MyPolyMesh>::CheckConsistency(pm));
    vcg::tri::UpdateIndexed<MyPolyMesh>::FromHalfEdges(pm );


    int ret = tri::io::ExporterOFF<MyPolyMesh>::Save(pm, outfile, tri::io::Mask::IOM_BITPOLYGONAL );
    if(ret != 0 )
    {
        cerr << "Error saving file" << endl;
        exit(1);
    }

    cout << "Simplification ended successfully!" << endl;

}
