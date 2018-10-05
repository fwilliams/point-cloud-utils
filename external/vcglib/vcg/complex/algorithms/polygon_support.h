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

#ifndef __VCGLIB_POLYGON_SUPPORT
#define __VCGLIB_POLYGON_SUPPORT

#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/space/planar_polygon_tessellation.h>

namespace vcg {
namespace tri {
    /// \ingroup trimesh

    /// \headerfile polygon_support.h vcg/complex/algorithms/polygon_support.h

    /// \brief This class is used convert between polygonal meshes and triangular meshes

    /**
    This class contains two members that allow to build a triangular mesh from a polygonal mesh
    and viceversa. In a trimesh, the generic polygons with n sides are codified represented by
    tagging the internal edge of the face as 'faux' with the SetF.
    */

    template <class TriMeshType,class PolyMeshType >
    class PolygonSupport{
      typedef typename TriMeshType::FaceIterator  TriFaceIterator;
      typedef typename PolyMeshType::FaceIterator PolyFaceIterator;
      typedef typename TriMeshType::VertexIterator  TriVertexIterator;
      typedef typename PolyMeshType::VertexIterator PolyVertexIterator;
      typedef typename TriMeshType::CoordType::ScalarType Scalar;

    public:
    /**
    Given a tri mesh (with per-face normals and FF connectivity),
    merges flat faces into larger polygons.
    The merging is done only by setting the faux bit
    **/
    static void MergeFlatFaces(TriMeshType & tm, double tolerance = 0.1E-4)
    {
        typedef typename TriMeshType::FaceType FaceType;
        Scalar minDist = 1 - Scalar(tolerance);
        for (TriFaceIterator fi=tm.face.begin(); fi!=tm.face.end(); fi++) {
            FaceType *fa = &*fi;
            for (int w=0; w<3; w++) {
                FaceType *fb = fa->FFp(w);
                if ( (fb>fa) && (fa->N()*fb->N() > minDist) ) {
                    fa->SetF( w );
                    fb->SetF( fa->FFi(w) ); // reciprocate
                }
            }
        }
    }

    /**
    Import a  trianglemesh from a polygon mesh
    **/
    static void ImportFromPolyMesh(TriMeshType & tm,  PolyMeshType & pm)
    {
      tri::RequirePolygonalMesh(pm);
      std::vector<typename PolyMeshType::CoordType> points;

      // the vertices are the same, simply import them
      PolyVertexIterator vi;
      TriVertexIterator tvi = Allocator<TriMeshType>::AddVertices(tm,pm.vert.size());
      int cnt = 0;
      for(tvi = tm.vert.begin(),vi = pm.vert.begin(); tvi != tm.vert.end(); ++tvi,++vi,++cnt)
        if(!(*vi).IsD()) (*tvi).ImportData(*vi); else tri::Allocator<TriMeshType>::DeleteVertex(tm,(*tvi));

      for(PolyFaceIterator fi = pm.face.begin(); fi != pm.face.end(); ++fi)
      {
        if(!((*fi).IsD())){
          points.clear();
          for(int i  = 0; i < (*fi).VN(); ++i) {
            typename	PolyMeshType::VertexType * v = (*fi).V(i);
            points.push_back(v->P());
          }
          std::vector<int> faces;
          TessellatePlanarPolygon3(points,faces);
          for(size_t i = 0; i<faces.size();i+=3){
            TriFaceIterator tfi = Allocator<TriMeshType>::AddFace(tm,
                  tri::Index(pm,(*fi).V( faces[i+0] )),
                  tri::Index(pm,(*fi).V( faces[i+1] )),
                  tri::Index(pm,(*fi).V( faces[i+2] )) );

            tfi->ImportData(*fi);
            // set the F flags
            if( (faces[i  ]+1)%points.size() != size_t(faces[i+1])) (*tfi).SetF(0);
            if( (faces[i+1]+1)%points.size() != size_t(faces[i+2])) (*tfi).SetF(1);
            if( (faces[i+2]+1)%points.size() != size_t(faces[i  ])) (*tfi).SetF(2);
          }
        }
      }
    }


    /**
    \brief Import a polygon mesh from a triangle mesh

    It assumes that the mesh has the faux edges bit set for a polygonal mesh and that have the FFAdjacency already computed.
    **/
    static void ImportFromTriMesh( PolyMeshType & pm,  TriMeshType & tm)
    {
      tri::RequirePolygonalMesh(pm);
      tri::RequireTriangularMesh(tm);

      tri::RequireCompactness(tm);
      tri::RequireFFAdjacency(tm);
      tri::UpdateFlags<TriMeshType>::FaceClearV(tm);
      // the vertices are the same, simply import them
        int cnt = 0;
        typename TriMeshType ::ConstVertexIterator tvi;
        typename PolyMeshType::VertexIterator vi  = tri::Allocator<PolyMeshType>::AddVertices(pm,tm.vert.size());
        for(tvi = tm.vert.begin(); tvi != tm.vert.end(); ++tvi,++vi,++cnt)
            (*vi).ImportData(*tvi);

        // convert the faces
        typename TriMeshType::FaceIterator tfi;
        face::JumpingPos<typename TriMeshType::FaceType> p;

        for( tfi = tm.face.begin(); tfi != tm.face.end(); ++tfi) if(!(*tfi).IsV())
        {
            std::vector<typename TriMeshType::VertexPointer> vs;// vertices of the polygon
            ExtractPolygon(&*tfi,vs);
            std::reverse(vs.begin(),vs.end());
            //now vs  contains all the vertices of the polygon (still in the trimesh)
            if (vs.size()==0)continue;
            typename PolyMeshType::FaceIterator pfi =  tri::Allocator<PolyMeshType>::AddFaces(pm,1);
            (*pfi).Alloc(vs.size());
            for(size_t i  = 0 ; i < vs.size(); ++i)
                (*pfi).V(i) = ( typename PolyMeshType::VertexType*)  & pm.vert[vs[i]-&(*tm.vert.begin())];
            if(tri::HasPerFaceColor(tm) && tri::HasPerFaceColor(pm)) pfi->C()=tfi->C();
            if(tri::HasPerFaceQuality(tm) && tri::HasPerFaceQuality(pm)) pfi->Q()=tfi->Q();
        }
    }
    /// \brief Collect tris and verts of a polygonal face marked by FauxEdges
    ///
    /// Given a face pointer, it builds a vector with all the face and vertex pointers
    /// around the polygonal face determined by the current FAUX-EDGE markings.
    /// It assumes that the mesh is 2Manifold and has FF adjacency already computed
    /// per face visited flag cleared.
    /// NOTE: All the faces touched are marked as visited and it assumes that you
    /// do not call this function on a visited face.
    static void ExtractPolygon(typename TriMeshType::FacePointer tfp,
                               std::vector<typename TriMeshType::VertexPointer> &vs,
                               std::vector<typename TriMeshType::FacePointer> &fs)
    {
        vs.clear();
        fs.clear();
        // find a non faux edge
        int se = -1;
        for(int i=0; i<3; i++) if (!( tfp->IsF(i))) { se = i; break;}

        // all faux edges return an empty vertex vector!
        if(se==-1) return;
        if(tfp->IsV()) return;

        // initialize a pos on the first non faux edge
        face::Pos<typename TriMeshType::FaceType> start(tfp,se,tfp->V(se));
        face::Pos<typename TriMeshType::FaceType> p(start);

        fs.push_back(p.F());
        p.F()->SetV();

        do
        {
            assert(!p.F()->IsF(p.E()));
            vs.push_back(p.V());
            p.FlipE();
            while( p.F()->IsF(p.E()) )
            {
                p.FlipF();
                if(!p.F()->IsV()) {
                  fs.push_back(p.F());
                  p.F()->SetV();
                }
                p.FlipE();
            }
            p.FlipV();
        } while(p!=start);
        //assert(vs.size() == fs.size()+2);
    }
    static void ExtractPolygon(typename TriMeshType::FacePointer tfp, std::vector<typename TriMeshType::VertexPointer> &vs)
    {
      std::vector<typename TriMeshType::FacePointer> fs;
      ExtractPolygon(tfp,vs,fs);
    }
}; // end of struct
}} // end namespace tri

#endif // __VCGLIB_TRI_CLIP
