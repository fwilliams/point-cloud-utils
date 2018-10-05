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

#ifndef __VCGLIB_CLEAN
#define __VCGLIB_CLEAN

// VCG headers
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/simplex/edge/topology.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/update/selection.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/space/triangle3.h>

namespace vcg {
namespace tri{

template <class ConnectedMeshType>
class ConnectedComponentIterator
{
public:
  typedef ConnectedMeshType MeshType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::ConstFaceIterator   ConstFaceIterator;
  typedef typename MeshType::FaceContainer  FaceContainer;

public:
  void operator ++()
  {
    FacePointer fpt=sf.top();
    sf.pop();
    for(int j=0;j<3;++j)
      if( !face::IsBorder(*fpt,j) )
      {
        FacePointer l=fpt->FFp(j);
        if( !tri::IsMarked(*mp,l) )
        {
          tri::Mark(*mp,l);
          sf.push(l);
        }
      }
  }

  void start(MeshType &m, FacePointer p)
  {
    tri::RequirePerFaceMark(m);
    mp=&m;
    while(!sf.empty()) sf.pop();
    UnMarkAll(m);
    tri::Mark(m,p);
    sf.push(p);
  }

  bool completed() {
    return sf.empty();
  }

  FacePointer operator *()
  {
    return sf.top();
  }
private:
  std::stack<FacePointer> sf;
  MeshType *mp;
};


///
/** \addtogroup trimesh */
/*@{*/
/// Class of static functions to clean//restore meshs.
template <class CleanMeshType>
class Clean
{

public:
  typedef CleanMeshType MeshType;
  typedef typename MeshType::VertexType           VertexType;
  typedef typename MeshType::VertexPointer        VertexPointer;
  typedef typename MeshType::VertexIterator       VertexIterator;
  typedef typename MeshType::ConstVertexIterator  ConstVertexIterator;
  typedef typename MeshType::EdgeIterator         EdgeIterator;
  typedef typename MeshType::EdgePointer          EdgePointer;
  typedef typename MeshType::CoordType            CoordType;
  typedef typename MeshType::ScalarType           ScalarType;
  typedef typename MeshType::FaceType             FaceType;
  typedef typename MeshType::FacePointer          FacePointer;
  typedef typename MeshType::FaceIterator         FaceIterator;
  typedef typename MeshType::ConstFaceIterator    ConstFaceIterator;
  typedef typename MeshType::FaceContainer        FaceContainer;
  typedef typename vcg::Box3<ScalarType>  Box3Type;

  typedef GridStaticPtr<FaceType, ScalarType > TriMeshGrid;

  /* classe di confronto per l'algoritmo di eliminazione vertici duplicati*/
  class RemoveDuplicateVert_Compare{
  public:
    inline bool operator()(VertexPointer const &a, VertexPointer const &b)
    {
        return ((*a).cP() == (*b).cP()) ? (a<b): ((*a).cP() < (*b).cP());
    }
  };


  /** This function removes all duplicate vertices of the mesh by looking only at their spatial positions.
    *  Note that it does not update any topology relation that could be affected by this like the VT or TT relation.
    *  the reason this function is usually performed BEFORE building any topology information.
    */
  static int RemoveDuplicateVertex( MeshType & m, bool RemoveDegenerateFlag=true)    // V1.0
  {
    if(m.vert.size()==0 || m.vn==0) return 0;

    std::map<VertexPointer, VertexPointer> mp;
    size_t i,j;
    VertexIterator vi;
    int deleted=0;
    int k=0;
    size_t num_vert = m.vert.size();
    std::vector<VertexPointer> perm(num_vert);
    for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi, ++k)
      perm[k] = &(*vi);

    RemoveDuplicateVert_Compare c_obj;

    std::sort(perm.begin(),perm.end(),c_obj);

    j = 0;
    i = j;
    mp[perm[i]] = perm[j];
    ++i;
    for(;i!=num_vert;)
    {
      if( (! (*perm[i]).IsD()) &&
          (! (*perm[j]).IsD()) &&
          (*perm[i]).P() == (*perm[j]).cP() )
      {
        VertexPointer t = perm[i];
        mp[perm[i]] = perm[j];
        ++i;
        Allocator<MeshType>::DeleteVertex(m,*t);
        deleted++;
      }
      else
      {
        j = i;
        ++i;
      }
    }

    for(FaceIterator fi = m.face.begin(); fi!=m.face.end(); ++fi)
      if( !(*fi).IsD() )
        for(k = 0; k < (*fi).VN(); ++k)
          if( mp.find( (typename MeshType::VertexPointer)(*fi).V(k) ) != mp.end() )
          {
            (*fi).V(k) = &*mp[ (*fi).V(k) ];
          }


    for(EdgeIterator ei = m.edge.begin(); ei!=m.edge.end(); ++ei)
      if( !(*ei).IsD() )
        for(k = 0; k < 2; ++k)
          if( mp.find( (typename MeshType::VertexPointer)(*ei).V(k) ) != mp.end() )
          {
            (*ei).V(k) = &*mp[ (*ei).V(k) ];
          }
    if(RemoveDegenerateFlag) RemoveDegenerateFace(m);
    if(RemoveDegenerateFlag && m.en>0) {
      RemoveDegenerateEdge(m);
      RemoveDuplicateEdge(m);
    }
    return deleted;
  }

  class SortedPair
  {
  public:
    SortedPair() {}
    SortedPair(unsigned int v0, unsigned int v1, EdgePointer _fp)
    {
      v[0]=v0;v[1]=v1;
      fp=_fp;
      if(v[0]>v[1]) std::swap(v[0],v[1]);
    }
    bool operator < (const SortedPair &p) const
    {
      return (v[1]!=p.v[1])?(v[1]<p.v[1]):
        (v[0]<p.v[0]);				}

      bool operator == (const SortedPair &s) const
      {
      if( (v[0]==s.v[0]) && (v[1]==s.v[1]) ) return true;
      return false;
    }

    unsigned int v[2];
    EdgePointer fp;
  };
  class SortedTriple
  {
  public:
    SortedTriple() {}
    SortedTriple(unsigned int v0, unsigned int v1, unsigned int v2,FacePointer _fp)
    {
      v[0]=v0;v[1]=v1;v[2]=v2;
      fp=_fp;
      std::sort(v,v+3);
    }
    bool operator < (const SortedTriple &p) const
    {
      return (v[2]!=p.v[2])?(v[2]<p.v[2]):
        (v[1]!=p.v[1])?(v[1]<p.v[1]):
          (v[0]<p.v[0]);				}

      bool operator == (const SortedTriple &s) const
      {
      if( (v[0]==s.v[0]) && (v[1]==s.v[1]) && (v[2]==s.v[2]) ) return true;
      return false;
    }

    unsigned int v[3];
    FacePointer fp;
  };


  /** This function removes all duplicate faces of the mesh by looking only at their vertex reference.
      So it should be called after unification of vertices.
      Note that it does not update any topology relation that could be affected by this like the VT or TT relation.
      the reason this function is usually performed BEFORE building any topology information.
     */
  static int RemoveDuplicateFace( MeshType & m)    // V1.0
  {
    std::vector<SortedTriple> fvec;
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if(!(*fi).IsD())
      {
        fvec.push_back(SortedTriple(    tri::Index(m,(*fi).V(0)),
                                        tri::Index(m,(*fi).V(1)),
                                        tri::Index(m,(*fi).V(2)),
                                        &*fi));
      }
    std::sort(fvec.begin(),fvec.end());
    int total=0;
    for(int i=0;i<int(fvec.size())-1;++i)
    {
      if(fvec[i]==fvec[i+1])
      {
        total++;
        tri::Allocator<MeshType>::DeleteFace(m, *(fvec[i].fp) );
      }
    }
    return total;
  }

  /** This function removes all duplicate faces of the mesh by looking only at their vertex reference.
            So it should be called after unification of vertices.
            Note that it does not update any topology relation that could be affected by this like the VT or TT relation.
            the reason this function is usually performed BEFORE building any topology information.
            */
  static int RemoveDuplicateEdge( MeshType & m)    // V1.0
  {
    if (m.en==0) return 0;
    std::vector<SortedPair> eVec;
    for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
      if(!(*ei).IsD())
      {
        eVec.push_back(SortedPair(	tri::Index(m,(*ei).V(0)), tri::Index(m,(*ei).V(1)), &*ei));
      }
    std::sort(eVec.begin(),eVec.end());
    int total=0;
    for(int i=0;i<int(eVec.size())-1;++i)
    {
      if(eVec[i]==eVec[i+1])
      {
        total++;
        tri::Allocator<MeshType>::DeleteEdge(m, *(eVec[i].fp) );
      }
    }
    return total;
  }

  static int CountUnreferencedVertex( MeshType& m)
  {
    return RemoveUnreferencedVertex(m,false);
  }


  /** This function removes that are not referenced by any face. The function updates the vn counter.
            @param m The mesh
            @return The number of removed vertices
            */
  static int RemoveUnreferencedVertex( MeshType& m, bool DeleteVertexFlag=true)   // V1.0
  {
    FaceIterator fi;
    EdgeIterator ei;
    VertexIterator vi;
    int referredBit = VertexType::NewBitFlag();

    int j;
    int deleted = 0;

    for(vi=m.vert.begin();vi!=m.vert.end();++vi)
      (*vi).ClearUserBit(referredBit);

    for(fi=m.face.begin();fi!=m.face.end();++fi)
      if( !(*fi).IsD() )
        for(j=0;j<(*fi).VN();++j)
          (*fi).V(j)->SetUserBit(referredBit);

    for(ei=m.edge.begin();ei!=m.edge.end();++ei)
      if( !(*ei).IsD() ){
        (*ei).V(0)->SetUserBit(referredBit);
        (*ei).V(1)->SetUserBit(referredBit);
      }

    for(vi=m.vert.begin();vi!=m.vert.end();++vi)
      if( (!(*vi).IsD()) && (!(*vi).IsUserBit(referredBit)))
      {
        if(DeleteVertexFlag) Allocator<MeshType>::DeleteVertex(m,*vi);
        ++deleted;
      }
    VertexType::DeleteBitFlag(referredBit);
    return deleted;
  }

  /**
      Degenerate vertices are vertices that have coords with invalid floating point values,
      All the faces incident on deleted vertices are also deleted
            */
  static int RemoveDegenerateVertex(MeshType& m)
  {
    VertexIterator vi;
    int count_vd = 0;

    for(vi=m.vert.begin(); vi!=m.vert.end();++vi)
      if(math::IsNAN( (*vi).P()[0]) ||
         math::IsNAN( (*vi).P()[1]) ||
         math::IsNAN( (*vi).P()[2]) )
      {
        count_vd++;
        Allocator<MeshType>::DeleteVertex(m,*vi);
      }

    FaceIterator fi;
    int count_fd = 0;

    for(fi=m.face.begin(); fi!=m.face.end();++fi)
      if(!(*fi).IsD())
        if( (*fi).V(0)->IsD() ||
            (*fi).V(1)->IsD() ||
            (*fi).V(2)->IsD() )
        {
          count_fd++;
          Allocator<MeshType>::DeleteFace(m,*fi);
        }
    return count_vd;
  }

  /**
      Degenerate faces are faces that are Topologically degenerate,
      i.e. have two or more vertex reference that link the same vertex
      (and not only two vertexes with the same coordinates).
      All Degenerate faces are zero area faces BUT not all zero area faces are degenerate.
      We do not take care of topology because when we have degenerate faces the
      topology calculation functions crash.
      */
  static int RemoveDegenerateFace(MeshType& m)
  {
    int count_fd = 0;

    for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
      if(!(*fi).IsD())
      {
        if((*fi).V(0) == (*fi).V(1) ||
           (*fi).V(0) == (*fi).V(2) ||
           (*fi).V(1) == (*fi).V(2) )
        {
          count_fd++;
          Allocator<MeshType>::DeleteFace(m,*fi);
        }
      }
    return count_fd;
  }

  static int RemoveDegenerateEdge(MeshType& m)
  {
    int count_ed = 0;

    for(EdgeIterator ei=m.edge.begin(); ei!=m.edge.end();++ei)
      if(!(*ei).IsD())
      {
        if((*ei).V(0) == (*ei).V(1) )
        {
          count_ed++;
          Allocator<MeshType>::DeleteEdge(m,*ei);
        }
      }
    return count_ed;
  }

  static int RemoveNonManifoldVertex(MeshType& m)
  {
    CountNonManifoldVertexFF(m,true);
    tri::UpdateSelection<MeshType>::FaceFromVertexLoose(m);
    int count_removed = 0;
    for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
      if(!(*fi).IsD() && (*fi).IsS())
        Allocator<MeshType>::DeleteFace(m,*fi);
    for(VertexIterator vi=m.vert.begin(); vi!=m.vert.end();++vi)
      if(!(*vi).IsD() && (*vi).IsS()) {
        ++count_removed;
        Allocator<MeshType>::DeleteVertex(m,*vi);
      }
    return count_removed;
  }

  static int SplitSelectedVertexOnEdgeMesh(MeshType& m)
  {
    tri::RequireCompactness(m);
    tri::UpdateFlags<MeshType>::VertexClearV(m);
    int count_split = 0;
    for(size_t i=0;i<m.edge.size();++i)
    {
      for(int j=0;j<2;++j)
      {
        VertexPointer vp = m.edge[i].V(j);
        if(vp->IsS())
        {
          if(!vp->IsV())
	    {
            m.edge[i].V(j) = &*(tri::Allocator<MeshType>::AddVertex(m,vp->P()));
	    ++count_split;
	    }
          else 
	    {
	      vp->SetV();
	    }
	  
        }
      }
    }
    return count_split;
  }


  static void SelectNonManifoldVertexOnEdgeMesh(MeshType &m)
  {
    tri::RequireCompactness(m);
    tri::UpdateSelection<MeshType>::VertexClear(m);
    std::vector<int> cnt(m.vn,0);

    for(size_t i=0;i<m.edge.size();++i)
    {
      cnt[tri::Index(m,m.edge[i].V(0))]++;
      cnt[tri::Index(m,m.edge[i].V(1))]++;
    }
    for(size_t i=0;i<m.vert.size();++i)
      if(cnt[i]>2) m.vert[i].SetS();
  }

  static void SelectCreaseVertexOnEdgeMesh(MeshType &m, ScalarType AngleRadThr)
  {
    tri::RequireCompactness(m);
    tri::RequireVEAdjacency(m);
    tri::UpdateTopology<MeshType>::VertexEdge(m);
    for(size_t i=0;i<m.vert.size();++i)
    {
      std::vector<VertexPointer> VVStarVec;
      edge::VVStarVE(&(m.vert[i]),VVStarVec);
      if(VVStarVec.size()==2)
      {
        CoordType v0 = m.vert[i].P() - VVStarVec[0]->P();
        CoordType v1 = m.vert[i].P() - VVStarVec[1]->P();
        float angle = M_PI-vcg::Angle(v0,v1);
        if(angle > AngleRadThr) m.vert[i].SetS();
      }
    }
  }


  /// Removal of faces that were incident on a non manifold edge.

  // Given a mesh with FF adjacency
  // it search for non manifold vertices and duplicate them.
  // Duplicated vertices are moved apart according to the move threshold param.
  // that is a percentage of the average vector from the non manifold vertex to the barycenter of the incident faces.

  static int SplitNonManifoldVertex(MeshType& m, ScalarType moveThreshold)
  {
    RequireFFAdjacency(m);
    typedef std::pair<FacePointer,int> FaceInt; // a face and the index of the vertex that we have to change
    //
    std::vector<std::pair<VertexPointer, std::vector<FaceInt> > >ToSplitVec;

    SelectionStack<MeshType> ss(m);
    ss.push();
    CountNonManifoldVertexFF(m,true);
    UpdateFlags<MeshType>::VertexClearV(m);
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)	if (!fi->IsD())
    {
      for(int i=0;i<3;i++)
        if((*fi).V(i)->IsS() && !(*fi).V(i)->IsV())
        {
          (*fi).V(i)->SetV();
          face::Pos<FaceType> startPos(&*fi,i);
          face::Pos<FaceType> curPos = startPos;
          std::set<FaceInt> faceSet;
          do
          {
            faceSet.insert(std::make_pair(curPos.F(),curPos.VInd()));
            curPos.NextE();
          } while (curPos != startPos);

          ToSplitVec.push_back(make_pair((*fi).V(i),std::vector<FaceInt>()));

          typename std::set<FaceInt>::const_iterator iii;

          for(iii=faceSet.begin();iii!=faceSet.end();++iii)
            ToSplitVec.back().second.push_back(*iii);
        }
    }
    ss.pop();
    // Second step actually add new vertices and split them.
    typename tri::Allocator<MeshType>::template PointerUpdater<VertexPointer> pu;
    VertexIterator firstVp = tri::Allocator<MeshType>::AddVertices(m,ToSplitVec.size(),pu);
    for(size_t i =0;i<ToSplitVec.size();++i)
    {
      //          qDebug("Splitting Vertex %i",ToSplitVec[i].first-&*m.vert.begin());
      VertexPointer np=ToSplitVec[i].first;
      pu.Update(np);
      firstVp->ImportData(*np);
      // loop on the face to be changed, and also compute the movement vector;
      CoordType delta(0,0,0);
      for(size_t j=0;j<ToSplitVec[i].second.size();++j)
      {
        FaceInt ff=ToSplitVec[i].second[j];
        ff.first->V(ff.second)=&*firstVp;
        delta+=Barycenter(*(ff.first))-np->cP();
      }
      delta /= ToSplitVec[i].second.size();
      firstVp->P() = firstVp->P() + delta * moveThreshold;
      firstVp++;
    }

    return ToSplitVec.size();
  }


  // Auxiliary function for sorting the non manifold faces according to their area. Used in  RemoveNonManifoldFace
  struct CompareAreaFP {
    bool operator ()(FacePointer const& f1, FacePointer const& f2) const {
      return DoubleArea(*f1) < DoubleArea(*f2);
    }
  };

  /// Removal of faces that were incident on a non manifold edge.
  static int RemoveNonManifoldFace(MeshType& m)
  {
    FaceIterator fi;
    int count_fd = 0;
    std::vector<FacePointer> ToDelVec;

    for(fi=m.face.begin(); fi!=m.face.end();++fi)
      if (!fi->IsD())
      {
        if ((!IsManifold(*fi,0))||
            (!IsManifold(*fi,1))||
            (!IsManifold(*fi,2)))
          ToDelVec.push_back(&*fi);
      }

    std::sort(ToDelVec.begin(),ToDelVec.end(),CompareAreaFP());

    for(size_t i=0;i<ToDelVec.size();++i)
    {
      if(!ToDelVec[i]->IsD())
      {
        FaceType &ff= *ToDelVec[i];
        if ((!IsManifold(ff,0))||
            (!IsManifold(ff,1))||
            (!IsManifold(ff,2)))
        {
          for(int j=0;j<3;++j)
            if(!face::IsBorder<FaceType>(ff,j))
              vcg::face::FFDetach<FaceType>(ff,j);

          Allocator<MeshType>::DeleteFace(m,ff);
          count_fd++;
        }
      }
    }
    return count_fd;
  }

  /* Remove the faces that are out of a given range of area  */
  static int RemoveFaceOutOfRangeArea(MeshType& m, ScalarType MinAreaThr=0, ScalarType MaxAreaThr=(std::numeric_limits<ScalarType>::max)(), bool OnlyOnSelected=false)
  {
    int count_fd = 0;
    MinAreaThr*=2;
    MaxAreaThr*=2;
    for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi){
      if(!(*fi).IsD())
        if(!OnlyOnSelected || (*fi).IsS())
        {
          const ScalarType doubleArea=DoubleArea<FaceType>(*fi);
          if((doubleArea<=MinAreaThr) || (doubleArea>=MaxAreaThr) )
          {
            Allocator<MeshType>::DeleteFace(m,*fi);
            count_fd++;
          }
        }
    }
    return count_fd;
  }

  static int RemoveZeroAreaFace(MeshType& m) { return RemoveFaceOutOfRangeArea(m,0);}

  

  /**
             * Is the mesh only composed by quadrilaterals?
             */
  static bool IsBitQuadOnly(const MeshType &m)
  {
    typedef typename MeshType::FaceType F;
    tri::RequirePerFaceFlags(m);
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD()) {
      unsigned int tmp = fi->Flags()&(F::FAUX0|F::FAUX1|F::FAUX2);
      if ( tmp != F::FAUX0 && tmp != F::FAUX1 && tmp != F::FAUX2) return false;
    }
    return true;
  }


  static bool IsFaceFauxConsistent(MeshType &m)
  {
    RequirePerFaceFlags(m);
    RequireFFAdjacency(m);
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
    {
      for(int z=0;z<(*fi).VN();++z)
      {
        FacePointer fp = fi->FFp(z);
        int zp = fi->FFi(z);
        if(fi->IsF(z) != fp->IsF(zp)) return false;
      }
    }
    return true;
  }

  /**
* Is the mesh only composed by triangles? (non polygonal faces)
*/
  static bool IsBitTriOnly(const MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if ( !fi->IsD()  &&  fi->IsAnyF() ) return false;
    }
    return true;
  }

  static bool IsBitPolygonal(const MeshType &m){
    return !IsBitTriOnly(m);
  }

  /**
   * Is the mesh only composed by quadrilaterals and triangles? (no pentas, etc)
   * It assumes that the bits are consistent. In that case there can be only a single faux edge.
   */
  static bool IsBitTriQuadOnly(const MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    typedef typename MeshType::FaceType F;
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD()) {
      unsigned int tmp = fi->cFlags()&(F::FAUX0|F::FAUX1|F::FAUX2);
      if ( tmp!=F::FAUX0 && tmp!=F::FAUX1 && tmp!=F::FAUX2 && tmp!=0 ) return false;
    }
    return true;
  }

  /**
   * How many quadrilaterals?
   * It assumes that the bits are consistent. In that case we count the tris with a single faux edge and divide by two.
   */
  static int CountBitQuads(const MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    typedef typename MeshType::FaceType F;
    int count=0;
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD()) {
      unsigned int tmp = fi->cFlags()&(F::FAUX0|F::FAUX1|F::FAUX2);
      if ( tmp==F::FAUX0 || tmp==F::FAUX1 || tmp==F::FAUX2) count++;
    }
    return count / 2;
  }

  /**
   * How many triangles? (non polygonal faces)
   */
  static int CountBitTris(const MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    int count=0;
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD()) {
      if (!(fi->IsAnyF())) count++;
    }
    return count;
  }

  /**
   * How many polygons of any kind? (including triangles)
   * it assumes that there are no faux vertexes (e.g vertices completely surrounded by faux edges)
   */
  static int CountBitPolygons(const MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    int count = 0;
    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD())  {
      if (fi->IsF(0)) count++;
      if (fi->IsF(1)) count++;
      if (fi->IsF(2)) count++;
    }
    return m.fn - count/2;
  }

  /**
  * The number of polygonal faces is
  *  FN - EN_f (each faux edge hides exactly one triangular face or in other words a polygon of n edges has n-3 faux edges.)
  * In the general case where a The number of polygonal faces is
  *	 FN - EN_f + VN_f
  *	where:
  *	 EN_f is the number of faux edges.
  *	 VN_f is the number of faux vertices (e.g vertices completely surrounded by faux edges)
  * as a intuitive proof think to a internal vertex that is collapsed onto a border of a polygon:
  * it deletes 2 faces, 1 faux edges and 1 vertex so to keep the balance you have to add back the removed vertex.
  */
  static int CountBitLargePolygons(MeshType &m)
  {
    tri::RequirePerFaceFlags(m);
    UpdateFlags<MeshType>::VertexSetV(m);
    // First loop Clear all referenced vertices
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if (!fi->IsD())
        for(int i=0;i<3;++i) fi->V(i)->ClearV();


    // Second Loop, count (twice) faux edges and mark all vertices touched by non faux edges
    // (e.g vertexes on the boundary of a polygon)
    int countE = 0;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if (!fi->IsD())  {
        for(int i=0;i<3;++i)
        {
          if (fi->IsF(i))
            countE++;
          else
          {
            fi->V0(i)->SetV();
            fi->V1(i)->SetV();
          }
        }
      }
    // Third Loop, count the number of referenced vertexes that are completely surrounded by faux edges.

    int countV = 0;
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if (!vi->IsD() && !vi->IsV()) countV++;

    return m.fn - countE/2 + countV ;
  }


  /**
  * Checks that the mesh has consistent per-face faux edges
  * (the ones that merges triangles into larger polygons).
  * A border edge should never be faux, and faux edges should always be
  * reciprocated by another faux edges.
  * It requires FF adjacency.
  */
  static bool HasConsistentPerFaceFauxFlag(const MeshType &m)
  {
    RequireFFAdjacency(m);
    RequirePerFaceFlags(m);

    for (ConstFaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
        for (int k=0; k<3; k++)
          if( ( fi->IsF(k) != fi->cFFp(k)->IsF(fi->cFFi(k)) ) ||
              ( fi->IsF(k) && face::IsBorder(*fi,k)) )
          {
            return false;
          }
    return true;
  }

  /**
   * Count the number of non manifold edges in a polylinemesh, e.g. the edges where there are more than 2 incident faces.
   *
   */
  static int CountNonManifoldEdgeEE( MeshType & m, bool SelectFlag=false)
  {
    MeshAssert<MeshType>::OnlyEdgeMesh(m);
    RequireEEAdjacency(m);
    tri::UpdateTopology<MeshType>::EdgeEdge(m);

    if(SelectFlag) UpdateSelection<MeshType>::VertexClear(m);

    int nonManifoldCnt=0;
    SimpleTempData<typename MeshType::VertContainer, int > TD(m.vert,0);

    // First Loop, just count how many faces are incident on a vertex and store it in the TemporaryData Counter.
    EdgeIterator ei;
    for (ei = m.edge.begin(); ei != m.edge.end(); ++ei)	if (!ei->IsD())
    {
      TD[(*ei).V(0)]++;
      TD[(*ei).V(1)]++;
    }

    tri::UpdateFlags<MeshType>::VertexClearV(m);
    // Second Loop, Check that each vertex have been seen 1 or 2 times.
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)	if (!vi->IsD())
    {
      if( TD[vi] >2 )
      {
        if(SelectFlag) (*vi).SetS();
        nonManifoldCnt++;
      }
    }
    return nonManifoldCnt;
  }

  /**
       * Count the number of non manifold edges in a mesh, e.g. the edges where there are more than 2 incident faces.
       *
       * Note that this test is not enough to say that a mesh is two manifold,
       * you have to count also the non manifold vertexes.
       */
  static int CountNonManifoldEdgeFF( MeshType & m, bool SelectFlag=false)
  {
    RequireFFAdjacency(m);
    int nmfBit[3];
    nmfBit[0]= FaceType::NewBitFlag();
    nmfBit[1]= FaceType::NewBitFlag();
    nmfBit[2]= FaceType::NewBitFlag();


    UpdateFlags<MeshType>::FaceClear(m,nmfBit[0]+nmfBit[1]+nmfBit[2]);

    if(SelectFlag){
      UpdateSelection<MeshType>::VertexClear(m);
      UpdateSelection<MeshType>::FaceClear(m);
    }

    int edgeCnt = 0;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!fi->IsD())
      {
        for(int i=0;i<3;++i)
          if(!IsManifold(*fi,i))
          {
            if(!(*fi).IsUserBit(nmfBit[i]))
            {
              ++edgeCnt;
              if(SelectFlag)
              {
                (*fi).V0(i)->SetS();
                (*fi).V1(i)->SetS();
              }
              // follow the ring of faces incident on edge i;
              face::Pos<FaceType> nmf(&*fi,i);
              do
              {
                if(SelectFlag) nmf.F()->SetS();
                nmf.F()->SetUserBit(nmfBit[nmf.E()]);
                nmf.NextF();
              }
              while(nmf.f != &*fi);
            }
          }
      }
    }
    return edgeCnt;
  }

  /** Count (and eventually select) non 2-Manifold vertexes of a mesh
       * e.g. the vertices with a non 2-manif. neighbourhood but that do not belong to not 2-manif edges.
       * typical situation two cones connected by one vertex.
       */
  static int CountNonManifoldVertexFF( MeshType & m, bool selectVert = true )
  {
    RequireFFAdjacency(m);
    if(selectVert) UpdateSelection<MeshType>::VertexClear(m);

    int nonManifoldCnt=0;
    SimpleTempData<typename MeshType::VertContainer, int > TD(m.vert,0);

    // First Loop, just count how many faces are incident on a vertex and store it in the TemporaryData Counter.
    FaceIterator fi;
    for (fi = m.face.begin(); fi != m.face.end(); ++fi)	if (!fi->IsD())
    {
      TD[(*fi).V(0)]++;
      TD[(*fi).V(1)]++;
      TD[(*fi).V(2)]++;
    }

    tri::UpdateFlags<MeshType>::VertexClearV(m);
    // Second Loop.
    // mark out of the game the vertexes that are incident on non manifold edges.
    for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!fi->IsD())
    {
      for(int i=0;i<3;++i)
        if (!IsManifold(*fi,i))  {
          (*fi).V0(i)->SetV();
          (*fi).V1(i)->SetV();
        }
    }
    // Third Loop, for safe vertexes, check that the number of faces that you can reach starting
    // from it and using FF is the same of the previously counted.
    for (fi = m.face.begin(); fi != m.face.end(); ++fi)	if (!fi->IsD())
    {
      for(int i=0;i<3;i++) if(!(*fi).V(i)->IsV()){
        (*fi).V(i)->SetV();
        face::Pos<FaceType> pos(&(*fi),i);

        int starSizeFF = pos.NumberOfIncidentFaces();

        if (starSizeFF != TD[(*fi).V(i)])
        {
          if(selectVert) (*fi).V(i)->SetS();
          nonManifoldCnt++;
        }
      }
    }
    return nonManifoldCnt;
  }
  /// Very simple test of water tightness. No boundary and no non manifold edges. 
  /// Assume that it is orientable. 
  /// It could be debated if a closed non orientable surface is watertight or not. 
  /// 
  /// The rationale of not testing orientability here is that 
  /// it requires FFAdj while this test do not require any adjacency.
  /// 
  static bool IsWaterTight(MeshType & m)
  {
    int edgeNum=0,edgeBorderNum=0,edgeNonManifNum=0;
    CountEdgeNum(m, edgeNum, edgeBorderNum,edgeNonManifNum);
    return  (edgeBorderNum==0) && (edgeNonManifNum==0);
  }

  static void CountEdgeNum( MeshType & m, int &total_e, int &boundary_e, int &non_manif_e )
  {
    std::vector< typename tri::UpdateTopology<MeshType>::PEdge > edgeVec;
    tri::UpdateTopology<MeshType>::FillEdgeVector(m,edgeVec,true);
    sort(edgeVec.begin(), edgeVec.end());		// Lo ordino per vertici
    total_e=0;
    boundary_e=0;
    non_manif_e=0;

    size_t f_on_cur_edge =1;
    for(size_t i=0;i<edgeVec.size();++i)
    {
      if(( (i+1) == edgeVec.size()) ||  !(edgeVec[i] == edgeVec[i+1]))
      {
        ++total_e;
        if(f_on_cur_edge==1)
          ++boundary_e;
        if(f_on_cur_edge>2)
          ++non_manif_e;
        f_on_cur_edge=1;
      }
      else
      {
        ++f_on_cur_edge;
      }
    } // end for
  }



  static int CountHoles( MeshType & m)
  {
    UpdateFlags<MeshType>::FaceClearV(m);
    int loopNum=0;
    for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi) if(!fi->IsD())
    {
        for(int j=0;j<3;++j)
        {
          if(!fi->IsV() && face::IsBorder(*fi,j))
          {
            face::Pos<FaceType> startPos(&*fi,j);
            face::Pos<FaceType> curPos=startPos;
            do
            {
              curPos.NextB();
              curPos.F()->SetV();
            }
            while(curPos!=startPos);          
            ++loopNum;
          }
        }
    }
    return loopNum;
  }

  /*
  Compute the set of connected components of a given mesh
  it fills a vector of pair < int , faceptr > with, for each connecteed component its size and a represnant
 */
  static int CountConnectedComponents(MeshType &m)
  {
    std::vector< std::pair<int,FacePointer> > CCV;
    return ConnectedComponents(m,CCV);
  }

  static int ConnectedComponents(MeshType &m, std::vector< std::pair<int,FacePointer> > &CCV)
  {
    tri::RequireFFAdjacency(m);
    CCV.clear();
    tri::UpdateFlags<MeshType>::FaceClearV(m);
    std::stack<FacePointer> sf;
    FacePointer fpt=&*(m.face.begin());
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      if(!((*fi).IsD()) && !(*fi).IsV())
      {
        (*fi).SetV();
        CCV.push_back(std::make_pair(0,&*fi));
        sf.push(&*fi);
        while (!sf.empty())
        {
          fpt=sf.top();
          ++CCV.back().first;
          sf.pop();
          for(int j=0;j<3;++j)
          {
            if( !face::IsBorder(*fpt,j) )
            {
              FacePointer l = fpt->FFp(j);
              if( !(*l).IsV() )
              {
                (*l).SetV();
                sf.push(l);
              }
            }
          }
        }
      }
    }
    return int(CCV.size());
  }

  static void ComputeValence( MeshType &m, typename MeshType::PerVertexIntHandle &h)
  {
    for(VertexIterator vi=m.vert.begin(); vi!= m.vert.end();++vi)
      h[vi]=0;

    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      if(!((*fi).IsD()))
        for(int j=0;j<fi->VN();j++)
          ++h[tri::Index(m,fi->V(j))];
    }
  }

  /**
            GENUS.

            A topologically invariant property of a surface defined as
            the largest number of non-intersecting simple closed curves that can be
            drawn on the surface without separating it.

      Roughly speaking, it is the number of holes in a surface.
            The genus g of a closed surface, also called the geometric genus, is related to the
            Euler characteristic by the relation $chi$ by $chi==2-2g$.

            The genus of a connected, orientable surface is an integer representing the maximum
            number of cuttings along closed simple curves without rendering the resultant
            manifold disconnected. It is equal to the number of handles on it.

            For general polyhedra the <em>Euler Formula</em> is:

                  V - E + F = 2 - 2G - B

            where V is the number of vertices, F is the number of faces, E is the
            number of edges, G is the genus and B is the number of <em>boundary polygons</em>.

            The above formula is valid for a mesh with one single connected component.
            By considering multiple connected components the formula becomes:

                  V - E + F = 2C - 2Gs - B   ->   2Gs = - ( V-E+F +B -2C)

            where C is the number of connected components and Gs is the sum of
            the genus of all connected components.

            Note that in the case of a mesh with boundaries the intuitive meaning of Genus is less intuitive that it could seem.
            A closed sphere, a sphere with one hole (e.g. a disk) and a sphere with two holes (e.g. a tube) all of them have Genus == 0

            */

  static int MeshGenus(int nvert,int nedges,int nfaces, int numholes, int numcomponents)
  {
    return -((nvert + nfaces - nedges + numholes - 2 * numcomponents) / 2);
  }

  static int MeshGenus(MeshType &m)
  {
    int nvert=m.vn;
    int nfaces=m.fn;
    int boundary_e,total_e,nonmanif_e;
    CountEdgeNum(m,total_e,boundary_e,nonmanif_e);
    int numholes=CountHoles(m);
    int numcomponents=CountConnectedComponents(m);
    int G=MeshGenus(nvert,total_e,nfaces,numholes,numcomponents);
    return G;
  }

  /**
             * Check if the given mesh is regular, semi-regular or irregular.
             *
             * Each vertex of a \em regular mesh has valence 6 except for border vertices
             * which have valence 4.
             *
             * A \em semi-regular mesh is derived from an irregular one applying
             * 1-to-4 subdivision recursively. (not checked for now)
             *
             * All other meshes are \em irregular.
             */
  static void IsRegularMesh(MeshType &m, bool &Regular, bool &Semiregular)
  {
    RequireVFAdjacency(m);
    Regular = true;

    VertexIterator vi;

    // for each vertex the number of edges are count
    for (vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    {
      if (!vi->IsD())
      {
        face::Pos<FaceType> he((*vi).VFp(), &*vi);
        face::Pos<FaceType> ht = he;

        int n=0;
        bool border=false;
        do
        {
          ++n;
          ht.NextE();
          if (ht.IsBorder())
            border=true;
        }
        while (ht != he);

        if (border)
          n = n/2;

        if ((n != 6)&&(!border && n != 4))
        {
          Regular = false;
          break;
        }
      }
    }

    if (!Regular)
      Semiregular = false;
    else
    {
      // For now we do not account for semi-regularity
      Semiregular = false;
    }
  }


  static bool IsCoherentlyOrientedMesh(MeshType &m)
  {
    RequireFFAdjacency(m);
    MeshAssert<MeshType>::FFAdjacencyIsInitialized(m);       
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if (!fi->IsD())
        for(int i=0;i<3;++i)
          if(!face::CheckOrientation(*fi,i))
            return false;

    return true;
  }

  static void OrientCoherentlyMesh(MeshType &m, bool &_IsOriented, bool &_IsOrientable)
  {
    RequireFFAdjacency(m);
    MeshAssert<MeshType>::FFAdjacencyIsInitialized(m);   
    bool IsOrientable = true;
    bool IsOriented = true;

    UpdateFlags<MeshType>::FaceClearV(m);
    std::stack<FacePointer> faces;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!fi->IsD() && !fi->IsV())
      {
        // each face put in the stack is selected (and oriented)
        fi->SetV();
        faces.push(&(*fi));
        while (!faces.empty())
        {
          FacePointer fp = faces.top();
          faces.pop();

          // make consistently oriented the adjacent faces
          for (int j = 0; j < 3; j++)
          {
            if (!face::IsBorder(*fp,j) && face::IsManifold<FaceType>(*fp, j))
            {
              FacePointer fpaux = fp->FFp(j);
              int iaux = fp->FFi(j);
              if (!CheckOrientation(*fpaux, iaux))
              {
                IsOriented = false;

                if (!fpaux->IsV())
                  face::SwapEdge<FaceType,true>(*fpaux, iaux);
                else
                {
                  IsOrientable = false;
                  break;
                }
              }
              if (!fpaux->IsV())
              {
                fpaux->SetV();
                faces.push(fpaux);
              }
            }
          }
        }
      }
      if (!IsOrientable)	break;
    }
    _IsOriented = IsOriented;
    _IsOrientable = IsOrientable;
  }


  /// Flip the orientation of the whole mesh flipping all the faces (by swapping the first two vertices)
  static void FlipMesh(MeshType &m, bool selected=false)
  {
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) if(!(*fi).IsD())
      if(!selected || (*fi).IsS())
      {
        face::SwapEdge<FaceType,false>((*fi), 0);
        if (HasPerWedgeTexCoord(m))
          std::swap((*fi).WT(0),(*fi).WT(1));
      }
  }
  /// Flip a mesh so that its normals are orented outside.
  /// Just for safety it uses a voting scheme.
  /// It assumes that
  /// mesh has already has coherent normals.
  /// mesh is watertight and signle component.
  static bool FlipNormalOutside(MeshType &m)
  {
    if(m.vert.empty()) return false;

    tri::UpdateNormal<MeshType>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<MeshType>::NormalizePerVertex(m);

    std::vector< VertexPointer > minVertVec;
    std::vector< VertexPointer > maxVertVec;

    // The set of directions to be choosen
    std::vector< CoordType > dirVec;
    dirVec.push_back(CoordType(1,0,0));
    dirVec.push_back(CoordType(0,1,0));
    dirVec.push_back(CoordType(0,0,1));
    dirVec.push_back(CoordType( 1, 1,1));
    dirVec.push_back(CoordType(-1, 1,1));
    dirVec.push_back(CoordType(-1,-1,1));
    dirVec.push_back(CoordType( 1,-1,1));
    for(size_t i=0;i<dirVec.size();++i)
    {
      Normalize(dirVec[i]);
      minVertVec.push_back(&*m.vert.begin());
      maxVertVec.push_back(&*m.vert.begin());
    }
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) if(!(*vi).IsD())
    {
      for(size_t i=0;i<dirVec.size();++i)
      {
        if( (*vi).cP().dot(dirVec[i]) < minVertVec[i]->P().dot(dirVec[i])) minVertVec[i] = &*vi;
        if( (*vi).cP().dot(dirVec[i]) > maxVertVec[i]->P().dot(dirVec[i])) maxVertVec[i] = &*vi;
      }
    }

    int voteCount=0;
    ScalarType angleThreshold = cos(math::ToRad(85.0));
    for(size_t i=0;i<dirVec.size();++i)
    {
      //          qDebug("Min vert along (%f %f %f) is %f %f %f",dirVec[i][0],dirVec[i][1],dirVec[i][2],minVertVec[i]->P()[0],minVertVec[i]->P()[1],minVertVec[i]->P()[2]);
      //          qDebug("Max vert along (%f %f %f) is %f %f %f",dirVec[i][0],dirVec[i][1],dirVec[i][2],maxVertVec[i]->P()[0],maxVertVec[i]->P()[1],maxVertVec[i]->P()[2]);
      if(minVertVec[i]->N().dot(dirVec[i]) > angleThreshold ) voteCount++;
      if(maxVertVec[i]->N().dot(dirVec[i]) < -angleThreshold ) voteCount++;
    }
    //        qDebug("votecount = %i",voteCount);
    if(voteCount < int(dirVec.size())/2) return false;
    FlipMesh(m);
    return true;
  }

  // Search and remove small single triangle folds
  // - a face has normal opposite to all other faces
  // - choose the edge that brings to the face f1 containing the vertex opposite to that edge.
  static int RemoveFaceFoldByFlip(MeshType &m, float normalThresholdDeg=175, bool repeat=true)
  {
    RequireFFAdjacency(m);
    RequirePerVertexMark(m);
    //Counters for logging and convergence
    int count, total = 0;

    do {
      tri::UpdateTopology<MeshType>::FaceFace(m);
      tri::UnMarkAll(m);
      count = 0;

      ScalarType NormalThrRad = math::ToRad(normalThresholdDeg);
      ScalarType eps = 0.0001; // this epsilon value is in absolute value. It is a distance from edge in baricentric coords.
      //detection stage
      for(FaceIterator fi=m.face.begin();fi!= m.face.end();++fi ) if(!(*fi).IsV())
      { Point3<ScalarType> NN = vcg::TriangleNormal((*fi)).Normalize();
        if( vcg::AngleN(NN,TriangleNormal(*(*fi).FFp(0)).Normalize()) > NormalThrRad &&
            vcg::AngleN(NN,TriangleNormal(*(*fi).FFp(1)).Normalize()) > NormalThrRad &&
            vcg::AngleN(NN,TriangleNormal(*(*fi).FFp(2)).Normalize()) > NormalThrRad )
        {
          (*fi).SetS();
          //(*fi).C()=Color4b(Color4b::Red);
          // now search the best edge to flip
          for(int i=0;i<3;i++)
          {
            Point3<ScalarType> &p=(*fi).P2(i);
            Point3<ScalarType> L;
            bool ret = vcg::InterpolationParameters((*(*fi).FFp(i)),TriangleNormal(*(*fi).FFp(i)),p,L);
            if(ret && L[0]>eps && L[1]>eps && L[2]>eps)
            {
              (*fi).FFp(i)->SetS();
              (*fi).FFp(i)->SetV();
              //(*fi).FFp(i)->C()=Color4b(Color4b::Green);
              if(face::CheckFlipEdge<FaceType>( *fi, i ))  {
                face::FlipEdge<FaceType>( *fi, i );
                ++count; ++total;
              }
            }
          }
        }
      }

      // tri::UpdateNormal<MeshType>::PerFace(m);
    }
    while( repeat && count );
    return total;
  }


  static int RemoveTVertexByFlip(MeshType &m, float threshold=40, bool repeat=true)
  {
    RequireFFAdjacency(m);
    RequirePerVertexMark(m);
    //Counters for logging and convergence
    int count, total = 0;

    do {
      tri::UpdateTopology<MeshType>::FaceFace(m);
      tri::UnMarkAll(m);
      count = 0;

      //detection stage
      for(unsigned int index = 0 ; index < m.face.size(); ++index )
      {
        FacePointer f = &(m.face[index]);    float sides[3]; CoordType dummy;
        sides[0] = Distance(f->P(0), f->P(1));
        sides[1] = Distance(f->P(1), f->P(2));
        sides[2] = Distance(f->P(2), f->P(0));
        // Find largest triangle side
        int i = std::find(sides, sides+3, std::max( std::max(sides[0],sides[1]), sides[2])) - (sides);
        if( tri::IsMarked(m,f->V2(i) )) continue;

        if( PSDist(f->P2(i),f->P(i),f->P1(i),dummy)*threshold <= sides[i] )
        {
          tri::Mark(m,f->V2(i));
          if(face::CheckFlipEdge<FaceType>( *f, i ))  {
            // Check if EdgeFlipping improves quality
            FacePointer g = f->FFp(i); int k = f->FFi(i);
            Triangle3<ScalarType> t1(f->P(i), f->P1(i), f->P2(i)), t2(g->P(k), g->P1(k), g->P2(k)),
                t3(f->P(i), g->P2(k), f->P2(i)), t4(g->P(k), f->P2(i), g->P2(k));

            if ( std::min( QualityFace(t1), QualityFace(t2) ) < std::min( QualityFace(t3), QualityFace(t4) ))
            {
              face::FlipEdge<FaceType>( *f, i );
              ++count; ++total;
            }
          }

        }
      }

      // tri::UpdateNormal<MeshType>::PerFace(m);
    }
    while( repeat && count );
    return total;
  }

  static int RemoveTVertexByCollapse(MeshType &m, float threshold=40, bool repeat=true)
  {
    RequirePerVertexMark(m);
    //Counters for logging and convergence
    int count, total = 0;

    do {
      tri::UnMarkAll(m);
      count = 0;

      //detection stage
      for(unsigned int index = 0 ; index < m.face.size(); ++index )
      {
        FacePointer f = &(m.face[index]);
        float sides[3];
        CoordType dummy;

        sides[0] = Distance(f->P(0), f->P(1));
        sides[1] = Distance(f->P(1), f->P(2));
        sides[2] = Distance(f->P(2), f->P(0));
        int i = std::find(sides, sides+3, std::max( std::max(sides[0],sides[1]), sides[2])) - (sides);
        if( tri::IsMarked(m,f->V2(i) )) continue;

        if( PSDist(f->P2(i),f->P(i),f->P1(i),dummy)*threshold <= sides[i] )
        {
          tri::Mark(m,f->V2(i));

          int j = Distance(dummy,f->P(i))<Distance(dummy,f->P1(i))?i:(i+1)%3;
          f->P2(i) = f->P(j);  tri::Mark(m,f->V(j));
          ++count; ++total;
        }
      }


      tri::Clean<MeshType>::RemoveDuplicateVertex(m);
      tri::Allocator<MeshType>::CompactFaceVector(m);
      tri::Allocator<MeshType>::CompactVertexVector(m);
    }
    while( repeat && count );

    return total;
  }

  static bool SelfIntersections(MeshType &m, std::vector<FaceType*> &ret)
  {
    RequirePerFaceMark(m);
    ret.clear();
    int referredBit = FaceType::NewBitFlag();
    tri::UpdateFlags<MeshType>::FaceClear(m,referredBit);

    TriMeshGrid gM;
    gM.Set(m.face.begin(),m.face.end());

    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!(*fi).IsD())
    {
      (*fi).SetUserBit(referredBit);
      Box3< ScalarType> bbox;
      (*fi).GetBBox(bbox);
      std::vector<FaceType*> inBox;
      vcg::tri::GetInBoxFace(m, gM, bbox,inBox);
      bool Intersected=false;
      typename std::vector<FaceType*>::iterator fib;
      for(fib=inBox.begin();fib!=inBox.end();++fib)
      {
        if(!(*fib)->IsUserBit(referredBit) && (*fib != &*fi) )
          if(Clean<MeshType>::TestFaceFaceIntersection(&*fi,*fib)){
            ret.push_back(*fib);
            if(!Intersected) {
              ret.push_back(&*fi);
              Intersected=true;
            }
          }
      }
      inBox.clear();
    }

    FaceType::DeleteBitFlag(referredBit);
    return (ret.size()>0);
  }

  /**
      This function simply test that the vn and fn counters be consistent with the size of the containers and the number of deleted simplexes.
      */
  static bool IsSizeConsistent(MeshType &m)
  {
    int DeletedVertNum=0;
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if((*vi).IsD()) DeletedVertNum++;

    int DeletedEdgeNum=0;
    for (EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
      if((*ei).IsD()) DeletedEdgeNum++;

    int DeletedFaceNum=0;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if((*fi).IsD()) DeletedFaceNum++;

    if(size_t(m.vn+DeletedVertNum) != m.vert.size()) return false;
    if(size_t(m.en+DeletedEdgeNum) != m.edge.size()) return false;
    if(size_t(m.fn+DeletedFaceNum) != m.face.size()) return false;

    return true;
  }

  /**
      This function simply test that all the faces have a consistent face-face topology relation.
      useful for checking that a topology modifying algorithm does not mess something.
      */
  static bool IsFFAdjacencyConsistent(MeshType &m)
  {
    RequireFFAdjacency(m);

    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
      {
        for(int i=0;i<3;++i)
          if(!FFCorrectness(*fi, i)) return false;
      }
    return true;
  }

  /**
      This function simply test that a mesh has some reasonable tex coord.
      */
  static bool HasConsistentPerWedgeTexCoord(MeshType &m)
  {
    tri::RequirePerFaceWedgeTexCoord(m);

    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
      { FaceType &f=(*fi);
        if( ! ( (f.WT(0).N() == f.WT(1).N()) && (f.WT(0).N() == (*fi).WT(2).N()) )  )
          return false; // all the vertices must have the same index.

        if((*fi).WT(0).N() <0) return false; // no undefined texture should be allowed
      }
    return true;
  }

  /**
  Simple check that there are no face with all collapsed tex coords.
  */
  static bool HasZeroTexCoordFace(MeshType &m)
  {
    tri::RequirePerFaceWedgeTexCoord(m);

    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if(!(*fi).IsD())
      {
        if( (*fi).WT(0).P() == (*fi).WT(1).P() && (*fi).WT(0).P() == (*fi).WT(2).P() ) return false;
      }
    return true;
  }


  /**
        This function test if two triangular faces of a mesh intersect.
        It assumes that the faces (as storage) are different (e.g different address)
        If the two faces are different but coincident (same set of vertexes) return true.
        if the faces share an edge no test is done.
        if the faces share only a vertex, the opposite edge is tested against the face
  */
  static	bool TestFaceFaceIntersection(FaceType *f0,FaceType *f1)
  {
    int sv = face::CountSharedVertex(f0,f1);
    if(sv==3) return true;
    if(sv==0) return (vcg::IntersectionTriangleTriangle<FaceType>((*f0),(*f1)));
    //  if the faces share only a vertex, the opposite edge (as a segment) is tested against the face
    //  to avoid degenerate cases where the two triangles have the opposite edge on a common plane
    //  we offset the segment to test toward the shared vertex
    if(sv==1)
    {
      int i0,i1; ScalarType a,b;
      face::FindSharedVertex(f0,f1,i0,i1);
      CoordType shP = f0->V(i0)->P()*0.5;
      if(vcg::IntersectionSegmentTriangle(Segment3<ScalarType>((*f0).V1(i0)->P()*0.5+shP,(*f0).V2(i0)->P()*0.5+shP), *f1, a, b) )
      {
        // a,b are the param coords of the intersection point of the segment.
        if(a+b>=1 || a<=EPSIL || b<=EPSIL ) return false;
        return true;
      }
      if(vcg::IntersectionSegmentTriangle(Segment3<ScalarType>((*f1).V1(i1)->P()*0.5+shP,(*f1).V2(i1)->P()*0.5+shP), *f0, a, b) )
      {
        // a,b are the param coords of the intersection point of the segment.
        if(a+b>=1 || a<=EPSIL || b<=EPSIL ) return false;
        return true;
      }

    }
    return false;
  }



  /**
      This function merge all the vertices that are closer than the given radius
*/
  static int MergeCloseVertex(MeshType &m, const ScalarType radius)
  {
    int mergedCnt=0;
    mergedCnt = ClusterVertex(m,radius);
    RemoveDuplicateVertex(m,true);
    return mergedCnt;
  }

  static int ClusterVertex(MeshType &m, const ScalarType radius)
  {
    if(m.vn==0) return 0;
    // some spatial indexing structure does not work well with deleted vertices...
    tri::Allocator<MeshType>::CompactVertexVector(m);
    typedef vcg::SpatialHashTable<VertexType, ScalarType> SampleSHT;
    SampleSHT sht;
    tri::EmptyTMark<MeshType> markerFunctor;
    std::vector<VertexType*> closests;
    int mergedCnt=0;
    sht.Set(m.vert.begin(), m.vert.end());
    UpdateFlags<MeshType>::VertexClearV(m);
    for(VertexIterator viv = m.vert.begin(); viv!= m.vert.end(); ++viv)
      if(!(*viv).IsD() && !(*viv).IsV())
      {
        (*viv).SetV();
        Point3<ScalarType> p = viv->cP();
        Box3<ScalarType> bb(p-Point3<ScalarType>(radius,radius,radius),p+Point3<ScalarType>(radius,radius,radius));
        GridGetInBox(sht, markerFunctor, bb, closests);
        // qDebug("Vertex %i has %i closest", &*viv - &*m.vert.begin(),closests.size());
        for(size_t i=0; i<closests.size(); ++i)
        {
          ScalarType dist = Distance(p,closests[i]->cP());
          if(dist < radius && !closests[i]->IsV())
          {
            //													printf("%f %f \n",dist,radius);
            mergedCnt++;
            closests[i]->SetV();
            closests[i]->P()=p;
          }
        }
      }
    return mergedCnt;
  }


  static std::pair<int,int>  RemoveSmallConnectedComponentsSize(MeshType &m, int maxCCSize)
  {
    std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
    int TotalCC=ConnectedComponents(m, CCV);
    int DeletedCC=0;

    ConnectedComponentIterator<MeshType> ci;
    for(unsigned int i=0;i<CCV.size();++i)
    {
      std::vector<typename MeshType::FacePointer> FPV;
      if(CCV[i].first<maxCCSize)
      {
        DeletedCC++;
        for(ci.start(m,CCV[i].second);!ci.completed();++ci)
          FPV.push_back(*ci);

        typename std::vector<typename MeshType::FacePointer>::iterator fpvi;
        for(fpvi=FPV.begin(); fpvi!=FPV.end(); ++fpvi)
          Allocator<MeshType>::DeleteFace(m,(**fpvi));
      }
    }
    return std::make_pair(TotalCC,DeletedCC);
  }


  /// Remove the connected components smaller than a given diameter
  // it returns a pair with the number of connected components and the number of deleted ones.
  static std::pair<int,int> RemoveSmallConnectedComponentsDiameter(MeshType &m, ScalarType maxDiameter)
  {
    std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
    int TotalCC=ConnectedComponents(m, CCV);
    int DeletedCC=0;
    tri::ConnectedComponentIterator<MeshType> ci;
    for(unsigned int i=0;i<CCV.size();++i)
    {
      Box3<ScalarType> bb;
      std::vector<typename MeshType::FacePointer> FPV;
      for(ci.start(m,CCV[i].second);!ci.completed();++ci)
      {
        FPV.push_back(*ci);
        bb.Add((*ci)->P(0));
        bb.Add((*ci)->P(1));
        bb.Add((*ci)->P(2));
      }
      if(bb.Diag()<maxDiameter)
      {
        DeletedCC++;
        typename std::vector<typename MeshType::FacePointer>::iterator fpvi;
        for(fpvi=FPV.begin(); fpvi!=FPV.end(); ++fpvi)
          tri::Allocator<MeshType>::DeleteFace(m,(**fpvi));
      }
    }
    return std::make_pair(TotalCC,DeletedCC);
  }

  /// Remove the connected components greater than a given diameter
  // it returns a pair with the number of connected components and the number of deleted ones.
  static std::pair<int,int> RemoveHugeConnectedComponentsDiameter(MeshType &m, ScalarType minDiameter)
  {
    std::vector< std::pair<int, typename MeshType::FacePointer> > CCV;
    int TotalCC=ConnectedComponents(m, CCV);
    int DeletedCC=0;
    tri::ConnectedComponentIterator<MeshType> ci;
    for(unsigned int i=0;i<CCV.size();++i)
    {
      Box3f bb;
      std::vector<typename MeshType::FacePointer> FPV;
      for(ci.start(m,CCV[i].second);!ci.completed();++ci)
      {
        FPV.push_back(*ci);
        bb.Add((*ci)->P(0));
        bb.Add((*ci)->P(1));
        bb.Add((*ci)->P(2));
      }
      if(bb.Diag()>minDiameter)
      {
        DeletedCC++;
        typename std::vector<typename MeshType::FacePointer>::iterator fpvi;
        for(fpvi=FPV.begin(); fpvi!=FPV.end(); ++fpvi)
          tri::Allocator<MeshType>::DeleteFace(m,(**fpvi));
      }
    }
    return std::make_pair(TotalCC,DeletedCC);
  }



  /**
  Select the folded faces using an angle threshold on the face normal.
  The face is selected if the dot product between the face normal and the normal of the plane fitted
  using the vertices of the one ring faces is below the cosThreshold.
  The cosThreshold requires a negative cosine value (a positive value is clamp to zero).
  */
  static void SelectFoldedFaceFromOneRingFaces(MeshType &m, ScalarType cosThreshold)
  {
    tri::RequireVFAdjacency(m);
    tri::RequirePerFaceNormal(m);
    tri::RequirePerVertexNormal(m);
    vcg::tri::UpdateSelection<MeshType>::FaceClear(m);
    vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
    vcg::tri::UpdateNormal<MeshType>::PerVertexNormalized(m);
    vcg::tri::UpdateTopology<MeshType>::VertexFace(m);
    if (cosThreshold > 0)
      cosThreshold = 0;

#pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < m.face.size(); i++)
    {
      std::vector<typename MeshType::VertexPointer> nearVertex;
      std::vector<typename MeshType::CoordType> point;
      typename MeshType::FacePointer f = &m.face[i];
      for (int j = 0; j < 3; j++)
      {
        std::vector<typename MeshType::VertexPointer> temp;
        vcg::face::VVStarVF<typename MeshType::FaceType>(f->V(j), temp);
              typename std::vector<typename MeshType::VertexPointer>::iterator iter = temp.begin();
        for (; iter != temp.end(); iter++)
        {
          if ((*iter) != f->V1(j) && (*iter) != f->V2(j))
          {
            nearVertex.push_back((*iter));
            point.push_back((*iter)->P());
          }
        }
        nearVertex.push_back(f->V(j));
        point.push_back(f->P(j));
      }

      if (point.size() > 3)
      {
        vcg::Plane3<typename MeshType::ScalarType> plane;
        vcg::FitPlaneToPointSet(point, plane);
        float avgDot = 0;
        for (int j = 0; j < nearVertex.size(); j++)
          avgDot += plane.Direction().dot(nearVertex[j]->N());
        avgDot /= nearVertex.size();
        typename MeshType::VertexType::NormalType normal;
        if (avgDot < 0)
          normal = -plane.Direction();
        else
          normal = plane.Direction();
        if (normal.dot(f->N()) < cosThreshold)
          f->SetS();
      }
    }
  }

}; // end class
/*@}*/

} //End Namespace Tri
} // End Namespace vcg
#endif
