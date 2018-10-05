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
#ifndef __VCG_TRI_UPDATE_SELECTION
#define __VCG_TRI_UPDATE_SELECTION

#include <queue>
#include <vcg/complex/algorithms/update/flag.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh
/// \brief A stack for saving and restoring selection.
/**
  This class is used to save the current selection onto a stack for later use.
  \todo it should be generalized to other attributes with a templated approach.
*/
template <class ComputeMeshType>
class SelectionStack
{
  typedef typename ComputeMeshType::template PerVertexAttributeHandle< bool > vsHandle;
  typedef typename ComputeMeshType::template PerEdgeAttributeHandle< bool >   esHandle;
  typedef typename ComputeMeshType::template PerFaceAttributeHandle< bool >   fsHandle;

public:
  SelectionStack(ComputeMeshType &m)
  {
    _m=&m;
  }

  bool push()
  {
    vsHandle vsH = Allocator<ComputeMeshType>::template AddPerVertexAttribute< bool >(*_m);
    esHandle esH = Allocator<ComputeMeshType>::template AddPerEdgeAttribute< bool >(*_m);
    fsHandle fsH = Allocator<ComputeMeshType>::template AddPerFaceAttribute< bool >  (*_m);
    typename ComputeMeshType::VertexIterator vi;
    for(vi = _m->vert.begin(); vi != _m->vert.end(); ++vi)
      if( !(*vi).IsD() ) vsH[*vi] = (*vi).IsS() ;

    typename ComputeMeshType::EdgeIterator ei;
    for(ei = _m->edge.begin(); ei != _m->edge.end(); ++ei)
      if( !(*ei).IsD() ) esH[*ei] = (*ei).IsS() ;

    typename ComputeMeshType::FaceIterator fi;
    for(fi = _m->face.begin(); fi != _m->face.end(); ++fi)
      if( !(*fi).IsD() ) fsH[*fi] = (*fi).IsS() ;

    vsV.push_back(vsH);
    esV.push_back(esH);
    fsV.push_back(fsH);
    return true;
  }

  bool popOr()
  {
    return pop(true);
  }

  bool pop(bool mergeFlag=false)
  {
    if(vsV.empty()) return false;
    vsHandle vsH = vsV.back();
    esHandle esH = esV.back();
    fsHandle fsH = fsV.back();
    if(! (Allocator<ComputeMeshType>::template IsValidHandle(*_m, vsH))) return false;

    typename ComputeMeshType::VertexIterator vi;
    for(vi = _m->vert.begin(); vi != _m->vert.end(); ++vi)
      if( !(*vi).IsD() )
      {
        if(vsH[*vi]) 
          (*vi).SetS();
        else
          if(!mergeFlag)
            (*vi).ClearS();
      }

    typename ComputeMeshType::EdgeIterator ei;
    for(ei = _m->edge.begin(); ei != _m->edge.end(); ++ei)
      if( !(*ei).IsD() )
      {
        if(esH[*ei]) 
          (*ei).SetS();
        else
          if(!mergeFlag)
            (*ei).ClearS();
      }
    typename ComputeMeshType::FaceIterator fi;
    for(fi = _m->face.begin(); fi != _m->face.end(); ++fi)
      if( !(*fi).IsD() )
      {  
        if(fsH[*fi]) 
          (*fi).SetS();
        else
          if(!mergeFlag)
            (*fi).ClearS();
      }

    Allocator<ComputeMeshType>::template DeletePerVertexAttribute<bool>(*_m,vsH);
    Allocator<ComputeMeshType>::template DeletePerEdgeAttribute<bool>(*_m,esH);
    Allocator<ComputeMeshType>::template DeletePerFaceAttribute<bool>(*_m,fsH);
    vsV.pop_back();
    esV.pop_back();
    fsV.pop_back();
    return true;
  }

private:
  ComputeMeshType *_m;
  std::vector<vsHandle> vsV;
  std::vector<esHandle> esV;
  std::vector<fsHandle> fsV;
};

/// \ingroup trimesh

/// \headerfile selection.h vcg/complex/algorithms/update/selection.h

/// \brief Management, updating and conditional computation of selections  (per-vertex, per-edge, and per-face).
/**
This class is used to compute or update the selected bit flag that can be stored in the vertex, edge or face component of a mesh. 
*/

template <class ComputeMeshType>
class UpdateSelection
{

public:
typedef ComputeMeshType MeshType;
typedef	typename MeshType::ScalarType			ScalarType;
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::EdgeIterator   EdgeIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;
typedef typename vcg::Box3<ScalarType>  Box3Type;

/// \brief This function select all the vertices.
static size_t VertexAll(MeshType &m)
{
  for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    if( !(*vi).IsD() )	(*vi).SetS();
  return m.vn;
}

/// \brief This function select all the edges.
static size_t EdgeAll(MeshType &m)
{
  for(EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
    if( !(*ei).IsD() )	(*ei).SetS();
  return m.fn;
}
/// \brief This function select all the faces.
static size_t FaceAll(MeshType &m)
{
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() )	(*fi).SetS();
  return m.fn;
}

/// \brief This function clear the selection flag for all the vertices.
static size_t VertexClear(MeshType &m)
{
  for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    if( !(*vi).IsD() )	(*vi).ClearS();
  return 0;
}

/// \brief This function clears the selection flag for all the edges.
static size_t EdgeClear(MeshType &m)
{
  for(EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
    if( !(*ei).IsD() )	(*ei).ClearS();
  return 0;
}

/// \brief This function clears the selection flag for all the faces.
static size_t FaceClear(MeshType &m)
{
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() )	(*fi).ClearS();
  return 0;
}

/// \brief This function clears the selection flag for all the elements of a mesh (vertices, edges, and faces).
static void Clear(MeshType &m)
{
  VertexClear(m);
  EdgeClear(m);
  FaceClear(m);
}

/// \brief This function returns the number of selected faces.
static size_t FaceCount(MeshType &m)
{
  size_t selCnt=0;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    if(!(*fi).IsD() && (*fi).IsS()) ++selCnt;
  return selCnt;
}

/// \brief This function returns the number of selected edges.
static size_t EdgeCount(MeshType &m)
{
  size_t selCnt=0;
  for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
    if(!(*ei).IsD() && (*ei).IsS()) ++selCnt;
  return selCnt;
}

/// \brief This function returns the number of selected vertices.
static size_t VertexCount(MeshType &m)
{
  size_t selCnt=0;
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
    if(!(*vi).IsD() && (*vi).IsS()) ++selCnt;
  return selCnt;
}

/// \brief This function inverts the selection flag for all the faces.
static size_t FaceInvert(MeshType &m)
{
  size_t selCnt=0;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if(!(*fi).IsD())
      {
        if((*fi).IsS()) (*fi).ClearS();
        else {
          (*fi).SetS();
          ++selCnt;
        }
      }
  return selCnt;
}

/// \brief This function inverts the selection flag for all the edges.
static size_t EdgeInvert(MeshType &m)
{
  size_t selCnt=0;
  for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
      if(!(*ei).IsD())
      {
        if((*ei).IsS()) (*ei).ClearS();
        else {
          (*ei).SetS();
          ++selCnt;
        }
      }
  return selCnt;
}

/// \brief This function inverts the selection flag for all the vertices.
static size_t VertexInvert(MeshType &m)
{
  size_t selCnt=0;
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
      if(!(*vi).IsD())
      {
        if((*vi).IsS()) (*vi).ClearS();
        else {
          (*vi).SetS();
          ++selCnt;
        }
      }
  return selCnt;
}

/// \brief Select all the vertices that are touched by at least a single selected faces
static size_t VertexFromFaceLoose(MeshType &m, bool preserveSelection=false)
{
  size_t selCnt=0;

  if(!preserveSelection) VertexClear(m);
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() && (*fi).IsS())
      for(int i = 0; i < (*fi).VN(); ++i)
        if( !(*fi).V(i)->IsS()) { (*fi).V(i)->SetS(); ++selCnt; }
  return selCnt;
}

/// \brief Select all the vertices that are touched by at least a single selected edge
static size_t VertexFromEdgeLoose(MeshType &m, bool preserveSelection=false)
{
  size_t selCnt=0;

  if(!preserveSelection) VertexClear(m);
  for(EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
    if( !(*ei).IsD() && (*ei).IsS())
    {
      if( !(*ei).V(0)->IsS()) { (*ei).V(0)->SetS(); ++selCnt; }
      if( !(*ei).V(1)->IsS()) { (*ei).V(1)->SetS(); ++selCnt; }
    }
  return selCnt;
}

/// \brief Select ONLY the vertices that are touched ONLY by selected faces
/** In other words this function will select all the vertices having all the faces incident on them selected.
*/
static size_t VertexFromFaceStrict(MeshType &m, bool preserveSelection=false)
{
  SelectionStack<MeshType> ss(m);
  if(preserveSelection) ss.push();  
  VertexFromFaceLoose(m);
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() && !(*fi).IsS())
      for(int i = 0; i < (*fi).VN(); ++i)
        (*fi).V(i)->ClearS();

  if(preserveSelection) ss.popOr();
  return VertexCount(m);
}

/// \brief Select ONLY the faces with ALL the vertices selected
static size_t FaceFromVertexStrict(MeshType &m, bool preserveSelection=false)
{
  SelectionStack<MeshType> ss(m);
  if(preserveSelection) ss.push();
  size_t selCnt=0;
  FaceClear(m);
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD())
    {
      bool selFlag=true;
      for(int i = 0; i < (*fi).VN(); ++i)
        if(!(*fi).V(i)->IsS())
          selFlag =false;
      if(selFlag)
      {
        (*fi).SetS();
        ++selCnt;
      }
    }
  
  if(preserveSelection) ss.popOr();
  return selCnt;
}

/// \brief Select all the faces with at least one selected vertex
static size_t FaceFromVertexLoose(MeshType &m, bool preserveSelection=false)
{
  size_t selCnt=0;
  if(!preserveSelection) FaceClear(m);
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD())
    {
      bool selVert=false;
      for(int i = 0; i < (*fi).VN(); ++i)
        if((*fi).V(i)->IsS()) 
            selVert=true;

      if(selVert) {
        (*fi).SetS();
        ++selCnt;
      }
    }
  return selCnt;
}

/// \brief This function select the vertices with the border flag set
static size_t VertexFromBorderFlag(MeshType &m, bool preserveSelection=false)
{
  size_t selCnt=0;
  if(!preserveSelection) VertexClear(m);
  for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    if( !(*vi).IsD() )
    {
      if((*vi).IsB() )
      {
        (*vi).SetS();
        ++selCnt;
      }
    }
  return selCnt;
}

/// \brief This function select the faces that have an edge with the border flag set.
static size_t FaceFromBorderFlag(MeshType &m, bool preserveSelection=false)
{
  tri::RequireTriangularMesh(m);
  size_t selCnt=0;
  if(!preserveSelection) FaceClear(m);
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() )
    {
      bool bordFlag=false;
      for(int i = 0; i < 3; ++i)
        if((*fi).IsB(i)) bordFlag=true;
      if(bordFlag)
      {
        (*fi).SetS();
        ++selCnt;
      }
    }
  return selCnt;
}

/// \brief This function select the faces that have an edge outside the given range.
/// You can skip the second parameter to choose all the edges smaller than a given lenght
static size_t FaceOutOfRangeEdge(MeshType &m, ScalarType MinEdgeThr, ScalarType MaxEdgeThr=(std::numeric_limits<ScalarType>::max)(), bool preserveSelection=false)
{
  if(!preserveSelection) FaceClear(m);
  size_t selCnt = 0;
  MinEdgeThr=MinEdgeThr*MinEdgeThr;
  MaxEdgeThr=MaxEdgeThr*MaxEdgeThr;
  for(FaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
    if(!(*fi).IsD())
      {
        for(int i=0;i<(*fi).VN();++i)
        {
          const ScalarType squaredEdge=SquaredDistance((*fi).V0(i)->cP(),(*fi).V1(i)->cP());
          if((squaredEdge<=MinEdgeThr) || (squaredEdge>=MaxEdgeThr) )
          {
            selCnt++;
            (*fi).SetS();
            break; // skip the rest of the edges of the tri
          }
        }
      }
      return selCnt;
}

/// \brief This function expand current selection to cover the whole connected component.
static size_t FaceConnectedFF(MeshType &m)
{
  // it also assumes that the FF adjacency is well computed.
  RequireFFAdjacency(m);
  UpdateFlags<MeshType>::FaceClearV(m);
  
  std::deque<FacePointer> visitStack;
  size_t selCnt=0;
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if( !(*fi).IsD() && (*fi).IsS() && !(*fi).IsV() )
      visitStack.push_back(&*fi);
  
  while(!visitStack.empty())
  {
    FacePointer fp = visitStack.front();
    visitStack.pop_front();
    assert(!fp->IsV());
    fp->SetV();
    for(int i=0;i<fp->VN();++i) {
      FacePointer ff = fp->FFp(i);
      if(! ff->IsS())
      {
        ff->SetS();
        ++selCnt;
        visitStack.push_back(ff);
        assert(!ff->IsV());
      }
    }
  }
  return selCnt;
}
/// \brief Select the faces whose quality is in the specified closed interval.
static size_t FaceFromQualityRange(MeshType &m,float minq, float maxq, bool preserveSelection=false)
{
  size_t selCnt=0;
  if(!preserveSelection) FaceClear(m);
  RequirePerFaceQuality(m);
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if(!(*fi).IsD())
      {
        if( (*fi).Q()>=minq &&  (*fi).Q()<=maxq )
          {
            (*fi).SetS();
            ++selCnt;
          }
      }
  return selCnt;
}

/// \brief Select the vertices whose quality is in the specified closed interval.
static size_t VertexFromQualityRange(MeshType &m,float minq, float maxq, bool preserveSelection=false)
{
  size_t selCnt=0;
  if(!preserveSelection) VertexClear(m);
  RequirePerVertexQuality(m);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
      if(!(*vi).IsD())
      {
        if( (*vi).Q()>=minq &&  (*vi).Q()<=maxq )
                    {
                        (*vi).SetS();
                        ++selCnt;
                    }
      }
  return selCnt;
}

/// \brief Select the vertices contained in the specified Box
static int VertexInBox( MeshType & m, const Box3Type &bb, bool preserveSelection=false)
{
  if(!preserveSelection) VertexClear(m);
  int selCnt=0;
  for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) if(!(*vi).IsD())
  {
    if(bb.IsIn((*vi).cP()) ) {
      (*vi).SetS();
      ++selCnt;
    }
  }
  return selCnt;
}


void VertexNonManifoldEdges(MeshType &m, bool preserveSelection=false)
{
  assert(HasFFTopology(m));

  if(!preserveSelection) VertexClear(m);
  for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)	if (!fi->IsD())
    {
      for(int i=0;i<fi->VN();++i)
      if(!IsManifold(*fi,i)){
        (*fi).V0(i)->SetS();
        (*fi).V1(i)->SetS();
        }
    }
}

}; // end class

}	// End namespace
}	// End namespace


#endif
