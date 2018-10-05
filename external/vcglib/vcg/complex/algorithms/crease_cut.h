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

#ifndef __VCG_CREASE_CUT
#define __VCG_CREASE_CUT
#include<vcg/simplex/face/jumping_pos.h>
#include<vcg/complex/algorithms/update/normal.h>
namespace vcg {
namespace tri {

/** \brief Open a mesh cutting all the edges where the two faces make an angle *larger* than the indicated threshold
 */

template<class MESH_TYPE>
void CreaseCut(MESH_TYPE &m, float angleRad)
{
  tri::UpdateFlags<MESH_TYPE>::FaceFauxSignedCrease(m, -angleRad, angleRad);
  CutMeshAlongNonFauxEdges(m);
}

/**
 * \brief Open a mesh along non-faux edges
 * 
 * Duplicate exisiting vertices so that non-faux edges become boundary edges. 
 * It assume FF topology and manifoldness.
 * The idea is that we scan faces around each vertex duplicating it each time we encounter a marked edge. 
 * 
 */
template<class MESH_TYPE>
void CutMeshAlongNonFauxEdges(MESH_TYPE &m)
{
  typedef typename MESH_TYPE::FaceIterator		FaceIterator;
  typedef typename MESH_TYPE::FaceType				FaceType;
  
  tri::Allocator<MESH_TYPE>::CompactVertexVector(m);
  tri::Allocator<MESH_TYPE>::CompactFaceVector(m);    
  tri::RequireFFAdjacency(m);
  
  tri::UpdateFlags<MESH_TYPE>::VertexClearV(m);
  std::vector<int> indVec(m.fn*3,-1);
  int newVertexCounter=m.vn;
  int startVn=m.vn;
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    for(int j=0;j<3;++j)
      if(!(*fi).V(j)->IsV() )  // foreach unvisited vertex we loop around it searching for creases.
      {
        (*fi).V(j)->SetV();
        
        face::JumpingPos<FaceType> iPos(&*fi,j,(*fi).V(j));
        size_t vertInd = Index(m, iPos.V()); 
        bool isBorderVertex = iPos.FindBorder();   // for border vertex we start from the border.
        face::JumpingPos<FaceType> startPos=iPos;
        if(!isBorderVertex)                        // for internal vertex we search the first crease and start from it
        {
          do {
              bool creaseFlag = !iPos.IsFaux();
              iPos.NextFE();
              if(creaseFlag) break;
          } while (startPos!=iPos);
          startPos=iPos;                       // the found crease become the new starting pos.
        }
        
        int locCreaseCounter=0;
        int curVertexCounter =vertInd;
        
        do { // The real Loop          
          size_t faceInd = Index(m,iPos.F());
          indVec[faceInd*3+ iPos.VInd()] = curVertexCounter;
          
          if(!iPos.IsFaux()) 
          { //qDebug("  Crease FOUND");
            ++locCreaseCounter;
            curVertexCounter=newVertexCounter;
            newVertexCounter++;
          }
          iPos.NextFE();
        } while (startPos!=iPos);
        if(locCreaseCounter>0 && (!isBorderVertex) ) newVertexCounter--;
        //printf("For vertex %i found %i creases\n",vertInd,locCreaseCounter);
      }
  } // end foreach face/vert 

  // Now the indVec vector contains for each the new index of each vertex (duplicated as necessary)
  // We do a second loop to copy split vertexes into new positions
  tri::Allocator<MESH_TYPE>::AddVertices(m,newVertexCounter-m.vn);
  
  tri::UpdateFlags<MESH_TYPE>::VertexClearV(m);
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    for(int j=0;j<3;++j) 
    {
      size_t faceInd = Index(m, *fi);
      size_t vertInd = Index(m, (*fi).V(j));
      int curVertexInd = indVec[faceInd*3+ j];
      assert(curVertexInd != -1);
      assert(curVertexInd < m.vn);
      if(curVertexInd < startVn) assert(size_t(curVertexInd) == vertInd);
      if(curVertexInd >= startVn)
      {
        m.vert[curVertexInd].ImportData(*((*fi).V(j)));
        (*fi).V(j) = & m.vert[curVertexInd];
      }
    }
}

} // end namespace tri
} // end namespace vcg
#endif

