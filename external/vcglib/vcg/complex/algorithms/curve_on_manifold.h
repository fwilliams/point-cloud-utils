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
#ifndef __VCGLIB_CURVE_ON_SURF_H
#define __VCGLIB_CURVE_ON_SURF_H

#include<vcg/complex/complex.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/math/histogram.h>
#include<vcg/space/distance3.h>
#include<eigenlib/Eigen/Core>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <wrap/io_trimesh/export_ply.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh
/// \brief A class for managing curves on a 2manifold.
/**
  This class is used to project/simplify/smooth polylines over a triangulated surface. 
  
*/

template <class MeshType>
class CoM
{
public:
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType     CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef Box3  <ScalarType>               Box3Type;
  typedef typename vcg::GridStaticPtr<FaceType, ScalarType> MeshGrid;  
  typedef typename vcg::GridStaticPtr<EdgeType, ScalarType> EdgeGrid;
  typedef typename face::Pos<FaceType> PosType;
  
  class Param 
  {
  public:
    
    ScalarType surfDistThr; // Distance between surface and curve; used in simplify and refine
    ScalarType polyDistThr; // Distance between the 
    ScalarType minRefEdgeLen;  // Minimal admitted Edge Lenght (used in refine: never make edge shorther than this value) 
    ScalarType maxSimpEdgeLen; // Minimal admitted Edge Lenght (used in simplify: never make edges longer than this value) 
    ScalarType maxSmoothDelta; // The maximum movement that is admitted during smoothing.
    ScalarType maxSnapThr;     // The maximum distance allowed when snapping a vertex of the polyline onto a mesh vertex
    ScalarType gridBailout;    // The maximum distance bailout used in grid sampling
    
    Param(MeshType &m) { Default(m);}
    
    void Default(MeshType &m)
    {
      surfDistThr   = m.bbox.Diag()/1000.0;
      polyDistThr   = m.bbox.Diag()/5000.0;
      minRefEdgeLen    = m.bbox.Diag()/16000.0;
      maxSimpEdgeLen    = m.bbox.Diag()/1000.0;
      maxSmoothDelta =   m.bbox.Diag()/100.0;
      maxSnapThr =       m.bbox.Diag()/10000.0;
      gridBailout =      m.bbox.Diag()/20.0;
    }
    void Dump() const
    {
      printf("surfDistThr    = %6.4f\n",surfDistThr   );
      printf("polyDistThr    = %6.4f\n",polyDistThr   );
      printf("minRefEdgeLen  = %6.4f\n",minRefEdgeLen    );
      printf("maxSimpEdgeLen = %6.4f\n",maxSimpEdgeLen    );
      printf("maxSmoothDelta = %6.4f\n",maxSmoothDelta);
    }
  };
  
  
  
  // The Data Members
  
  MeshType &base; 
  MeshGrid uniformGrid;
  
  Param par; 
  CoM(MeshType &_m) :base(_m),par(_m){}
 
  
  //******** CUT TREE GENERATION 

    
  // Given two points return true if on the base mesh there exist an edge with that two coords
  // if return true the pos indicate the found edge. 
  bool ExistEdge(KdTree<ScalarType> &kdtree, CoordType &p0, CoordType &p1, PosType &fpos)
  {
    ScalarType locEps = SquaredDistance(p0,p1)/100000.0;
    
    VertexType *v0=0,*v1=0;
    unsigned int veInd;
    ScalarType sqdist;
    kdtree.doQueryClosest(p0,veInd,sqdist);
    if(sqdist<locEps) 
      v0 = &base.vert[veInd];
    kdtree.doQueryClosest(p1,veInd,sqdist);
    if(sqdist<locEps) 
      v1 = &base.vert[veInd];
    if(v0 && v1)
    {
      face::VFIterator<FaceType> vfi(v0);
      while(!vfi.End())
      { 
        if(vfi.V1()==v1)
        {
          fpos = PosType(vfi.F(),vfi.I(), v0);
          return true;
        }
        if(vfi.V2()==v1)
        {
          fpos = PosType(vfi.F(),(vfi.I()+1)%3, v1);
          return true;
        }
        ++vfi;
      }      
    }
    return false;
  }
  
  bool OptimizeTree(MeshType &t)
  {
    tri::Allocator<MeshType>::CompactEveryVector(t);
    tri::UpdateTopology<MeshType>::VertexEdge(t);
    tri::UpdateTopology<MeshType>::VertexFace(base);
    VertexConstDataWrapper<MeshType > vdw(base);
    KdTree<ScalarType> kdtree(vdw);
    
    // First simple loop that search for 2->1 moves. 
    for(VertexIterator vi=t.vert.begin();vi!=t.vert.end();++vi)
    {
      std::vector<VertexType *> starVec;
      edge::VVStarVE(&*vi,starVec);
      if(starVec.size()==2)
      {
        PosType pos;
        if(ExistEdge(kdtree,starVec[0]->P(),starVec[1]->P(),pos))
          edge::VEEdgeCollapse(t,&*vi);
      }
    }
    return (t.en < t.edge.size());
  }
  
  int findNonVisitedEdgesDuringRetract(VertexType * vp, EdgeType * &ep)
  {
    std::vector<EdgeType *> starVec;
    edge::VEStarVE(&*vp,starVec);
    int cnt =0;
    for(size_t i=0;i<starVec.size();++i) {
      if(!starVec[i]->IsV()) {
        cnt++;
        ep = starVec[i];
      }
    }  
    return cnt;
  }
  
  void Retract(MeshType &t)
  {
    tri::Clean<MeshType>::RemoveDuplicateVertex(t);
    printf("Retracting a tree of %i edges and %i vertices\n",t.en,t.vn);
    tri::UpdateTopology<MeshType>::VertexEdge(t);
  
    std::stack<VertexType *> vertStack;
  
    // Put on the stack all the vertex with just a single incident edge. 
    for(VertexIterator vi=t.vert.begin();vi!=t.vert.end();++vi)
    {
      std::vector<EdgeType *> starVec;
      edge::VEStarVE(&*vi,starVec);
      if(starVec.size()==1)
        vertStack.push(&*vi);
    }
  
    tri::UpdateFlags<MeshType>::EdgeClearV(t);
    tri::UpdateFlags<MeshType>::VertexClearV(t);
    
    int unvisitedEdgeNum = t.en;
    while((!vertStack.empty()) && (unvisitedEdgeNum > 2) )
    {
      VertexType *vp = vertStack.top();
      vertStack.pop();
      vp->C()=Color4b::Blue;
      EdgeType *ep=0;
      int eCnt =  findNonVisitedEdgesDuringRetract(vp,ep);
      if(eCnt==1) // We have only one non visited edge over vp
      {
        assert(!ep->IsV());
        ep->SetV();
        --unvisitedEdgeNum;
        VertexType *otherVertP;
        if(ep->V(0)==vp) otherVertP = ep->V(1);
        else otherVertP = ep->V(0);
        vertStack.push(otherVertP);
      }
    }
    assert(unvisitedEdgeNum >0);
    
    for(size_t i =0; i<t.edge.size();++i)
      if(t.edge[i].IsV()) tri::Allocator<MeshType>::DeleteEdge(t,t.edge[i]);
    assert(t.en >0);
    tri::Clean<MeshType>::RemoveUnreferencedVertex(t);
    tri::Allocator<MeshType>::CompactEveryVector(t);
  }
  
  void CleanSpuriousDanglingEdges(MeshType &poly)
  {    
    EdgeGrid edgeGrid;
    edgeGrid.Set(poly.edge.begin(), poly.edge.end());    
    std::vector< std::pair<VertexType *, VertexType *> > mergeVec;
    UpdateFlags<MeshType>::FaceClearV(base);
    UpdateTopology<MeshType>::FaceFace(base);
    for(int i=0;i<base.fn;++i)
    {
      FaceType *fp=&base.face[i];
      if(!fp->IsV())
      {
        for(int j=0;j<3;++j)
          if(face::IsBorder(*fp,j))
          {
            face::Pos<FaceType> startPos(fp,int(j));
            assert(startPos.IsBorder());
            face::Pos<FaceType> curPos=startPos;
            face::Pos<FaceType> prevPos=startPos;
            int edgeCnt=0;
            do
            {
              prevPos=curPos;
              curPos.F()->SetV();
              curPos.NextB();
              edgeCnt++;
              if(Distance(curPos.V()->P(),prevPos.VFlip()->P())<0.000000001f)
              {
                Point3f endp = curPos.VFlip()->P();
                float minDist=par.gridBailout;
                Point3f closestP;
                EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,endp,par.gridBailout,minDist,closestP);        
                if(minDist > par.polyDistThr){
                  mergeVec.push_back(std::make_pair(curPos.V(),prevPos.VFlip()));
                  printf("Vertex %i and %i should be merged\n",tri::Index(base,curPos.V()),tri::Index(base,prevPos.VFlip()));
                }
              }
              
            } while(curPos!=startPos);
            printf("walked along a border of %i edges\n",edgeCnt);
            break;
          }
      }
    }
    printf("Found %i vertex pairs to be merged\n",mergeVec.size());
    for(int k=0;k<mergeVec.size();++k)
    {
      VertexType *vdel=mergeVec[k].first;
      VertexType *vmerge=mergeVec[k].second;
      for(int i=0;i<base.fn;++i)
      {
        FaceType *fp=&base.face[i];
        for(int j=0;j<3;++j)
          if(fp->V(j)==vdel) fp->V(j)=vmerge;        
      }
      Allocator<MeshType>::DeleteVertex(base,*vdel);
    }
    Allocator<MeshType>::CompactEveryVector(base);
    
   
      for(int i=0;i<base.fn;++i)
      {
        FaceType *fp=&base.face[i];
        for(int j=0;j<3;++j)
          if(fp->V0(j)==fp->V1(j)) 
        {
            Allocator<MeshType>::DeleteFace(base,*fp);
            printf("Deleted face %i\n",tri::Index(base,fp));
        }
        
      }
      Allocator<MeshType>::CompactEveryVector(base);
  }  
  
  
  // \brief This function build a cut tree. 
  //
  // First we make a bread first FF face visit. 
  // Each time that we encounter a visited face we cut add to the tree the edge 
  // that brings to the already visited face.
  // this structure build a dense graph and we retract this graph retracting each 
  // leaf until we remains with just the loops that cuts the object. 
  
  void BuildVisitTree(MeshType &dualMesh)
  {
    tri::UpdateTopology<MeshType>::FaceFace(base);
    tri::UpdateFlags<MeshType>::FaceClearV(base);
    
    std::vector<face::Pos<FaceType> > visitStack; // the stack contain the pos on the 'starting' face. 
    
    
    MeshType treeMesh; // this mesh will contain the set of the non traversed edge 
    
    base.face[0].SetV();
    for(int i=0;i<3;++i)
      visitStack.push_back(PosType(&(base.face[0]),i,base.face[0].V(i)));
  
    int cnt=1;
    
    while(!visitStack.empty())
    {
      std::swap(visitStack.back(),visitStack[rand()%visitStack.size()]);
      PosType c = visitStack.back();
      visitStack.pop_back();
      assert(c.F()->IsV());
      c.F()->C() = Color4b::ColorRamp(0,base.fn,cnt);
      c.FlipF();
      if(!c.F()->IsV())
      {
        tri::Allocator<MeshType>::AddEdge(treeMesh,Barycenter(*(c.FFlip())),Barycenter(*(c.F())));
        ++cnt;
        c.F()->SetV();
        c.FlipE();c.FlipV();
        visitStack.push_back(c);
        c.FlipE();c.FlipV();
        visitStack.push_back(c);
      }
      else
      {
        tri::Allocator<MeshType>::AddEdge(dualMesh,c.V()->P(),c.VFlip()->P());
      }
    }
    assert(cnt==base.fn);
    
    Retract(dualMesh);
} 
  
  /*
   * Given a mesh <m> and a polyline <e> whose edges are a subset of edges of m
   * This function marks the edges of m as non-faux when they coincide with the polyline ones. 
   * 
   * Use this function toghether with the CutMeshAlongCrease function to actually cut the mesh. 
   */
  
  void MarkFauxEdgeWithPolyLine(MeshType &m, MeshType &e)
  {
    tri::UpdateFlags<MeshType>::FaceSetF(m);
    tri::UpdateTopology<MeshType>::VertexEdge(e);
    VertexConstDataWrapper<MeshType > vdw(e);
    KdTree<ScalarType> edgeTree(vdw);
    
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      for(int i=0;i<3;++i)
      {
        ScalarType locEps = SquaredDistance(fi->P0(i), fi->P1(i))/100000.0;
        unsigned int veInd;
        ScalarType sqdist;
        edgeTree.doQueryClosest(fi->P(i),veInd,sqdist);
        if(sqdist<locEps)
        {
          std::vector<VertexPointer> starVecVp;
          edge::VVStarVE(&(e.vert[veInd]),starVecVp);      
          for(size_t j=0;j<starVecVp.size();++j)
          {
            if(SquaredDistance(starVecVp[j]->P(),fi->P1(i))< locEps)
              fi->ClearF(i);
          }
        }        
      }
    }
    
  }
  
  
  
  
  void SnapPolyline(MeshType &t)
  {
    printf("Selected vert num %i\n",tri::UpdateSelection<MeshType>::VertexCount(t));
    
    tri::UpdateFlags<MeshType>::VertexClearS(base);
    tri::UpdateFlags<MeshType>::VertexClearS(t);
    tri::Allocator<MeshType>::CompactEveryVector(t);
    VertexConstDataWrapper<MeshType > vdw(base);
    KdTree<ScalarType> kdtree(vdw);

    for(VertexIterator vi=t.vert.begin();vi!=t.vert.end();++vi)
    {
      unsigned int veInd;
      ScalarType sqdist;
      kdtree.doQueryClosest(vi->P(),veInd,sqdist);
      VertexPointer vp = &base.vert[veInd];
      if(!vp->IsS())
      {
        ScalarType dist = sqrt(sqdist);
        std::vector<VertexPointer> starVecVp;        
        face::VVStarVF<FaceType>(vp,starVecVp);
        ScalarType minEdgeLen = std::numeric_limits<ScalarType>::max();
        for(int i=0;i<starVecVp.size();++i)
          minEdgeLen = std::min(Distance(vp->P(),starVecVp[i]->P()),minEdgeLen);
        if(dist < minEdgeLen/3.0)
        {
          vi->P() = vp->P();
          vi->SetS();
          vp->SetS();
        }
      }
    }
    printf("Selected vert num %i\n",tri::UpdateSelection<MeshType>::VertexCount(t));    
  }
  
  
  
  
  float MinDistOnEdge(Point3f samplePnt, EdgeGrid &edgeGrid, MeshType &poly, Point3f &closestPoint)
  {
      float polyDist;
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,par.gridBailout,polyDist,closestPoint);        
      return polyDist;    
  }
  
  // Given an edge of a mesh, supposedly intersecting the polyline, 
  // we search on it the closest point to the polyline
  static float MinDistOnEdge(VertexType *v0,VertexType *v1, EdgeGrid &edgeGrid, MeshType &poly, Point3f &closestPoint)
  {
    float minPolyDist = std::numeric_limits<ScalarType>::max();
    const float sampleNum = 50;
    const float maxDist = poly.bbox.Diag()/10.0;
    for(float k = 0;k<sampleNum+1;++k)
    {
      float polyDist;
      Point3f closestPPoly;
      Point3f samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,maxDist,polyDist,closestPPoly);        
      
      if(polyDist < minPolyDist)
      {
        minPolyDist = polyDist;
        closestPoint = samplePnt;
//        closestPoint = closestPPoly;
      }
    }
    return minPolyDist;    
  }
  
  
  /**
   * @brief ExtractVertex
   * must extract an unambiguous representation of a vertex 
   * to be used with attribute_seam.h
   * 
   */
  static inline void ExtractVertex(const MeshType & srcMesh, const FaceType & f, int whichWedge, const MeshType & dstMesh, VertexType & v)
  {
      (void)srcMesh;
      (void)dstMesh;
      // This is done to preserve every single perVertex property
      // perVextex Texture Coordinate is instead obtained from perWedge one.
      v.ImportData(*f.cV(whichWedge));
      v.C() = f.cC();
  }
  
  static inline bool CompareVertex(const MeshType & m, const VertexType & vA, const VertexType & vB)
  {
      (void)m;
      
      if(vA.C() == Color4b(Color4b::Red) && vB.C() == Color4b(Color4b::Blue) ) return false;
      if(vA.C() == Color4b(Color4b::Blue) && vB.C() == Color4b(Color4b::Red) ) return false;
      return true;      
  }
  
  
  
  static Point3f QLerp(VertexType *v0, VertexType *v1)
  {
    
    float qSum = fabs(v0->Q())+fabs(v1->Q());      
    float w0 = (qSum - fabs(v0->Q()))/qSum;
    float w1 = (qSum - fabs(v1->Q()))/qSum;      
    return v0->P()*w0 + v1->P()*w1;      
  }
  
  class QualitySign
  {
  public:
      EdgeGrid &edgeGrid;
      MeshType &poly;
      CoM &com;
      QualitySign(EdgeGrid &_e,MeshType &_poly, CoM &_com):edgeGrid(_e),poly(_poly),com(_com) {};
      bool operator()(face::Pos<FaceType> ep) const
      {
        VertexType *v0 = ep.V();
        VertexType *v1 = ep.VFlip();
        if(v0->Q() * v1->Q() < 0)
        { 
          Point3f pp = QLerp(v0,v1);
          Point3f closestP;
          if(com.MinDistOnEdge(pp,edgeGrid,poly,closestP)<com.par.polyDistThr) return true;
          float minDist = com.MinDistOnEdge(v0,v1,edgeGrid,poly,closestP);
          if(minDist < com.par.polyDistThr) return true;
        }
        return false;
      }
  };
  
  struct QualitySignSplit : public std::unary_function<face::Pos<FaceType> ,  Point3f>
  {
    EdgeGrid &edgeGrid;
    MeshType &poly;
    CoM &com;
    vector<int> &newVertVec; 
    
    QualitySignSplit(EdgeGrid &_e,MeshType &_p, CoM &_com, vector<int> &_vec):edgeGrid(_e),poly(_p),com(_com),newVertVec(_vec) {};
      void operator()(VertexType &nv, face::Pos<FaceType> ep)
      {
        VertexType *v0 = ep.V();
        VertexType *v1 = ep.VFlip();
        Point3f pp = QLerp(v0,v1);
        Point3f closestP;
        com.MinDistOnEdge(pp,edgeGrid,poly,closestP);
        
//        float minDist = MinDistOnEdge(v0,v1,edgeGrid,poly,closestP);
        nv.P()=closestP;        
        nv.Q()=0;
        newVertVec.push_back(tri::Index(com.base,&nv));
//        nv.P() = CoS::QLerp(v0,v1);        
      }
      Color4b WedgeInterp(Color4b &c0, Color4b &c1)
      {
          Color4b cc;
          cc.lerp(c0,c1,0.5f);
          return Color4b::Red;
      }
      TexCoord2f WedgeInterp(TexCoord2f &t0, TexCoord2f &t1)
      {
          TexCoord2f tmp;
          assert(t0.n()== t1.n());
          tmp.n()=t0.n();
          tmp.t()=(t0.t()+t1.t())/2.0;
          return tmp;
      }
  };
  
  class EdgePointPred
  {
  public:
      CoM &com;
      KdTree<ScalarType> &kdtree;
      
      EdgePointPred(CoM &_com, KdTree<ScalarType> &_kdtree):com(_com),kdtree(_kdtree) {};
      bool operator()(face::Pos<FaceType> ep) const
      {
        CoordType p0 = ep.V()->P();
        CoordType p1 = ep.VFlip()->P();
        float stepNum=100.0;
        float locEps = Distance(p0,p1)/stepNum;
        for(float j=0;j<stepNum;++j)
        {
          CoordType qp = (p0*j+p1*(stepNum-j))/stepNum;
          unsigned int veInd;
          ScalarType sqdist;
          kdtree.doQueryClosest(qp,veInd,sqdist);
          if(sqrt(sqdist)<locEps) return true;
        }
        return false;
      }
  };
  
  struct EdgePointSplit : public std::unary_function<face::Pos<FaceType> ,  Point3f>
  {
    CoM &com;
    KdTree<ScalarType> &kdtree;
    MeshType &poly;
    
    EdgePointSplit(CoM &_com, KdTree<ScalarType> &_kdtree, MeshType &_poly):com(_com),kdtree(_kdtree),poly(_poly) {};
      void operator()(VertexType &nv, face::Pos<FaceType> ep)
      {
        CoordType p0 = ep.V()->P();
        CoordType p1 = ep.VFlip()->P();
        float stepNum=100.0;
        float locEps = Distance(p0,p1)/stepNum;
        for(float j=0;j<stepNum;++j)
        {
          CoordType qp = (p0*j+p1*(stepNum-j))/stepNum;
          unsigned int veInd;
          ScalarType sqdist;
          kdtree.doQueryClosest(qp,veInd,sqdist);
          if(sqrt(sqdist)<locEps) 
            nv.P() = poly.vert[veInd].P();           
        }
        return;
      }
      Color4b WedgeInterp(Color4b &c0, Color4b &c1)
      {
          Color4b cc;
          cc.lerp(c0,c1,0.5f);
          return Color4b::Red;
      }
      TexCoord2f WedgeInterp(TexCoord2f &t0, TexCoord2f &t1)
      {
          TexCoord2f tmp;
          assert(t0.n()== t1.n());
          tmp.n()=t0.n();
          tmp.t()=(t0.t()+t1.t())/2.0;
          return tmp;
      }
  };
  
  
  
  
  void DumpPlaneMesh(MeshType &poly, std::vector<Plane3f> &planeVec, int i =0)
  {
    MeshType full;
    for(int i=0;i<planeVec.size();++i)
    {
      float radius=edge::Length(poly.edge[i])/2.0;
      MeshType t;
      OrientedDisk(t,16,edge::Center(poly.edge[i]),planeVec[i].Direction(),radius);
      tri::Append<MeshType,MeshType>::Mesh(full,t);
    }
    char buf[100];
    sprintf(buf,"planes%03i.ply",i);
    tri::io::ExporterPLY<MeshType>::Save(full,buf);  
  }

  Plane3f ComputeEdgePlane(VertexType *v0, VertexType *v1)
  {
    Plane3f pl;
    Point3f mid = (v0->P()+v1->P())/2.0;
    Point3f delta = (v1->P()-v0->P()).Normalize();
    Point3f n0 =  (v0->N() ^ delta).Normalize();
    Point3f n1 =  (v1->N() ^ delta).Normalize();
    Point3f n = (n0+n1)/2.0;
    pl.Init(mid,n);
    return pl;
  }


  void ComputePlaneField(MeshType &poly, EdgeGrid &edgeGrid, int ind)
  {
    // First Compute per-edge planes
    std::vector<Plane3f> planeVec(poly.en);
    for(int i=0;i<poly.en;++i)
    {
      planeVec[i] = ComputeEdgePlane(poly.edge[i].V(0), poly.edge[i].V(1));
    }
    
    DumpPlaneMesh(poly,planeVec,ind);
    edgeGrid.Set(poly.edge.begin(), poly.edge.end());    
    
    for(VertexIterator vi=base.vert.begin();vi!=base.vert.end();++vi)
    {
      Point3<ScalarType> p = vi->P();
      
      float minDist=par.gridBailout;
      Point3f closestP;
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,p,par.gridBailout,minDist,closestP);        
      if(cep)
      {
        int ind = tri::Index(poly,cep);        
        vi->Q() = SignedDistancePlanePoint(planeVec[ind],p);
        Point3f proj = planeVec[ind].Projection(p);
        
//        if(Distance(proj,closestP)>edge::Length(*cep))
//        {
//          vi->Q() =1;
//          vi->SetS();
//        }
      }
      else {
        vi->Q() =1;
      }
    } 
  }

  
  void CutAlongPolyLineUsingField(MeshType &poly,EdgeGrid &edgeGrid,std::vector<int> &newVertVec)
{
  QualitySign qsPred(edgeGrid,poly,*this);
  QualitySignSplit qsSplit(edgeGrid,poly,*this,newVertVec);
  tri::UpdateTopology<MeshType>::FaceFace(base);
  tri::RefineE(base,qsSplit,qsPred);    
  tri::UpdateTopology<MeshType>::FaceFace(base);
    
  
  for(int i=0;i<base.fn;++i)
  {
    FaceType *fp = &base.face[i];
    if(!fp->IsD())
    {
    for(int j=0;j<3;++j)
    {
      if(Distance(fp->P0(j),fp->P1(j)) < par.polyDistThr)
      {
        if(face::FFLinkCondition(*fp,j))
        {
//            if(fp->V0(j)->Q()==0) fp->V1(j)->Q()=0;
//            face::FFEdgeCollapse(base,*fp,j);        
            break;
        }
      }
    }
    }
  }
  tri::Allocator<MeshType>::CompactEveryVector(base);
  
  for(int i=0;i<base.fn;++i)
  {
    FaceType *fp = &base.face[i];
    if( (fp->V(0)->Q()==0) &&
        (fp->V(1)->Q()==0) &&
        (fp->V(2)->Q()==0) )
    {
      ScalarType maxDist = 0;
      int maxInd = -1;
      for(int j=0;j<3;++j)
      {
        Point3f closestPt;
        ScalarType d = MinDistOnEdge(fp->P(j),edgeGrid,poly,closestPt);
        if(d>maxDist)
        { 
          maxDist= d;
          maxInd=j;
        }
      }        
//      assert(maxInd!=-1);
//      if(maxInd>=0 && maxDist > par.surfDistThr)
//        fp->V(maxInd)->Q() = maxDist;
    }
  }  
  
    for(int i=0;i<base.fn;++i)
    {
      FaceType *fp = &base.face[i];
      if( (fp->V(0)->Q()>=0) &&
          (fp->V(1)->Q()>=0) &&
          (fp->V(2)->Q()>=0) )
        fp->C() = Color4b::Blue;
      if( (fp->V(0)->Q()<=0) &&
          (fp->V(1)->Q()<=0) &&
          (fp->V(2)->Q()<=0) )
        fp->C() = Color4b::Red;
      if( (fp->V(0)->Q()==0) &&
          (fp->V(1)->Q()==0) &&
          (fp->V(2)->Q()==0) )
        fp->C() = Color4b::Green;
      
      if( (fp->V(0)->Q()>0) &&
          (fp->V(1)->Q()>0) &&
          (fp->V(2)->Q()>0) )
        fp->C() = Color4b::White;
      if( (fp->V(0)->Q()<0) &&
          (fp->V(1)->Q()<0) &&
          (fp->V(2)->Q()<0) )
        fp->C() = Color4b::White;
  }
    tri::AttributeSeam::SplitVertex(base, ExtractVertex, CompareVertex);
    
 }
  
  void WalkAlongPolyLine(MeshType &poly, std::vector<VertexType *> &ptVec)
  {
    // Search a starting vertex
    VertexType *startVert;
    for(int i=0;i<base.vn;++i)
    {
      if(Distance(base.vert[i].P(),ptVec[0]->P()) < par.polyDistThr)
      {
        startVert = &base.vert[i];
        break;
      }
    }
    tri::UpdateTopology<MeshType>::VertexFace(base);
    tri::UpdateTopology<MeshType>::FaceFace(base);
    
    
    
  }
  
    
    
    
//  }
// } 
  
/**
 *
 * 
 */

  void CutWithPolyLine(MeshType &poly)
  {    
    std::vector<int> newVertVec;
    SnapPolyline(poly, &newVertVec);
    tri::io::ExporterPLY<MeshType>::Save(poly,"poly_snapped.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);    
    
    DecomposeNonManifoldPolyline(poly);
    tri::io::ExporterPLY<MeshType>::Save(poly,"poly_manif.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);    
    std::vector< std::vector< int> > ccVec;
    BuildConnectedComponentVectors(poly,ccVec);
    printf("PolyLine of %i edges decomposed into %i manifold components\n",poly.en,ccVec.size());
    Reorient(poly,ccVec);
    char buf[1024];
    for(int i=0;i<ccVec.size();++i)
//      for(int i=0;i<10;++i)
      {
        MeshType subPoly;
        ExtractSubMesh(poly,ccVec[i],subPoly);
        std::vector< VertexType *> ptVec;
        FindTerminalPoints(subPoly,ptVec);
        printf("Component %i (%i edges) has %i terminal points\n",i,subPoly.en, ptVec.size());fflush(stdout);
        SplitMeshWithPoints(base,ptVec,newVertVec);
//        sprintf(buf,"CuttingPoly%02i.ply",i);
//        tri::io::ExporterPLY<MeshType>::Save(subPoly, buf,tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);        
        EdgeGrid edgeGrid;
        ComputePlaneField(subPoly, edgeGrid,i);
        sprintf(buf,"PlaneField%02i.ply",i);
        tri::io::ExporterPLY<MeshType>::Save(base,buf,tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );
        CutAlongPolyLineUsingField(subPoly,edgeGrid,newVertVec);
        sprintf(buf,"PlaneCut%02i.ply",i);
        tri::io::ExporterPLY<MeshType>::Save(base,buf,tri::io::Mask::IOM_FACECOLOR + tri::io::Mask::IOM_VERTQUALITY );  
      }
    
//    printf("Added %i vertices\n",newVertVec.size());
//    for(int i=0;i<newVertVec.size();++i)
//        base.vert[newVertVec[i]].C()=Color4b::Red;
    tri::io::ExporterPLY<MeshType>::Save(base,"base_cut.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );  
    CleanSpuriousDanglingEdges(poly);
    tri::io::ExporterPLY<MeshType>::Save(base,"base_cut_clean.ply",tri::io::Mask::IOM_VERTCOLOR + tri::io::Mask::IOM_VERTQUALITY );  
    
  }

  
  
  
  void SnapPolyline(MeshType &poly, std::vector<int> *newVertVec)
  {
    const float maxDist = base.bbox.Diag()/100.0;
    const ScalarType interpEps = 0.0001;
    int vertSnapCnt=0;
    int edgeSnapCnt=0;
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      float closestDist;
      Point3f closestP,closestN,ip;
      FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,vi->P(),maxDist, closestDist, closestP, closestN,ip);
      assert(f);
      VertexType *closestVp=0;
      int indIp = -1;
      ScalarType minDist = std::numeric_limits<ScalarType>::max();
      ScalarType minIp = minDist;
      for(int i=0;i<3;++i)
      {
        if(Distance(vi->P(),f->P(i))<minDist)
        {
          minDist = Distance(vi->P(),f->P(i));
          closestVp = f->V(i);
        }
        if(minIp > ip[i]) 
        {
          indIp = i;
          minIp=ip[i];
        }          
      }
      assert(closestVp && (indIp!=-1));
      
      
      if(minDist < par.maxSnapThr) {  // First Case: Snap to vertex;
        vi->P() = closestVp->P();
        vertSnapCnt++;
        if(newVertVec)
        newVertVec->push_back(tri::Index(base,closestVp));
      } else {
        if(minIp < interpEps) {       // Second Case: Snap to edge;
          ScalarType T1 = ip[(indIp+1)%3];
          ScalarType T2 = ip[(indIp+2)%3];
          vi->P() = (f->V1(indIp)->P() * T1 + f->V2(indIp)->P() * T2)/(T1+T2);
          edgeSnapCnt++;
        }
      }
    }
    printf("Snapped %i onto vert and %i onto edges\n",vertSnapCnt, edgeSnapCnt);
  }
 

  
  /*
   * Make an edge mesh 1-manifold by splitting all the
   * vertexes that have more than two incident edges
   * 
   * It performs the split in three steps. 
   * - First it collects and counts the vertices to be splitten. 
   * - Then it adds the vertices to the mesh and 
   * - lastly it updates the poly with the newly added vertices. 
   * 
   * singSplitFlag allows to ubersplit each singularity in a number of vertex of the same order of its degree. 
   * This is not really necessary but helps the management of sharp turns in the poly mesh.
   * \todo add corner detection and split.
   */
  
  void DecomposeNonManifoldPolyline(MeshType &poly, bool singSplitFlag = true)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    std::vector<int> degreeVec(poly.vn, 0);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int neededVert=0;
    int delta;
    if(singSplitFlag) delta = 1;
                 else delta = 2;
      
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      std::vector<EdgeType *> starVec;
      edge::VEStarVE(&*vi,starVec);
      degreeVec[tri::Index(poly, *vi)] = starVec.size();
      if(starVec.size()>2)
        neededVert += starVec.size()-delta;
    }
    printf("DecomposeNonManifold Adding %i vert to a polyline of %i vert\n",neededVert,poly.vn);
    VertexIterator firstVi = tri::Allocator<MeshType>::AddVertices(poly,neededVert);
    
    for(size_t i=0;i<degreeVec.size();++i)
    {
      if(degreeVec[i]>2)
      {
        std::vector<EdgeType *> edgeStarVec;
        edge::VEStarVE(&(poly.vert[i]),edgeStarVec);
        assert(edgeStarVec.size() == degreeVec[i]);
        for(size_t j=delta;j<edgeStarVec.size();++j)
        {
          EdgeType *ep = edgeStarVec[j];
          int ind; // index of the vertex to be changed
          if(tri::Index(poly,ep->V(0)) == i) ind = 0;
              else ind = 1;
  
          ep->V(ind) = &*firstVi;
          ep->V(ind)->P() = poly.vert[i].P();
          ep->V(ind)->N() = poly.vert[i].N();
          ++firstVi;
        }
      }
    }
    assert(firstVi == poly.vert.end());
  }
  
  /*
   * Given a manifold edgemesh it returns the boundary/terminal endpoints.
   */
  static void FindTerminalPoints(MeshType &poly, std::vector<VertexType *> &vec)
  {
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      if(edge::VEDegree<EdgeType>(&*vi)==1)
         vec.push_back(&*vi);
    }
  }
  
  // This function will decompose the input edge mesh into a set of 
  // connected components.
  // the vector will contain, for each connected component, a vector with all the edge indexes. 
  void BuildConnectedComponentVectors(MeshType &poly, std::vector< std::vector< int> > &ccVec)
  {
    UpdateTopology<MeshType>::VertexEdge(poly);
    for(size_t i=0;i<poly.vn;++i)
    {
      assert(edge::VEDegree<EdgeType>(&(poly.vert[i])) <=2);
    }
    
    tri::UpdateTopology<MeshType>::EdgeEdge(poly);
    tri::UpdateFlags<MeshType>::EdgeClearV(poly);

    int visitedEdgeNum=0 ;
    int ccCnt=0;
    EdgeIterator eIt = poly.edge.begin();  
    
    while(visitedEdgeNum < poly.en)
    {
      ccVec.resize(ccVec.size()+1);
      while(eIt->IsV()) ++eIt;
//      printf("Starting component from edge %i\n",tri::Index(poly,&*eIt));
      assert(eIt != poly.edge.end());
      edge::Pos<EdgeType> startPos(&*eIt,0);
      edge::Pos<EdgeType> curPos(&*eIt,0);
      do
      {
//        printf("(%i %i %i)-",tri::Index(poly,curPos.VFlip()), tri::Index(poly,curPos.E()) ,tri::Index(poly,curPos.V()));
        curPos.NextE();
      } 
      while(curPos!=startPos && !curPos.IsBorder()) ;
      
      curPos.FlipV();
      assert(!curPos.IsBorder());
      do
      {
//         printf("<%i %i %i>-",tri::Index(poly,curPos.VFlip()), tri::Index(poly,curPos.E()) ,tri::Index(poly,curPos.V()));
        curPos.E()->SetV();
        visitedEdgeNum++;
        ccVec[ccCnt].push_back(tri::Index(poly,curPos.E()));        
        curPos.NextE();      
      } while(!curPos.E()->IsV());
      printf("Completed component %i of %i edges\n",ccCnt, ccVec[ccCnt].size());
      ccCnt++;      
    }    
  }
  
  // This function will decompose the input edge mesh into a set of 
  // connected components.
  // the vector will contain, for each connected component, a vector with all the edge indexes. 
  void BuildConnectedComponentVectorsOld(MeshType &poly, std::vector< std::vector< int> > &ccVec)
  {
    tri::UpdateTopology<MeshType>::EdgeEdge(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    tri::UpdateFlags<MeshType>::EdgeClearV(poly);
    
    int visitedEdgeNum=0 ;
    int ccCnt=0;
    
    EdgeIterator eIt = poly.edge.begin();  
    while(visitedEdgeNum < poly.en)
    {
      ccVec.resize(ccVec.size()+1);
      while((eIt != poly.edge.end()) && eIt->IsV()) ++eIt;
      EdgeType *startE = &*eIt;
      
      EdgeType *curEp = &*eIt;
      int     curEi = 0;
      printf("Starting Visit of connected Component %i from edge %i\n",ccCnt,tri::Index(poly,*eIt));
      while( (curEp->EEp(curEi) != startE) &&
             (curEp->EEp(curEi) != curEp) )
      {
        EdgeType *nextEp = curEp->EEp(curEi);
        int nextEi =  (curEp->EEi(curEi)+1)%2;
        curEp = nextEp;
        curEi = nextEi;
      }   
  
      curEp->SetV();
      curEi = (curEi +1)%2; // Flip the visit direction!
      visitedEdgeNum++;
      ccVec[ccCnt].push_back(tri::Index(poly,curEp));        
      while(! curEp->EEp(curEi)->IsV())
      {
        EdgeType *nextEp = curEp->EEp(curEi);
        int nextEi =  (curEp->EEi(curEi)+1)%2;
        curEp->SetV();
        curEp = nextEp;
        curEi = nextEi;
        curEp->V(curEi)->C() = Color4b::Scatter(30,ccCnt%30);
        if(!curEp->IsV()) {
          ccVec[ccCnt].push_back(tri::Index(poly,curEp));      
          visitedEdgeNum++;      
        }          
      }
      printf("Completed visit of component of size %i\n",ccVec[ccCnt].size());
      ccCnt++;    
    }
    printf("en %i - VisitedEdgeNum = %i\n",poly.en, visitedEdgeNum);
    
  }
  
  void ExtractSubMesh(MeshType &poly, std::vector<int> &ind, MeshType &subPoly)
  {
    subPoly.Clear();
    std::vector<int> remap(poly.vert.size(),-1);
    for(size_t i=0;i<ind.size();++i)
    {
      int v0 = tri::Index(poly,poly.edge[ind[i]].V(0));
      int v1 = tri::Index(poly,poly.edge[ind[i]].V(1));
      if(remap[v0]==-1) 
      {
        remap[v0]=subPoly.vn;
        tri::Allocator<MeshType>::AddVertex(subPoly,poly.vert[v0].P(),poly.vert[v0].N());      
      }
      if(remap[v1]==-1) 
      {
        remap[v1]=subPoly.vn;
        tri::Allocator<MeshType>::AddVertex(subPoly,poly.vert[v1].P(),poly.vert[v1].N());      
      }
      tri::Allocator<MeshType>::AddEdge(subPoly, &subPoly.vert[remap[v0]], &subPoly.vert[remap[v1]]);   
    }  
  }
  
  // It takes a vector of vector of connected components and cohorently reorient each one of them.
  // it usese the EE adjacency and requires that the input edgemesh is 1manifold. 
  
  void Reorient(MeshType &poly, std::vector< std::vector< int> > &ccVec)
  {
    UpdateTopology<MeshType>::VertexEdge(poly);
    for(size_t i=0;i<poly.vn;++i)
    {
      assert(edge::VEDegree<EdgeType>(&(poly.vert[i])) <=2);
    }
        UpdateTopology<MeshType>::EdgeEdge(poly);
    
    for(size_t i=0;i<ccVec.size();++i)
    {
      std::vector<bool> toFlipVec(ccVec[i].size(),false);
      
      for(int j=0;j<ccVec[i].size();++j)
      {
        EdgeType *cur  = & poly.edge[ccVec[i][j]];
        EdgeType *prev;
        if(j==0) 
        {
          if(cur->EEp(0) == cur) 
            prev = cur; // boundary
          else 
            prev = & poly.edge[ccVec[i].back()]; // cc is a loop
        }
        else prev = & poly.edge[ccVec[i][j-1]];
        
        if(cur->EEp(0) != prev)
        {
          toFlipVec[j] = true;
          assert(cur->EEp(1) == prev || j==0);
        }      
      }
      for(int j=0;j<ccVec[i].size();++j)
        if(toFlipVec[j]) 
          std::swap(poly.edge[ccVec[i][j]].V(0), poly.edge[ccVec[i][j]].V(1));    
    }
  }
  
  
  /*
   * Given a mesh and vector of vertex pointers it splits the mesh with the given points.
   * To avoid degeneracies it snaps the splitting points that came nearest than a given
   * threshold to the edge or the vertex of the closest triangle. When an input vertex
   * is snapped the passed vertex position is modiried accordingly
   * (this is the reason for having a vector of vertex pointers instead just a vector of points)
   *
   */
  void SplitMeshWithPoints(MeshType &m, std::vector<VertexType *> &vec, std::vector<int> &newVertVec)
  {
    int faceToAdd=0;
    int vertToAdd=0;
    
    // For each splitting point we save the index of the face to be splitten and the "kind" of split to do:
    // 3 -> means classical 1 to 3 face split
    // 2 -> means edge split.
    // 0 -> means no need of a split (e.g. the point is coincident with a vertex)
    
    std::vector< std::pair<int,int> > toSplitVec(vec.size(), std::make_pair(0,0));
    MeshGrid uniformGrid;
    uniformGrid.Set(m.face.begin(), m.face.end());
   
    for(size_t i =0; i<vec.size();++i)
    {
      Point3f newP = vec[i]->P();
      float closestDist;
      Point3f closestP;
      FaceType *f = vcg::tri::GetClosestFaceBase(m,uniformGrid,newP,par.gridBailout, closestDist, closestP);
      assert(f);
      VertexType *closestVp=0;
      ScalarType minDist = std::numeric_limits<ScalarType>::max();
      for(int i=0;i<3;++i) {
        if(Distance(newP,f->P(i))<minDist)
        {
          minDist = Distance(newP,f->P(i));
          closestVp = f->V(i);
        }
      }
      assert(closestVp);
      if(minDist < par.maxSnapThr)  {
           vec[i]->P() = closestVp->P();
      }
        else
      {
        toSplitVec[i].first = tri::Index(m,f);
        toSplitVec[i].second = 3;
        faceToAdd += 2;
        vertToAdd += 1;
      }
    }
//    printf("Splitting with %i points: adding %i faces and %i vertices\n",vec.size(), faceToAdd,vertToAdd);
    FaceIterator newFi = tri::Allocator<MeshType>::AddFaces(m,faceToAdd);
    VertexIterator newVi = tri::Allocator<MeshType>::AddVertices(m,vertToAdd);
    
    tri::UpdateColor<MeshType>::PerFaceConstant(m,Color4b::White);
    
    for(size_t i =0; i<vec.size();++i)
    {
      if(toSplitVec[i].second==3)
      {        
        FaceType *fp0 = &m.face[toSplitVec[i].first];
        FaceType *fp1 = &*newFi; newFi++;
        FaceType *fp2 = &*newFi; newFi++;
        VertexType *vp = &*(newVi++);
        newVertVec.push_back(tri::Index(base,vp));
        vp->P() = vec[i]->P();
        VertexType *vp0 = fp0->V(0);
        VertexType *vp1 = fp0->V(1);
        VertexType *vp2 = fp0->V(2);
        
        fp0->V(0) = vp0; fp0->V(1) = vp1; fp0->V(2) = vp;
        fp1->V(0) = vp1; fp1->V(1) = vp2; fp1->V(2) = vp;
        fp2->V(0) = vp2; fp2->V(1) = vp0; fp2->V(2) = vp;
        
        fp0->C() = Color4b::Green;
        fp1->C() = Color4b::Green;
        fp2->C() = Color4b::Green;
      }
    }
  }
  
  
  void Init()
  {
    // Construction of the uniform grid
    uniformGrid.Set(base.face.begin(), base.face.end());
    UpdateNormal<MeshType>::PerFaceNormalized(base);
    UpdateTopology<MeshType>::FaceFace(base);    
    uniformGrid.Set(base.face.begin(), base.face.end());    
  }
  
  
  void Simplify( MeshType &poly)
  {
    int startEn = poly.en;
    Distribution<ScalarType> hist;
    for(int i =0; i<poly.en;++i) 
      hist.Add(edge::Length(poly.edge[i]));
        
    UpdateTopology<MeshType>::VertexEdge(poly);
    
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if( (starVecVp.size()==2) && (!poly.vert[i].IsS()))
      {      
        float newSegLen = Distance(starVecVp[0]->P(), starVecVp[1]->P());
        Segment3f seg(starVecVp[0]->P(),starVecVp[1]->P());
        float segDist;
        Point3f closestPSeg;
        SegmentPointDistance(seg,poly.vert[i].cP(),closestPSeg,segDist);
        Point3f fp,fn;
        float maxSurfDist = MaxSegDist(starVecVp[0], starVecVp[1],fp,fn);
        
        if(maxSurfDist < par.surfDistThr && (newSegLen < par.maxSimpEdgeLen) )
        {
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));
          
        }
      }
    }
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    printf("Simplify %5i -> %5i (total len %5.2f)\n",startEn,poly.en,hist.Sum());
  }
  
  void EvaluateHausdorffDistance(MeshType &poly, Distribution<ScalarType> &dist)
  {
    dist.Clear();
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    tri::UpdateQuality<MeshType>::VertexConstant(poly,0);
    for(int i =0; i<poly.edge.size();++i)
    {      
      Point3f farthestP, farthestN;      
      float maxDist = MaxSegDist(poly.edge[i].V(0),poly.edge[i].V(1), farthestP, farthestN, &dist);      
      poly.edge[i].V(0)->Q()+= maxDist;
      poly.edge[i].V(1)->Q()+= maxDist;
    }
    for(int i=0;i<poly.vn;++i)
    {
      ScalarType deg = edge::VEDegree<EdgeType>(&poly.vert[i]);
      poly.vert[i].Q()/=deg;
    }
    tri::UpdateColor<MeshType>::PerVertexQualityRamp(poly,0,dist.Max());    
  }
  
  
  
  // Given a segment find the maximum distance from it to the original surface. 
  float MaxSegDist(VertexType *v0, VertexType *v1, Point3f &farthestPointOnSurf, Point3f &farthestN, Distribution<ScalarType> *dist=0)
  {
    float maxSurfDist = 0;
    const float sampleNum = 10;
    const float maxDist = base.bbox.Diag()/10.0;
    for(float k = 1;k<sampleNum;++k)
    {
      float surfDist;
      Point3f closestPSurf;
      Point3f samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,samplePnt,maxDist, surfDist, closestPSurf);        
      if(dist)
        dist->Add(surfDist);
      assert(f);
      if(surfDist > maxSurfDist)
      {
        maxSurfDist = surfDist;
        farthestPointOnSurf = closestPSurf;
        farthestN = f->N();
      }
    }
    return maxSurfDist;
  }
  
  void Refine(MeshType &poly, bool uniformFlag = false)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    int startEdgeSize = poly.en;
    for(int i =0; i<startEdgeSize;++i)
    {
      EdgeType &ei = poly.edge[i];
      if(edge::Length(ei)>par.minRefEdgeLen)  
      {      
        Point3f farthestP, farthestN;
        float maxDist = MaxSegDist(ei.V(0),ei.V(1),farthestP, farthestN);
        if(maxDist > par.surfDistThr)  
        {
          edge::VEEdgeSplit(poly, &ei, farthestP, farthestN); 
        }
        else if(uniformFlag)
        {
         edge::VEEdgeSplit(poly,&ei,(ei.P(0)+ei.P(1))/2.0,(ei.V(0)->N()+ei.V(1)->N())/2.0); 
        }
      }
    }
//    tri::Allocator<MeshType>::CompactEveryVector(poly);
    printf("Refine %i -> %i\n",startEdgeSize,poly.en);fflush(stdout);
  }
  
  // Take a poly and find the position of edge -edge crossing.
  void RefineBaseMesh(MeshType &poly)
  {
    std::vector<CoordType> newPtVec;
    for(int i=0;i<poly.edge.size();++i)
    {
      {
        Point3f p0 = poly.edge[i].P(0);
        Point3f p1 = poly.edge[i].P(1);
        float stepNum = 10;
        FaceType *lastFace=0; 
        Point3f lastqp;
        Segment3f seg(p0, p1);        
        for(float j=1;j<stepNum;++j) 
        {
          Point3f qp = (p0*j + p1*(stepNum-j))/stepNum;
          const float maxDist = base.bbox.Diag()/20.0;
          float surfDist;
          Point3f closestPSurf, normSur, ip;
          FaceType *fp = vcg::tri::GetClosestFaceBase(base,uniformGrid,qp,maxDist, surfDist, closestPSurf, normSur, ip);                  
          if((fp != lastFace) && (lastFace!=0)  && (fp!=0)  )
          {
            for(int ei =0 ; ei<3;++ei)
            {
              if(fp->FFp(ei) == lastFace)
              {
                Point3f ps0,ps1;
                float segDist;
                bool par;
                Segment3f faceSeg(fp->P0(ei),fp->P1(ei));
                float angDeg = fabs(math::ToDeg(Angle(faceSeg.Direction(),seg.Direction())));
                
                SegmentSegmentDistance(seg,faceSeg,segDist,par,ps0,ps1);
                if(!par && (angDeg > 5 && angDeg<175)) newPtVec.push_back(ps1);             
              }
            }
          }
          lastFace=fp;        
          lastqp = qp;
        }
      }      
    }  
    MeshType meshPt; tri::BuildMeshFromCoordVector(meshPt,newPtVec);
    tri::io::ExporterPLY<MeshType>::Save(meshPt,"newPoints.ply");  

    VertexConstDataWrapper<MeshType > vdw(meshPt);
    KdTree<ScalarType> kdtree(vdw);
    
    EdgePointPred epPred(*this,kdtree);
    EdgePointSplit epSplit(*this,kdtree,meshPt);
    tri::UpdateTopology<MeshType>::FaceFace(base);
    tri::RefineE(base,epSplit,epPred);    
    tri::UpdateTopology<MeshType>::FaceFace(base);
    tri::io::ExporterPLY<MeshType>::Save(base,"newbase.ply");  
  }
  
  
  /**
   * @brief SmoothProject
   * @param poly
   * @param iterNum
   * @param smoothWeight  [0..1] range;  
   * @param projectWeight [0..1] range;
   * 
   * The very important function to adapt a polyline onto the base mesh
   * The projection process must be done slowly to guarantee some empirical convergence...
   * For each iteration it choose a new position of each vertex of the polyline. 
   * The new position is a blend between the smoothed position, the closest point on the surface and the original position. 
   * You need a good balance...
   * after each iteration the polyline is refined and simplified. 
   */
  void SmoothProject(MeshType &poly, int iterNum, ScalarType smoothWeight, ScalarType projectWeight)
  {
    tri::RequireCompactness(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    printf("SmoothProject: Selected vert num %i\n",tri::UpdateSelection<MeshType>::VertexCount(poly));
    assert(poly.en>0 && base.fn>0);
    for(int k=0;k<iterNum;++k)
    {
      std::vector<Point3f> posVec(poly.vn,Point3f(0,0,0));
      std::vector<int>     cntVec(poly.vn,0);
  
      for(int i =0; i<poly.en;++i)
      {
        for(int j=0;j<2;++j)
        {
          int vertInd = tri::Index(poly,poly.edge[i].V(j));
          posVec[vertInd] += poly.edge[i].V1(j)->P();
          cntVec[vertInd] += 1;
        }
      }
      
      const float maxDist = base.bbox.Diag()/10.0;
      for(int i=0; i<poly.vn; ++i)
        if(!poly.vert[i].IsS())
        {
          Point3f smoothPos = (poly.vert[i].P() + posVec[i])/float(cntVec[i]+1);
          
          Point3f newP = poly.vert[i].P()*(1.0-smoothWeight) + smoothPos *smoothWeight;
          
          Point3f delta =  newP - poly.vert[i].P();
          if(delta.Norm() > par.maxSmoothDelta) 
          {
            newP =  poly.vert[i].P() + ( delta / delta.Norm()) * maxDist*0.5;
          }
          
          float minDist;
          Point3f closestP;
          FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,newP,maxDist, minDist, closestP);
          assert(f);
          poly.vert[i].P() = newP*(1.0-projectWeight) +closestP*projectWeight;
          poly.vert[i].N() = f->N();
        }
      
//      Refine(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    Refine(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    Simplify(poly);
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      int dupVertNum = Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(dupVertNum) {
        printf("****REMOVED %i Duplicated\n",dupVertNum);
        tri::Allocator<MeshType>::CompactEveryVector(poly);
        tri::UpdateTopology<MeshType>::VertexEdge(poly);
      }
    }
  }  
   
  
};
} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_CURVE_ON_SURF_H
