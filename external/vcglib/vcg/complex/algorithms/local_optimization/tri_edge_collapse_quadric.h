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

#ifndef __VCG_TRIMESHCOLLAPSE_QUADRIC__
#define __VCG_TRIMESHCOLLAPSE_QUADRIC__

#include<vcg/math/quadric.h>
#include<vcg/simplex/face/pos.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include<vcg/complex/algorithms/local_optimization.h>
#include<vcg/complex/algorithms/stat.h>


namespace vcg{
namespace tri{




/**
  This class describe Quadric based collapse operation.

    Requirements:

    Vertex
    must have:
   incremental mark
   VF topology

    must have:
        members

      QuadricType Qd();

            ScalarType W() const;
                A per-vertex Weight that can be used in simplification
                lower weight means that error is lowered,
                standard: return W==1.0

            void Merge(MESH_TYPE::vertex_type const & v);
                Merges the attributes of the current vertex with the ones of v
                (e.g. its weight with the one of the given vertex, the color ect).
                Standard: void function;

      OtherWise the class should be templated with a static helper class that helps to retrieve these functions.
      If the vertex class exposes these functions a default static helper class is provided.

*/
        //**Helper CLASSES**//
        template <class VERTEX_TYPE>
        class QInfoStandard
        {
        public:
      QInfoStandard(){}
      static void Init(){}
      static math::Quadric<double> &Qd(VERTEX_TYPE &v) {return v.Qd();}
      static math::Quadric<double> &Qd(VERTEX_TYPE *v) {return v->Qd();}
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE * /*v*/) {return 1.0;}
      static typename VERTEX_TYPE::ScalarType W(VERTEX_TYPE &/*v*/) {return 1.0;}
      static void Merge(VERTEX_TYPE & /*v_dest*/, VERTEX_TYPE const & /*v_del*/){}
        };


class TriEdgeCollapseQuadricParameter : public BaseParameterClass
{
public:
  double    BoundaryWeight;
  double    CosineThr;
  bool      FastPreserveBoundary;
  bool      NormalCheck;
  double    NormalThrRad;
  bool      OptimalPlacement;
  bool      PreserveTopology;
  bool      PreserveBoundary;
  double    QuadricEpsilon;
  bool      QualityCheck;
  bool      QualityQuadric; // During the initialization manage all the edges as border edges adding a set of additional quadrics that are useful mostly for keeping face aspect ratio good.
  double    QualityThr;     // Collapsed that generate faces with quality LOWER than this value are penalized. So 
  bool      QualityWeight;
  double    QualityWeightFactor;
  double    ScaleFactor;
  bool      ScaleIndependent;
  bool      UseArea;
  bool      UseVertexWeight;

  void SetDefaultParams()
  {
    BoundaryWeight=.5;
    CosineThr=cos(M_PI/2);
    FastPreserveBoundary=false;
    NormalCheck=false;
    NormalThrRad=M_PI/2;
    OptimalPlacement=true;
    PreserveBoundary = false;
    PreserveTopology = false;
    QuadricEpsilon =1e-15;
    QualityCheck=true;
    QualityQuadric=false;
    QualityThr=.3;        // higher the value -> better the quality of the accepted triangles
    QualityWeight=false;
    QualityWeightFactor=100.0;
    ScaleFactor=1.0;
    ScaleIndependent=true;
    UseArea=true;
    UseVertexWeight=false;
  }

  TriEdgeCollapseQuadricParameter() {this->SetDefaultParams();}
};


template<class TriMeshType, class VertexPair, class MYTYPE, class HelperType = QInfoStandard<typename TriMeshType::VertexType> >
class TriEdgeCollapseQuadric: public TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE>
{
public:
  typedef typename vcg::tri::TriEdgeCollapse< TriMeshType, VertexPair, MYTYPE > TEC;
  typedef typename TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapType HeapType;
  typedef typename TriEdgeCollapse<TriMeshType, VertexPair, MYTYPE>::HeapElem HeapElem;
  typedef typename TriMeshType::CoordType CoordType;
  typedef typename TriMeshType::ScalarType ScalarType;
  typedef typename TriMeshType::FaceType FaceType;
  typedef typename TriMeshType::VertexType VertexType;
  typedef typename TriMeshType::VertexIterator VertexIterator;
  typedef typename TriMeshType::FaceIterator FaceIterator;
  typedef typename vcg::face::VFIterator<FaceType> VFIterator;
  typedef  math::Quadric< double > QuadricType;
  typedef TriEdgeCollapseQuadricParameter QParameter;
  typedef HelperType QH;
  
  
  // puntatori ai vertici che sono stati messi non-w per preservare il boundary
  static std::vector<typename TriMeshType::VertexPointer>  & WV(){
    static std::vector<typename TriMeshType::VertexPointer> _WV; return _WV;
  }
  
  inline TriEdgeCollapseQuadric(){}
  
  inline TriEdgeCollapseQuadric(const VertexPair &p, int i, BaseParameterClass *pp)
  {
    this->localMark = i;
    this->pos=p;
    this->_priority = ComputePriority(pp);
  }
  

    inline bool IsFeasible(BaseParameterClass *_pp){
      QParameter *pp=(QParameter *)_pp;
      if(!pp->PreserveTopology) return true;

      bool res = ( EdgeCollapser<TriMeshType, VertexPair>::LinkConditions(this->pos) );
      if(!res) ++( TEC::FailStat::LinkConditionEdge() );
      return res;
    }

    CoordType ComputePosition(BaseParameterClass *_pp)
    {
      QParameter *pp=(QParameter *)_pp;
      CoordType newPos = (this->pos.V(0)->P()+this->pos.V(1)->P())/2.0;
      if(pp->OptimalPlacement)
      {
        if((QH::Qd(this->pos.V(0)).Apply(newPos) + QH::Qd(this->pos.V(1)).Apply(newPos)) > 200.0*pp->QuadricEpsilon)
          newPos = ComputeMinimal();
      }      
      else newPos=this->pos.V(1)->P();      
      return newPos;
    }

    void Execute(TriMeshType &m, BaseParameterClass *_pp)
  {
    QH::Qd(this->pos.V(1))+=QH::Qd(this->pos.V(0));
    EdgeCollapser<TriMeshType,VertexPair>::Do(m, this->pos, ComputePosition(_pp)); // v0 is deleted and v1 take the new position
  }



    // Final Clean up after the end of the simplification process
    static void Finalize(TriMeshType &m, HeapType& /*h_ret*/, BaseParameterClass *_pp)
    {
      QParameter *pp=(QParameter *)_pp;

      // If we had the boundary preservation we should clean up the writable flags
      if(pp->FastPreserveBoundary)
      {
        typename 	TriMeshType::VertexIterator  vi;
          for(vi=m.vert.begin();vi!=m.vert.end();++vi)
          if(!(*vi).IsD()) (*vi).SetW();
      }
      if(pp->PreserveBoundary)
      {
        typename 	std::vector<typename TriMeshType::VertexPointer>::iterator wvi;
        for(wvi=WV().begin();wvi!=WV().end();++wvi)
          if(!(*wvi)->IsD()) (*wvi)->SetW();
      }
    }

  static void Init(TriMeshType &m, HeapType &h_ret, BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;

  typename 	TriMeshType::VertexIterator  vi;
  typename 	TriMeshType::FaceIterator  pf;

  pp->CosineThr=cos(pp->NormalThrRad);

  vcg::tri::UpdateTopology<TriMeshType>::VertexFace(m);
  vcg::tri::UpdateFlags<TriMeshType>::FaceBorderFromVF(m);

  if(pp->FastPreserveBoundary)
    {
      for(pf=m.face.begin();pf!=m.face.end();++pf)
      if( !(*pf).IsD() && (*pf).IsW() )
        for(int j=0;j<3;++j)
          if((*pf).IsB(j))
          {
            (*pf).V(j)->ClearW();
            (*pf).V1(j)->ClearW();
          }
      }

  if(pp->PreserveBoundary)
    {
      WV().clear();
      for(pf=m.face.begin();pf!=m.face.end();++pf)
      if( !(*pf).IsD() && (*pf).IsW() )
        for(int j=0;j<3;++j)
          if((*pf).IsB(j))
          {
            if((*pf).V(j)->IsW())  {(*pf).V(j)->ClearW(); WV().push_back((*pf).V(j));}
            if((*pf).V1(j)->IsW()) {(*pf).V1(j)->ClearW();WV().push_back((*pf).V1(j));}
          }
      }

    InitQuadric(m,pp);

  // Initialize the heap with all the possible collapses
    if(IsSymmetric(pp))
    { // if the collapse is symmetric (e.g. u->v == v->u)
      for(vi=m.vert.begin();vi!=m.vert.end();++vi)
        if(!(*vi).IsD() && (*vi).IsRW())
            {
                vcg::face::VFIterator<FaceType> x;
                for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x){
                  x.V1()->ClearV();
                  x.V2()->ClearV();
                }
                for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++x )
                {
                  assert(x.F()->V(x.I())==&(*vi));
                  if((x.V0()<x.V1()) && x.V1()->IsRW() && !x.V1()->IsV()){
                        x.V1()->SetV();
                        h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
                        }
                  if((x.V0()<x.V2()) && x.V2()->IsRW()&& !x.V2()->IsV()){
                        x.V2()->SetV();
                        h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp )));
                      }
                }
            }
    }
        else
    { // if the collapse is A-symmetric (e.g. u->v != v->u)
            for(vi=m.vert.begin();vi!=m.vert.end();++vi)
                if(!(*vi).IsD() && (*vi).IsRW())
                    {
                        vcg::face::VFIterator<FaceType> x;
                        UnMarkAll(m);
                        for( x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F()!=0; ++ x)
                        {
                            assert(x.F()->V(x.I())==&(*vi));
                            if(x.V()->IsRW() && x.V1()->IsRW() && !IsMarked(m,x.F()->V1(x.I()))){
                    h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V1()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
                                        }
                            if(x.V()->IsRW() && x.V2()->IsRW() && !IsMarked(m,x.F()->V2(x.I()))){
                    h_ret.push_back( HeapElem( new MYTYPE( VertexPair (x.V(),x.V2()),TriEdgeCollapse< TriMeshType,VertexPair,MYTYPE>::GlobalMark(),_pp)));
                                    }
                        }
                    }
      }
}
  static float HeapSimplexRatio(BaseParameterClass *_pp) {return IsSymmetric(_pp)?5.0f:9.0f;}
  static bool IsSymmetric(BaseParameterClass *_pp) {return ((QParameter *)_pp)->OptimalPlacement;}
  static bool IsVertexStable(BaseParameterClass *_pp) {return !((QParameter *)_pp)->OptimalPlacement;}


/** Evaluate the priority (error) for an edge collapse
  *
  * It simulate the collapse and compute the quadric error 
  * generated by this collapse. This error is weighted with 
  * - aspect ratio of involved triangles
  * - normal variation
  */
  ScalarType ComputePriority(BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    std::vector<CoordType> onVec; // vector with incident faces original normals 
    VertexType * v[2];
    v[0] = this->pos.V(0);
    v[1] = this->pos.V(1);
    
    if(pp->NormalCheck){ // Compute maximal normal variation
      // store the old normals for non-collapsed face in v0
      for(VFIterator x(v[0]); !x.End(); ++x )	 // for all faces in v0
        if( x.V1()!=v[1] && x.V2()!=v[1] ) // skip faces with v1
          onVec.push_back(TriangleNormal(*x.F()).Normalize());
      // store the old normals for non-collapsed face in v1
      for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
        if( x.V1()!=v[0] && x.V2()!=v[0] ) // skip faces with v0
          onVec.push_back(TriangleNormal(*x.F()).Normalize());
    }
    
    //// Move the two vertexes into new position (storing the old ones)
    CoordType OldPos0=v[0]->P();
    CoordType OldPos1=v[1]->P();
    CoordType newPos = ComputePosition(_pp);      
    v[0]->P() = v[1]->P() = newPos;
    
    //// Rescan faces and compute quality and difference between normals
    int i=0;
    double MinCos  = std::numeric_limits<double>::max(); // minimo coseno di variazione di una normale della faccia
    // (e.g. max angle) Mincos varia da 1 (normali coincidenti) a
    // -1 (normali opposte);
    double MinQual = std::numeric_limits<double>::max();
    for(VFIterator x(v[0]); !x.End(); ++x )  // for all faces in v0
      if( x.V1()!=v[1] && x.V2()!=v[1] )     // skiping faces with v1
      {
        if(pp->NormalCheck){
          CoordType nn=NormalizedTriangleNormal(*x.F());
          double ndiff=nn.dot(onVec[i++]);
          MinCos=std::min(MinCos,ndiff);
        }
        if(pp->QualityCheck){ 
          double qt= QualityFace(*x.F());
          MinQual=std::min(MinQual,qt);
        }
      }
    for(VFIterator x(v[1]); !x.End(); ++x )	 // for all faces in v1
      if( x.V1()!=v[0] && x.V2()!=v[0] ) // skip faces with v0
      {
        if(pp->NormalCheck){
          CoordType nn=NormalizedTriangleNormal(*x.F());
          double ndiff=nn.dot(onVec[i++]);
          MinCos=std::min(MinCos,ndiff);
        }
        if(pp->QualityCheck){
          double qt= QualityFace(*x.F());
          MinQual=std::min(MinQual,qt);
        }
      }
    
    QuadricType qq=QH::Qd(v[0]);
    qq+=QH::Qd(v[1]);

    double QuadErr = pp->ScaleFactor*qq.Apply(Point3d::Construct(v[1]->P()));
    
    assert(!math::IsNAN(QuadErr));
    // All collapses involving triangles with quality larger than <QualityThr> have no penalty;
    if(MinQual>pp->QualityThr) MinQual=pp->QualityThr;
    
    if(pp->NormalCheck){     
      // All collapses where the normal vary less than <NormalThr> (e.g. more than CosineThr)
      // have no penalty
      if(MinCos>pp->CosineThr) MinCos=pp->CosineThr;
      MinCos=fabs((MinCos+1)/2.0); // Now it is in the range 0..1 with 0 very dangerous!
    }

    QuadErr= std::max(QuadErr,pp->QuadricEpsilon);
    if(QuadErr <= pp->QuadricEpsilon) 
    {
      QuadErr = - 1/Distance(OldPos0,OldPos1);  
    }

    if( pp->UseVertexWeight ) QuadErr *= (QH::W(v[1])+QH::W(v[0]))/2;
    
    ScalarType error;
    if(!pp->QualityCheck && !pp->NormalCheck) error = (ScalarType)(QuadErr);
    if( pp->QualityCheck && !pp->NormalCheck) error = (ScalarType)(QuadErr / MinQual);
    if(!pp->QualityCheck &&  pp->NormalCheck) error = (ScalarType)(QuadErr / MinCos);
    if( pp->QualityCheck &&  pp->NormalCheck) error = (ScalarType)(QuadErr / (MinQual*MinCos));
    
    // Restore old position of v0 and v1
    v[0]->P()=OldPos0;
    v[1]->P()=OldPos1;
    
    this->_priority = error;
    return this->_priority;
  }
  
//
//static double MaxError() {return 1e100;}
//
  inline void AddCollapseToHeap(HeapType & h_ret, VertexType *v0, VertexType *v1, BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;    
    h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v0,v1), this->GlobalMark(),_pp)));
    std::push_heap(h_ret.begin(),h_ret.end());
    if(!IsSymmetric(pp)){
      h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v1,v0), this->GlobalMark(),_pp)));
      std::push_heap(h_ret.begin(),h_ret.end());
    }
  }
  
  inline  void UpdateHeap(HeapType & h_ret, BaseParameterClass *_pp)
  {
    this->GlobalMark()++;
    VertexType *v[2];
    v[0]= this->pos.V(0);
    v[1]= this->pos.V(1);
    v[1]->IMark() = this->GlobalMark();

    // First loop around the surviving vertex to unmark the Visit flags
    for(VFIterator vfi(v[1]); !vfi.End(); ++vfi ) {
      vfi.V1()->ClearV();
      vfi.V2()->ClearV();
    }

    // Second Loop
    for(VFIterator vfi(v[1]); !vfi.End(); ++vfi ) {
      if( !(vfi.V1()->IsV()) && vfi.V1()->IsRW())
      {
        vfi.V1()->SetV();
        AddCollapseToHeap(h_ret,vfi.V0(),vfi.V1(),_pp);
      }
      if(  !(vfi.V2()->IsV()) && vfi.V2()->IsRW())
      {
        vfi.V2()->SetV();
        AddCollapseToHeap(h_ret,vfi.V2(),vfi.V0(),_pp);
      }
      if(vfi.V1()->IsRW() && vfi.V2()->IsRW() )
        AddCollapseToHeap(h_ret,vfi.V1(),vfi.V2(),_pp);
    } // end second loop around surviving vertex.
  }

  static void InitQuadric(TriMeshType &m,BaseParameterClass *_pp)
  {
    QParameter *pp=(QParameter *)_pp;
    QH::Init();
    //	m.ClearFlags();
    for(VertexIterator pv=m.vert.begin();pv!=m.vert.end();++pv)		// Azzero le quadriche
      if( ! (*pv).IsD() && (*pv).IsW())
        QH::Qd(*pv).SetZero();    
    
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if( !(*fi).IsD() && (*fi).IsR() )
        if((*fi).V(0)->IsR() &&(*fi).V(1)->IsR() &&(*fi).V(2)->IsR())
        {
          Plane3<ScalarType,false> facePlane;
          facePlane.SetDirection( ( (*fi).V(1)->cP() - (*fi).V(0)->cP() ) ^  ( (*fi).V(2)->cP() - (*fi).V(0)->cP() ));
          if(!pp->UseArea)
            facePlane.Normalize();
          facePlane.SetOffset( facePlane.Direction().dot((*fi).V(0)->cP()));                   

          QuadricType q;
          q.ByPlane(facePlane);          
          
          // The basic < add face quadric to each vertex > loop
          for(int j=0;j<3;++j)
            if( (*fi).V(j)->IsW() )
              QH::Qd((*fi).V(j)) += q;
          
          for(int j=0;j<3;++j)
            if( (*fi).IsB(j) || pp->QualityQuadric )
            {
              Plane3<ScalarType,false> borderPlane; 
              QuadricType bq;
              // Border quadric record the squared distance from the plane orthogonal to the face and passing 
              // through the edge. 
              borderPlane.SetDirection(facePlane.Direction() ^ ( (*fi).V1(j)->cP() - (*fi).V(j)->cP() ).normalized());
              if(  (*fi).IsB(j) ) borderPlane.SetDirection(borderPlane.Direction()* (ScalarType)(pp->BoundaryWeight ));        // amplify border planes
              else                borderPlane.SetDirection(borderPlane.Direction()* (ScalarType)(pp->BoundaryWeight/100.0));   // and consider much less quadric for quality
              borderPlane.SetOffset(borderPlane.Direction().dot((*fi).V(j)->cP()));
              bq.ByPlane(borderPlane);
              
              if( (*fi).V (j)->IsW() )	QH::Qd((*fi).V (j)) += bq;
              if( (*fi).V1(j)->IsW() )	QH::Qd((*fi).V1(j)) += bq;
            }
        }
    
    if(pp->ScaleIndependent)
    {
      vcg::tri::UpdateBounding<TriMeshType>::Box(m);
      //Make all quadric independent from mesh size
      pp->ScaleFactor = 1e8*pow(1.0/m.bbox.Diag(),6); // scaling factor
    }

    if(pp->QualityWeight) // we map quality range into a squared 01 and than this into the 1..QualityWeightFactor range
    {
      float minQ, maxQ;
      tri::Stat<TriMeshType>::ComputePerVertexQualityMinMax(m,minQ,maxQ);      
      for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
        if( ! (*vi).IsD() && (*vi).IsW())
        {
          const double quality01squared = pow((double)((vi->Q()-minQ)/(maxQ-minQ)),2.0);
          QH::Qd(*vi) *= 1.0 + quality01squared * (pp->QualityWeightFactor-1.0); 
        }
    }
  }
  


//
//
//
//
//
//
//static void InitMesh(MESH_TYPE &m){
//	pp->CosineThr=cos(pp->NormalThr);
//  InitQuadric(m);
//	//m.Topology();
//	//OldInitQuadric(m,UseArea);
//	}
//
 CoordType ComputeMinimal()
{
   typename TriMeshType::VertexType * v[2];
   v[0] = this->pos.V(0);
   v[1] = this->pos.V(1);
   QuadricType q=QH::Qd(v[0]);
   q+=QH::Qd(v[1]);
   
   Point3<QuadricType::ScalarType> x;
   
   bool rt=q.Minimum(x);
   if(!rt) { // if the computation of the minimum fails we choose between the two edge points and the middle one.
     Point3<QuadricType::ScalarType> x0=Point3d::Construct(v[0]->P());
     Point3<QuadricType::ScalarType> x1=Point3d::Construct(v[1]->P());
     x.Import((v[0]->P()+v[1]->P())/2);
     double qvx=q.Apply(x);
     double qv0=q.Apply(x0);
     double qv1=q.Apply(x1);
     if(qv0<qvx) x=x0;
     if(qv1<qvx && qv1<qv0) x=x1;
   }
   
   return CoordType::Construct(x);
 }
//
//

};
        } // namespace tri
    } // namespace vcg
#endif
