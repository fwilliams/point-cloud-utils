#ifndef CLOTH_MESH_H
#define CLOTH_MESH_H

#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <vector>

using namespace vcg;
          
class CVertex;
class CEdge;
class CFace;

class CEdge {
  public:
   CVertex *v[2];         
   CFace *f;
   bool operator<(const CEdge& t) const {
      if(v[0] < t.v[0]) return true;
      if(v[0] > t.v[0]) return false;
      return v[1] < t.v[1];
   }
   bool operator==(const CEdge& t) const {
      return v[0] == t.v[0] && v[1] == t.v[1];
   }
};

class CVertex: public
      VertexSimp2<CVertex, CEdge, CFace,
                  vcg::vert::Coord3f, vert::Normal3f, vert::BitFlags, vert::Mark,
                  vert::VFAdj, 
                  vert::Qualityf> {
  public:     
   float color;
};

class CFace: public FaceSimp2 <CVertex, CEdge, CFace, face::VertexRef, 
                               face::BitFlags, face::VFAdj, face::FFAdj > {};
                       
class CMesh: public tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

#endif
