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

#ifndef __VCGLIB_PLATONIC
#define __VCGLIB_PLATONIC

#include<vcg/math/base.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/position.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/polygon_support.h>
#include<vcg/complex/algorithms/smooth.h>


namespace vcg {
namespace tri {
/** \addtogroup trimesh */
//@{
    /**
        A set of functions that builds meshes
        that represent surfaces of platonic solids,
                and other simple shapes.

                 The 1st parameter is the mesh that will
                be filled with the solid.
        */
template <class TetraMeshType>
void Tetrahedron(TetraMeshType &in)
{
 typedef typename TetraMeshType::CoordType CoordType;
 typedef typename TetraMeshType::VertexPointer  VertexPointer;
 typedef typename TetraMeshType::VertexIterator VertexIterator;
 typedef typename TetraMeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<TetraMeshType>::AddVertices(in,4);
 Allocator<TetraMeshType>::AddFaces(in,4);

 VertexPointer ivp[4];
 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1.0, 1.0, 1.0); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (-1.0, 1.0,-1.0); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (-1.0,-1.0, 1.0); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType ( 1.0,-1.0,-1.0);

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1];
}


/// builds a Dodecahedron,
/// (each pentagon is composed of 5 triangles)
template <class DodMeshType>
void Dodecahedron(DodMeshType & in)
{
 typedef DodMeshType MeshType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;
 typedef typename MeshType::ScalarType     ScalarType;
 const int N_penta=12;
 const int N_points=62;

 int penta[N_penta*3*3]=
    {20,11, 18,  18, 11,  8,  8, 11,  4,
        13,23,  4,  4, 23,  8,  8, 23, 16,
    13, 4, 30, 30,  4, 28, 28, 4,  11,
    16,34,  8,  8, 34, 18, 18, 34, 36,
    11,20, 28, 28, 20, 45, 45, 20, 38,
    13,30, 23, 23, 30, 41, 41, 30, 47,
    16,23, 34, 34, 23, 50, 50, 23, 41,
    20,18, 38, 38, 18, 52, 52, 18, 36,
    30,28, 47, 47, 28, 56, 56, 28, 45,
    50,60, 34, 34, 60, 36, 36, 60, 52,
    45,38, 56, 56, 38, 60, 60, 38, 52,
    50,41, 60, 60, 41, 56, 56, 41, 47 };
   //A B   E                D       C
  const ScalarType p=(1.0 + math::Sqrt(5.0)) / 2.0;
  const ScalarType p2=p*p;
  const ScalarType p3=p*p*p;
    ScalarType vv[N_points*3]=
    {
   0, 0, 2*p2,     p2, 0, p3,      p, p2, p3,
   0, p, p3,       -p, p2, p3,     -p2, 0, p3,
   -p, -p2, p3,    0,   -p, p3,    p,  -p2, p3,
   p3,  p, p2,     p2,  p2, p2,    0,   p3, p2,
   -p2, p2, p2,    -p3, p, p2,     -p3, -p, p2,
   -p2, -p2, p2,   0, -p3, p2,     p2, -p2, p2,
   p3, -p, p2,     p3, 0, p,       p2, p3, p,
   -p2, p3, p,     -p3, 0, p,      -p2, -p3, p,
   p2, -p3, p,     2*p2, 0, 0,     p3, p2, 0,
   p, p3, 0,       0, 2*p2, 0,     -p, p3, 0,
   -p3, p2, 0,     -2*p2, 0, 0,    -p3, -p2, 0,
   -p, -p3, 0,     0, -2*p2, 0,    p, -p3, 0,
   p3, -p2, 0,     p3, 0, -p,      p2, p3, -p,
   -p2, p3, -p,    -p3, 0, -p,     -p2, -p3, -p,
   p2, -p3, -p,    p3, p, -p2,     p2, p2, -p2,
   0, p3, -p2,     -p2, p2, -p2,   -p3, p, -p2,
   -p3, -p, -p2,   -p2, -p2, -p2,  0, -p3, -p2,
   p2, -p2, -p2,   p3, -p, -p2,    p2, 0, -p3,
   p, p2, -p3,     0, p, -p3,      -p, p2, -p3,
   -p2, 0, -p3,    -p, -p2, -p3,   0, -p, -p3,
   p, -p2, -p3,    0, 0, -2*p2
    };
    in.Clear();
    //in.face.clear();
  Allocator<DodMeshType>::AddVertices(in,20+12);
  Allocator<DodMeshType>::AddFaces(in, 5*12); // five pentagons, each made by 5 tri

  int h,i,j,m=0;

    bool used[N_points];
    for (i=0; i<N_points; i++) used[i]=false;

    int reindex[20+12 *10];
    ScalarType xx,yy,zz, sx,sy,sz;

    int order[5]={0,1,8,6,2};
    int added[12];

    VertexIterator vi=in.vert.begin();

    for (i=0; i<12; i++) {
        sx=sy=sz=0;
        for (int j=0; j<5; j++) {
            h= penta[ i*9 + order[j]  ]-1;
          xx=vv[h*3];yy=vv[h*3+1];zz=vv[h*3+2]; sx+=xx; sy+=yy; sz+=zz;
            if (!used[h]) {
                (*vi).P()=CoordType( xx, yy, zz ); vi++;
                used[h]=true;
                reindex[ h ] = m++;
            }
        }
        (*vi).P()=CoordType( sx/5.0, sy/5.0, sz/5.0 ); 	vi++;
        added[ i ] = m++;
    }

  std::vector<VertexPointer> index(in.vn);

    for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &(*vi);

  FaceIterator fi=in.face.begin();

    for (i=0; i<12; i++) {
        for (j=0; j<5; j++){
          (*fi).V(0)=index[added[i] ];
        (*fi).V(1)=index[reindex[penta[i*9 + order[j      ] ] -1 ] ];
        (*fi).V(2)=index[reindex[penta[i*9 + order[(j+1)%5] ] -1 ] ];
        if (HasPerFaceFlags(in)) {
        // tag faux edges
        (*fi).SetF(0);
        (*fi).SetF(2);
      }
          fi++;
        }
    }
}

template <class OctMeshType>
void Octahedron(OctMeshType &in)
{
 typedef OctMeshType MeshType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<OctMeshType>::AddVertices(in,6);
 Allocator<OctMeshType>::AddFaces(in,8);

 VertexPointer ivp[6];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 0, 0); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType ( 0, 1, 0); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType ( 0, 0, 1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (-1, 0, 0); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType ( 0,-1, 0); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType ( 0, 0,-1);

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1];
}

template <class IcoMeshType>
void Icosahedron(IcoMeshType &in)
{
 typedef IcoMeshType MeshType;
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

  ScalarType L=ScalarType((math::Sqrt(5.0)+1.0)/2.0);
    CoordType vv[12]={
    CoordType ( 0, L, 1),
    CoordType ( 0, L,-1),
    CoordType ( 0,-L, 1),
    CoordType ( 0,-L,-1),

    CoordType ( L, 1, 0),
    CoordType ( L,-1, 0),
    CoordType (-L, 1, 0),
    CoordType (-L,-1, 0),

    CoordType ( 1, 0, L),
    CoordType (-1, 0, L),
    CoordType ( 1, 0,-L),
    CoordType (-1, 0,-L)
    };

    int ff[20][3]={
        {1,0,4},{0,1,6},{2,3,5},{3,2,7},
        {4,5,10},{5,4,8},{6,7,9},{7,6,11},
        {8,9,2},{9,8,0},{10,11,1},{11,10,3},
        {0,8,4},{0,6,9},{1,4,10},{1,11,6},
        {2,5,8},{2,9,7},{3,10,5},{3,7,11}
    };


  in.Clear();
  Allocator<IcoMeshType>::AddVertices(in,12);
  Allocator<IcoMeshType>::AddFaces(in,20);
  VertexPointer ivp[12];

  VertexIterator vi;
  int i;
  for(i=0,vi=in.vert.begin();vi!=in.vert.end();++i,++vi){
    (*vi).P()=vv[i];
      ivp[i]=&*vi;
  }

 FaceIterator fi;
 for(i=0,fi=in.face.begin();fi!=in.face.end();++i,++fi){
        (*fi).V(0)=ivp[ff[i][0]];
        (*fi).V(1)=ivp[ff[i][1]];
        (*fi).V(2)=ivp[ff[i][2]];
    }
}

template <class MeshType>
void Hexahedron(MeshType &in)
{
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<MeshType>::AddVertices(in,8);
 Allocator<MeshType>::AddFaces(in,12);

 VertexPointer ivp[8];

 VertexIterator vi=in.vert.begin();

 ivp[7]=&*vi;(*vi).P()=CoordType (-1,-1,-1); ++vi;
 ivp[6]=&*vi;(*vi).P()=CoordType ( 1,-1,-1); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType (-1, 1,-1); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType ( 1, 1,-1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (-1,-1, 1); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType ( 1,-1, 1); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (-1, 1, 1); ++vi;
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 1, 1);

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[6];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[5];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[4];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[2];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[1];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[3];

  if (HasPerFaceFlags(in)) {
    FaceIterator fi=in.face.begin();
    for (int k=0; k<12; k++) {
      (*fi).SetF(1); fi++;
    }
  }

}

template <class MeshType>
void Square(MeshType &in)
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceIterator   FaceIterator;

  in.Clear();
  Allocator<MeshType>::AddVertices(in,4);
  Allocator<MeshType>::AddFaces(in,2);

  VertexPointer ivp[4];

  VertexIterator vi=in.vert.begin();
  ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 0, 0); ++vi;
  ivp[1]=&*vi;(*vi).P()=CoordType ( 0, 1, 0); ++vi;
  ivp[2]=&*vi;(*vi).P()=CoordType (-1, 0, 0); ++vi;
  ivp[3]=&*vi;(*vi).P()=CoordType ( 0,-1, 0);

  FaceIterator fi=in.face.begin();
  (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
  (*fi).V(0)=ivp[2];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[0];

  if (HasPerFaceFlags(in)) {
    FaceIterator fi=in.face.begin();
    for (int k=0; k<2; k++) {
      (*fi).SetF(2); fi++;
    }
  }
}

template <class MeshType>
void SphericalCap(MeshType &in, float angleRad, const int subdiv = 3 )
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexIterator VertexIterator;
  in.Clear();
  tri::Allocator<MeshType>::AddVertex(in,CoordType(0,0,0));
  for(int i=0;i<6;++i)
    tri::Allocator<MeshType>::AddVertex(in,CoordType(cos(math::ToRad(i*60.0)),sin(math::ToRad(i*60.0)),0));

  for(int i=0;i<6;++i)
    tri::Allocator<MeshType>::AddFace(in,&(in.vert[0]),&(in.vert[1+i]),&(in.vert[1+(i+1)%6]));

  tri::UpdateTopology<MeshType>::FaceFace(in);
  for(int i=0;i<subdiv;++i)
  {
    tri::Refine(in, MidPoint<MeshType>(&in));

    tri::UpdateFlags<MeshType>::FaceBorderFromFF(in);
    tri::UpdateFlags<MeshType>::VertexBorderFromFaceBorder(in);

    for(int i=0;i<in.vn;++i)
      if(in.vert[i].IsB())
        in.vert[i].P().Normalize();

    tri::UpdateSelection<MeshType>::VertexFromBorderFlag(in);
    tri::UpdateSelection<MeshType>::VertexInvert(in);
    tri::Smooth<MeshType>::VertexCoordLaplacian(in,10,true);
  }

  float angleHalfRad = angleRad /2.0f;
  float width = sin(angleHalfRad);
  tri::UpdatePosition<MeshType>::Scale(in,width);
  tri::Allocator<MeshType>::CompactEveryVector(in);
  for(VertexIterator vi=in.vert.begin(); vi!=in.vert.end();++vi)
  {
    float cosVi =  vi->P().Norm();
    float angVi = asin (cosVi);
    vi->P()[2] = cos(angVi) -  cos(angleHalfRad);
  }
}

// this function build a sphere starting from a eventually not empty mesh.
// If the mesh is not empty it is 'spherified' and used as base for the subdivision process.
// otherwise an icosahedron is used.
template <class MeshType>
void Sphere(MeshType &in, const int subdiv = 3 )
{
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;
    if(in.vn==0 && in.fn==0) Icosahedron(in);

    for(VertexIterator vi = in.vert.begin(); vi!=in.vert.end();++vi)
        vi->P().Normalize();

    for(int i = 0 ; i < subdiv; ++i)
    {
      MeshType newM;
      for(FaceIterator fi=in.face.begin();fi!=in.face.end();++fi)
      {
        CoordType me01 =  (fi->P(0)+fi->P(1))/2.0;
        CoordType me12 =  (fi->P(1)+fi->P(2))/2.0;
        CoordType me20 =  (fi->P(2)+fi->P(0))/2.0;
        tri::Allocator<MeshType>::AddFace(newM,me01,me12,me20);
        tri::Allocator<MeshType>::AddFace(newM,fi->P(0),me01,me20);
        tri::Allocator<MeshType>::AddFace(newM,fi->P(1),me12,me01);
        tri::Allocator<MeshType>::AddFace(newM,fi->P(2),me20,me12);
      }
      tri::Clean<MeshType>::RemoveDuplicateVertex(newM);
      tri::Append<MeshType,MeshType>::MeshCopy(in,newM);

        for(VertexIterator vi = in.vert.begin(); vi != in.vert.end(); ++vi)
            vi->P().Normalize();
    }
}


    /// r1 = raggio 1, r2 = raggio2, h = altezza (asse y)
template <class MeshType>
void Cone( MeshType& in,
          const typename MeshType::ScalarType r1,
          const typename MeshType::ScalarType r2,
          const typename MeshType::ScalarType h,
          const int SubDiv = 36  )
{
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

    int i,b1,b2;
  in.Clear();
  int VN,FN;
    if(r1==0 || r2==0) {
        VN=SubDiv+2;
        FN=SubDiv*2;
    }	else {
        VN=SubDiv*2+2;
        FN=SubDiv*4;
    }

  Allocator<MeshType>::AddVertices(in,VN);
  Allocator<MeshType>::AddFaces(in,FN);
    VertexPointer  *ivp = new VertexPointer[VN];

  VertexIterator vi=in.vert.begin();
  ivp[0]=&*vi;(*vi).P()=CoordType ( 0,-h/2.0,0 ); ++vi;
  ivp[1]=&*vi;(*vi).P()=CoordType ( 0, h/2.0,0 ); ++vi;

    b1 = b2 = 2;
 int cnt=2;
    if(r1!=0)
    {
        for(i=0;i<SubDiv;++i)
        {
            double a = math::ToRad(i*360.0/SubDiv);
            ivp[cnt]=&*vi; (*vi).P()= CoordType(r1*cos(a), -h/2.0, r1*sin(a)); ++vi;++cnt;
        }
        b2 += SubDiv;
    }

    if(r2!=0)
    {
        for(i=0;i<SubDiv;++i)
        {
            double a = math::ToRad(i*360.0/SubDiv);
            ivp[cnt]=&*vi; (*vi).P()= CoordType( r2*cos(a), h/2.0, r2*sin(a)); ++vi;++cnt;
        }
    }

  FaceIterator fi=in.face.begin();

  if(r1!=0) for(i=0;i<SubDiv;++i,++fi)	{
      (*fi).V(0)=ivp[0];
      (*fi).V(1)=ivp[b1+i];
      (*fi).V(2)=ivp[b1+(i+1)%SubDiv];
        }

    if(r2!=0) for(i=0;i<SubDiv;++i,++fi) {
      (*fi).V(0)=ivp[1];
      (*fi).V(2)=ivp[b2+i];
      (*fi).V(1)=ivp[b2+(i+1)%SubDiv];
        }

    if(r1==0) for(i=0;i<SubDiv;++i,++fi)
        {
      (*fi).V(0)=ivp[0];
      (*fi).V(1)=ivp[b2+i];
      (*fi).V(2)=ivp[b2+(i+1)%SubDiv];
        }
  if(r2==0)	for(i=0;i<SubDiv;++i,++fi){
      (*fi).V(0)=ivp[1];
      (*fi).V(2)=ivp[b1+i];
      (*fi).V(1)=ivp[b1+(i+1)%SubDiv];
        }

    if(r1!=0 && r2!=0)for(i=0;i<SubDiv;++i)
        {
      (*fi).V(0)=ivp[b1+i];
      (*fi).V(1)=ivp[b2+i];
      (*fi).V(2)=ivp[b2+(i+1)%SubDiv];
      ++fi;
      (*fi).V(0)=ivp[b1+i];
      (*fi).V(1)=ivp[b2+(i+1)%SubDiv];
      (*fi).V(2)=ivp[b1+(i+1)%SubDiv];
      ++fi;
        }
}


template <class MeshType >
void Box(MeshType &in, const typename MeshType::BoxType & bb )
{
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<MeshType>::AddVertices(in,8);
 VertexPointer ivp[8];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.min[1],bb.min[2]); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.min[1],bb.min[2]); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.max[1],bb.min[2]); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.max[1],bb.min[2]); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.min[1],bb.max[2]); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.min[1],bb.max[2]); ++vi;
 ivp[6]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.max[1],bb.max[2]); ++vi;
 ivp[7]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.max[1],bb.max[2]);

 Allocator<MeshType>::AddFace(in,ivp[2],ivp[1],ivp[0]);
 Allocator<MeshType>::AddFace(in,ivp[1],ivp[2],ivp[3]);
 Allocator<MeshType>::AddFace(in,ivp[4],ivp[2],ivp[0]);
 Allocator<MeshType>::AddFace(in,ivp[2],ivp[4],ivp[6]);
 Allocator<MeshType>::AddFace(in,ivp[1],ivp[4],ivp[0]);
 Allocator<MeshType>::AddFace(in,ivp[4],ivp[1],ivp[5]);
 Allocator<MeshType>::AddFace(in,ivp[6],ivp[5],ivp[7]);
 Allocator<MeshType>::AddFace(in,ivp[5],ivp[6],ivp[4]);
 Allocator<MeshType>::AddFace(in,ivp[3],ivp[6],ivp[7]);
 Allocator<MeshType>::AddFace(in,ivp[6],ivp[3],ivp[2]);
 Allocator<MeshType>::AddFace(in,ivp[5],ivp[3],ivp[7]);
 Allocator<MeshType>::AddFace(in,ivp[3],ivp[5],ivp[1]);

 if (HasPerFaceFlags(in)) {
    FaceIterator fi=in.face.begin();
    for (int k=0; k<12; k++) {
      (*fi).SetF(0); fi++;
    }
  }

}

// Torus
template <class MeshType>
void Torus(MeshType &m, float hRingRadius, float vRingRadius, int hRingDiv=24, int vRingDiv=12 )
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::ScalarType ScalarType;
  typedef Matrix44<ScalarType> Matrix44x;
  m.Clear();
  ScalarType angleStepV = (2.0f*M_PI)/vRingDiv;
  ScalarType angleStepH = (2.0f*M_PI)/hRingDiv;

  Allocator<MeshType>::AddVertices(m,(vRingDiv+1)*(hRingDiv+1));
  for(int i=0;i<hRingDiv+1;++i)
  {
    Matrix44x RotM; RotM.SetRotateRad(float(i%hRingDiv)*angleStepH,CoordType(0,0,1));
    for(int j=0;j<vRingDiv+1;++j)
    {
      CoordType p;
      p[0]= vRingRadius*cos(float(j%vRingDiv)*angleStepV) + hRingRadius;
      p[1] = 0;
      p[2]= vRingRadius*sin(float(j%vRingDiv)*angleStepV);

      m.vert[i*(vRingDiv+1)+j].P() = RotM*p;
    }
  }
  FaceGrid(m,vRingDiv+1,hRingDiv+1);
  tri::Clean<MeshType>::RemoveDuplicateVertex(m);
  tri::Allocator<MeshType>::CompactEveryVector(m);

}

/// Auxilary functions for superquadric surfaces
/// Used by SuperToroid and SuperEllipsoid
template <class ScalarType>
static  ScalarType _SQfnC(ScalarType a, ScalarType b){
  return math::Sgn(cos(a))*pow(fabs(cos(a)),b);
};
template <class ScalarType>
static ScalarType _SQfnS(ScalarType a, ScalarType b){
  return math::Sgn(sin(a))*pow(fabs(sin(a)),b);
};


/**
 * SuperToroid
 * 
 * Generate a  a supertoroid, e.g. a member of a family of doughnut-like surfaces 
 * (technically, a topological torus) whose shape is defined by mathematical formulas 
 * similar to those that define the superquadrics. 
 */
template <class MeshType>
void SuperToroid(MeshType &m, float hRingRadius, float vRingRadius, float vSquareness, float hSquareness, int hRingDiv=24, int vRingDiv=12 )
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::ScalarType ScalarType;
  m.Clear();
  ScalarType angleStepV = (2.0f*M_PI)/vRingDiv;
  ScalarType angleStepH = (2.0f*M_PI)/hRingDiv;
  
  ScalarType u,v;
  int count;
  Allocator<MeshType>::AddVertices(m,(vRingDiv+1)*(hRingDiv+1));
  for(int i=0;i<hRingDiv+1;++i)
  {
    u=float(i%hRingDiv)*angleStepH;
    count=0;
    for(int j=vRingDiv;j>=0;--j)
    {
      CoordType p;
      v=float(j%vRingDiv)*angleStepV;
      p[0]= (hRingRadius+vRingRadius*_SQfnC(u,vSquareness))*_SQfnC(v,hSquareness);;
      p[1]= (hRingRadius+vRingRadius*_SQfnC(u,vSquareness))*_SQfnS(v,hSquareness);
      p[2] = vRingRadius*_SQfnS(u,vSquareness);
      m.vert[i*(vRingDiv+1)+count].P() = p;
      count++;
    }
  }
  FaceGrid(m,vRingDiv+1,hRingDiv+1);
  tri::Clean<MeshType>::RemoveDuplicateVertex(m);
  tri::Allocator<MeshType>::CompactEveryVector(m);

}
/**
 * Generate a SuperEllipsoid eg  a solid whose horizontal sections are super-ellipses (Lamé curves)
 * with the same exponent r, and whose vertical sections through the center are super-ellipses with 
 * the same exponent t.
 */
template <class MeshType>
void SuperEllipsoid(MeshType &m, float rFeature, float sFeature, float tFeature, int hRingDiv=24, int vRingDiv=12 )
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::ScalarType ScalarType;
  m.Clear();
  ScalarType angleStepV = (2.0f*M_PI)/vRingDiv;
  ScalarType angleStepH = (1.0f*M_PI)/hRingDiv;
  float u;
  float v;
  Allocator<MeshType>::AddVertices(m,(vRingDiv+1)*(hRingDiv+1));
  for(int i=0;i<hRingDiv+1;++i)
  {
    //u=ScalarType(i%hRingDiv)*angleStepH + angleStepH/2.0;
    u=i*angleStepH;
    for(int j=0;j<vRingDiv+1;++j)
    {
      CoordType p;
      v=ScalarType(j%vRingDiv)*angleStepV;
      p[0] = _SQfnC(v,2/rFeature)*_SQfnC(u,2/rFeature);
      p[1] = _SQfnC(v,2/sFeature)*_SQfnS(u,2/sFeature);
      p[2] = _SQfnS(v,2/tFeature);
      m.vert[i*(vRingDiv+1)+j].P() = p;
    }
  }
  FaceGrid(m,vRingDiv+1,hRingDiv+1);
  tri::Clean<MeshType>::MergeCloseVertex(m,ScalarType(angleStepV*angleStepV*0.001));
  tri::Allocator<MeshType>::CompactEveryVector(m);
  bool oriented, orientable;
  tri::UpdateTopology<MeshType>::FaceFace(m);
  tri::Clean<MeshType>::OrientCoherentlyMesh(m,oriented,orientable);  
  tri::UpdateSelection<MeshType>::Clear(m);
}

/** This function build a mesh starting from a vector of generic coords (InCoordType) and indexes (InFaceIndexType)
 *  InCoordsType needs to have a [] access method for accessing the three coordinates
 *  and similarly the InFaceIndexType requires [] access method for accessing the three indexes
 */

template <class MeshType, class InCoordType, class InFaceIndexType >
void BuildMeshFromCoordVectorIndexVector(MeshType & in, const std::vector<InCoordType> & v, const std::vector<InFaceIndexType> & f)
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;

  in.Clear();
  Allocator<MeshType>::AddVertices(in,v.size());
  Allocator<MeshType>::AddFaces(in,f.size());

  for(size_t i=0;i<v.size();++i)
  {
    const InCoordType &vv = v[i];
    in.vert[i].P() = CoordType( vv[0],vv[1],vv[2]);
  }

  std::vector<VertexPointer> index(in.vn);
  VertexIterator j;
  int k;
  for(k=0,j=in.vert.begin();j!=in.vert.end();++j,++k)
    index[k] = &*j;

  for(size_t i=0;i<f.size();++i)
  {
    const InFaceIndexType &ff= f[i];
    assert( ff[0]>=0 );
    assert( ff[1]>=0 );
    assert( ff[2]>=0 );
    assert( ff[0]<in.vn );
    assert( ff[1]<in.vn );
    assert( ff[2]<in.vn );
    in.face[i].V(0) = &in.vert[ ff[0] ];
    in.face[i].V(1) = &in.vert[ ff[0] ];
    in.face[i].V(2) = &in.vert[ ff[0] ];
  }

  tri::UpdateBounding<MeshType>::Box(in);
}


template <class MeshType,class V>
void BuildMeshFromCoordVector( MeshType & in, const V & v)
{
  std::vector<Point3i> dummyfaceVec;
  BuildMeshFromCoordVectorIndexVector(in,v,dummyfaceVec);
}


template <class TriMeshType,class EdgeMeshType >
void BuildFromNonFaux(TriMeshType &in, EdgeMeshType &out)
{
  tri::RequireCompactness(in);
  std::vector<typename tri::UpdateTopology<TriMeshType>::PEdge> edgevec;
  tri::UpdateTopology<TriMeshType>::FillUniqueEdgeVector(in, edgevec, false);
  out.Clear();
  for(size_t i=0;i<in.vert.size();++i)
    tri::Allocator<EdgeMeshType>::AddVertex(out, in.vert[i].P());
  tri::UpdateFlags<EdgeMeshType>::VertexClearV(out);

  for(size_t i=0;i<edgevec.size();++i)
  {
    int i0 = tri::Index(in,edgevec[i].v[0]);
    int i1 = tri::Index(in,edgevec[i].v[1]);
    out.vert[i0].SetV();
    out.vert[i1].SetV();
    tri::Allocator<EdgeMeshType>::AddEdge(out,&out.vert[i0],&out.vert[i1]);
    if(in.vert[i0].IsS()) out.vert[i0].SetS();
    if(in.vert[i1].IsS()) out.vert[i1].SetS();
  }

  for(size_t i=0;i<out.vert.size();++i)
    if(!out.vert[i].IsV()) tri::Allocator<EdgeMeshType>::DeleteVertex(out,out.vert[i]);

  tri::Allocator<EdgeMeshType>::CompactEveryVector(out);
}

// Build a regular grid mesh as a typical height field mesh
// x y are the position on the grid scaled by wl and hl (at the end x is in the range 0..wl and y is in 0..hl)
// z is taken from the <data> array
// Once generated the vertex positions it uses the FaceGrid function to generate the faces;

template <class MeshType>
void Grid(MeshType & in, int w, int h, float wl, float hl, float *data=0)
{
  typedef typename MeshType::CoordType CoordType;

  in.Clear();
  Allocator<MeshType>::AddVertices(in,w*h);

  float wld=wl/float(w-1);
  float hld=hl/float(h-1);
  float zVal=0;
  for(int i=0;i<h;++i)
    for(int j=0;j<w;++j)
    {
      if(data) zVal=data[i*w+j];
      in.vert[i*w+j].P()=CoordType ( j*wld, i*hld, zVal) ;
    }
  FaceGrid(in,w,h);
}


// Build a regular grid mesh of faces as a typical height field mesh
// Vertexes are assumed to be already be allocated.

template <class MeshType>
void FaceGrid(MeshType & in, int w, int h)
{
    assert(in.vn == (int)in.vert.size()); // require a compact vertex vector
    assert(in.vn >= w*h); // the number of vertices should match the number of expected grid vertices

    Allocator<MeshType>::AddFaces(in,(w-1)*(h-1)*2);

//   i+0,j+0 -- i+0,j+1
//      |   \     |
//      |    \    |
//      |     \   |
//      |      \  |
//   i+1,j+0 -- i+1,j+1
//
  for(int i=0;i<h-1;++i)
    for(int j=0;j<w-1;++j)
    {
      in.face[2*(i*(w-1)+j)+0].V(0) = &(in.vert[(i+1)*w+j+1]);
      in.face[2*(i*(w-1)+j)+0].V(1) = &(in.vert[(i+0)*w+j+1]);
      in.face[2*(i*(w-1)+j)+0].V(2) = &(in.vert[(i+0)*w+j+0]);

      in.face[2*(i*(w-1)+j)+1].V(0) = &(in.vert[(i+0)*w+j+0]);
      in.face[2*(i*(w-1)+j)+1].V(1) = &(in.vert[(i+1)*w+j+0]);
      in.face[2*(i*(w-1)+j)+1].V(2) = &(in.vert[(i+1)*w+j+1]);
    }

  if (HasPerFaceFlags(in)) {
    for (int k=0; k<(h-1)*(w-1)*2; k++) {
      in.face[k].SetF(2);
    }
  }

}


// Build a regular grid mesh of faces as the resulto of a sparsely regularly sampled height field.
// Vertexes are assumed to be already be allocated, but not all the grid vertexes are present.
// For this purpose vector with a grid of indexes is also passed. 
// Negative indexes in this vector means that there is no vertex.

template <class MeshType>
void SparseFaceGrid(MeshType & in, const std::vector<int> &grid, int w, int h)
{
    tri::RequireCompactness(in);
    assert(in.vn <= w*h); // the number of vertices should match the number of expected grid vertices

//	    V0       V1
//   i+0,j+0 -- i+0,j+1
//      |   \     |
//      |    \    |
//      |     \   |
//      |      \  |
//   i+1,j+0 -- i+1,j+1
//     V2        V3


  for(int i=0;i<h-1;++i)
    for(int j=0;j<w-1;++j)
    {
      int V0i= grid[(i+0)*w+j+0];
            int V1i= grid[(i+0)*w+j+1];
            int V2i= grid[(i+1)*w+j+0];
            int V3i= grid[(i+1)*w+j+1];

            int ndone=0;
            bool quad = (V0i>=0 && V1i>=0 && V2i>=0 && V3i>=0 ) && tri::HasPerFaceFlags(in);

            if(V0i>=0 && V2i>=0 && V3i>=0 )
            {
                typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
                f->V(0)=&(in.vert[V3i]);
                f->V(1)=&(in.vert[V2i]);
                f->V(2)=&(in.vert[V0i]);
                if (quad) f->SetF(2);
                ndone++;
            }
            if(V0i>=0 && V1i>=0 && V3i>=0 )
            {
                typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
                f->V(0)=&(in.vert[V0i]);
                f->V(1)=&(in.vert[V1i]);
                f->V(2)=&(in.vert[V3i]);
                if (quad) f->SetF(2);
                ndone++;
            }

            if (ndone==0) { // try diag the other way
         if(V2i>=0 && V0i>=0 && V1i>=0 )
            {
                typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
                f->V(0)=&(in.vert[V2i]);
                f->V(1)=&(in.vert[V0i]);
                f->V(2)=&(in.vert[V1i]);
                ndone++;
             }
             if(V1i>=0 && V3i>=0 && V2i>=0 )
             {
                typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
                f->V(0)=&(in.vert[V1i]);
                f->V(1)=&(in.vert[V3i]);
                f->V(2)=&(in.vert[V2i]);
                ndone++;
             }
      }
    }
}
template <class MeshType>
void Annulus(MeshType & m, float externalRadius, float internalRadius, int slices)
{
  m.Clear();
  typename MeshType::VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertices(m,slices*2);

  for ( int j = 0; j < slices; ++j)
  {
    float x = cos(	2.0 * M_PI / slices * j);
    float y = sin(	2.0 * M_PI / slices * j);

    (*vi).P() = typename MeshType::CoordType(x,y,0)*internalRadius;
    ++vi;
    (*vi).P() = typename MeshType::CoordType(x,y,0)*externalRadius;
    ++vi;
  }
  typename MeshType::FaceIterator fi ;
  for ( int j = 0; j < slices; ++j)
  {
    fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
    (*fi).V(0) = &m.vert[ ((j+0)*2+0)%(slices*2) ];
    (*fi).V(1) = &m.vert[ ((j+1)*2+1)%(slices*2) ];
    (*fi).V(2) = &m.vert[ ((j+0)*2+1)%(slices*2) ];

    fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
    (*fi).V(0) = &m.vert[ ((j+1)*2+0)%(slices*2) ];
    (*fi).V(1) = &m.vert[ ((j+1)*2+1)%(slices*2) ];
    (*fi).V(2) = &m.vert[ ((j+0)*2+0)%(slices*2) ];
  }
}

template <class MeshType>
void OrientedAnnulus(MeshType & m, Point3f center, Point3f norm, float externalRadius, float internalRadius, int slices)
{
  Annulus(m,externalRadius,internalRadius, slices);
  float angleRad = Angle(Point3f(0,0,1),norm);
  Point3f axis = Point3f(0,0,1)^norm;

  Matrix44f rotM;
  rotM.SetRotateRad(angleRad,axis);
  tri::UpdatePosition<MeshType>::Matrix(m,rotM);
  tri::UpdatePosition<MeshType>::Translate(m,center);
}


template <class MeshType>
void Disk(MeshType & m, int slices)
{
  m.Clear();
  typename MeshType::VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertices(m,slices+1);
  (*vi).P() = typename MeshType::CoordType(0,0,0);
  ++vi;

  for ( int j = 0; j < slices; ++j)
  {
    float x = cos(	2.0 * M_PI / slices * j);
    float y = sin(	2.0 * M_PI / slices * j);

    (*vi).P() = typename MeshType::CoordType(x,y,0);
    ++vi;
  }
  typename MeshType::FaceIterator fi ;
  for ( int j = 0; j < slices; ++j)
  {
    int a =  1+(j+0)%slices;
    int b =  1+(j+1)%slices;
    fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
    (*fi).V(0) = &m.vert[ 0 ];
    (*fi).V(1) = &m.vert[ a ];
    (*fi).V(2) = &m.vert[ b ];
  }
}

template <class MeshType>
void OrientedDisk(MeshType &m, int slices, typename MeshType::CoordType center, typename MeshType::CoordType norm, float radius)
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType  CoordType;

    Disk(m,slices);
    tri::UpdatePosition<MeshType>::Scale(m,radius);
    ScalarType angleRad = Angle(CoordType(0,0,1),norm);
    CoordType axis = CoordType(0,0,1)^norm;

    Matrix44<ScalarType> rotM;
    rotM.SetRotateRad(angleRad,axis);
    tri::UpdatePosition<MeshType>::Matrix(m,rotM);
    tri::UpdatePosition<MeshType>::Translate(m,center);
}

template <class MeshType>
void OrientedEllipticPrism(MeshType & m, const typename MeshType::CoordType origin, const typename MeshType::CoordType end, float radius, float xScale, float yScale,bool capped, int slices=32, int stacks=4 )
{
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef Matrix44<typename MeshType::ScalarType> Matrix44x;
  Cylinder(slices,stacks,m,capped);
  tri::UpdatePosition<MeshType>::Translate(m,CoordType(0,1,0));
  tri::UpdatePosition<MeshType>::Scale(m,CoordType(1,0.5f,1));
  tri::UpdatePosition<MeshType>::Scale(m,CoordType(xScale,1.0f,yScale));

  float height = Distance(origin,end);
  tri::UpdatePosition<MeshType>::Scale(m,CoordType(radius,height,radius));
  CoordType norm = end-origin;
  ScalarType angleRad = Angle(CoordType(0,1,0),norm);
  CoordType axis = CoordType(0,1,0)^norm;
  Matrix44x rotM;
  rotM.SetRotateRad(angleRad,axis);
  tri::UpdatePosition<MeshType>::Matrix(m,rotM);
  tri::UpdatePosition<MeshType>::Translate(m,origin);

}

template <class MeshType>
void OrientedCylinder(MeshType & m, const typename MeshType::CoordType origin, const typename MeshType::CoordType end, float radius, bool capped, int slices=32, int stacks=4 )
{
  OrientedEllipticPrism(m,origin,end,radius,1.0f,1.0f,capped,slices,stacks);
}


template <class MeshType>
void Cylinder(int slices, int stacks, MeshType & m, bool capped=false)
{
  m.Clear();
  typename MeshType::VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertices(m,slices*(stacks+1));
  for ( int i = 0; i < stacks+1; ++i)
    for ( int j = 0; j < slices; ++j)
    {
      float x,y,h;
      x = cos(	2.0 * M_PI / slices * j);
      y = sin(	2.0 * M_PI / slices * j);
      h = 2 * i / (float)(stacks) - 1;

      (*vi).P() = typename MeshType::CoordType(x,h,y);
      ++vi;
    }

  for ( int j = 0; j < stacks; ++j)
    for ( int i = 0; i < slices; ++i)
    {
      int a,b,c,d;
      a =  (j+0)*slices + i;
      b =  (j+1)*slices + i;
      c =  (j+1)*slices + (i+1)%slices;
      d =  (j+0)*slices + (i+1)%slices;
      if(((i+j)%2) == 0){
        vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ a ], &m.vert[ b ], &m.vert[ c ]);
        vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ c ], &m.vert[ d ], &m.vert[ a ]);
      }
      else{
        vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ b ], &m.vert[ c ], &m.vert[ d ]);
        vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ d ], &m.vert[ a ], &m.vert[ b ]);
      }
    }

  if(capped)
  {
    tri::Allocator<MeshType>::AddVertex(m,typename MeshType::CoordType(0,-1,0));
    tri::Allocator<MeshType>::AddVertex(m,typename MeshType::CoordType(0, 1,0));
    int base = 0;
    for ( int i = 0; i < slices; ++i)
       vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ m.vn-2 ], &m.vert[ base+i ], &m.vert[ base+(i+1)%slices ]);
    base = (stacks)*slices;
    for ( int i = 0; i < slices; ++i)
       vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ m.vn-1 ], &m.vert[ base+(i+1)%slices ], &m.vert[ base+i ]);
  }
  if (HasPerFaceFlags(m)) {
    for (typename MeshType::FaceIterator fi=m.face.begin(); fi!=m.face.end(); fi++) {
      (*fi).SetF(2);
    }
  }
}



class _SphFace;
class _SphVertex;
struct _SphUsedTypes : public UsedTypes<	Use<_SphVertex>   ::AsVertexType,
                                        Use<_SphFace>     ::AsFaceType>{};

class _SphVertex  : public Vertex<_SphUsedTypes,  vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class _SphFace    : public Face< _SphUsedTypes,   face::VertexRef, face::Normal3f, face::BitFlags, face::FFAdj > {};
class _SphMesh    : public tri::TriMesh< vector<_SphVertex>, vector<_SphFace>   > {};


template <class MeshType>
void BuildPrismFaceShell(MeshType &mIn, MeshType &mOut, float height=0, float inset=0, bool smoothFlag=true  )
{
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::FacePointer FacePointer;
  typedef typename MeshType::CoordType CoordType;
  if(height==0) height = mIn.bbox.Diag()/100.0f;
  if(inset==0) inset = mIn.bbox.Diag()/200.0f;
  tri::UpdateTopology<MeshType>::FaceFace(mIn);
  tri::UpdateFlags<MeshType>::FaceClearV(mIn);
  for(size_t i=0;i<mIn.face.size();++i) if(!mIn.face[i].IsV())
  {
    MeshType faceM;
    std::vector<VertexPointer> vertVec;
    std::vector<FacePointer> faceVec;
    tri::PolygonSupport<MeshType,MeshType>::ExtractPolygon(&(mIn.face[i]),vertVec,faceVec);
    size_t vn = vertVec.size();

    CoordType nf(0,0,0);
    for(size_t j=0;j<faceVec.size();++j)
      nf+=faceVec[j]->N().Normalize() * DoubleArea(*faceVec[j]);
    nf.Normalize();
    nf = nf*height/2.0f;

    CoordType bary(0,0,0);
    for(size_t j=0;j<faceVec.size();++j)
      bary+= Barycenter(*faceVec[j]);
    bary/=float(faceVec.size());

    // Add vertices (alternated top and bottom)
    tri::Allocator<MeshType>::AddVertex(faceM, bary-nf);
    tri::Allocator<MeshType>::AddVertex(faceM, bary+nf);
    for(size_t j=0;j<vn;++j){
      CoordType delta = (vertVec[j]->P() - bary);
      delta.Normalize();
      delta = delta*inset;
      tri::Allocator<MeshType>::AddVertex(faceM, vertVec[j]->P()-delta-nf);
      tri::Allocator<MeshType>::AddVertex(faceM, vertVec[j]->P()-delta+nf);
    }

    // Build top and bottom faces
    for(size_t j=0;j<vn;++j)
      tri::Allocator<MeshType>::AddFace(faceM, 0, 2+(j+0)*2, 2+((j+1)%vn)*2 );
    for(size_t j=0;j<vn;++j)
      tri::Allocator<MeshType>::AddFace(faceM, 1, 3+((j+1)%vn)*2, 3+(j+0)*2 );

    // Build side strip
    for(size_t j=0;j<vn;++j){
      size_t j0=j;
      size_t j1=(j+1)%vn;
      tri::Allocator<MeshType>::AddFace(faceM, 2+ j0*2 + 0 , 2+ j0*2+1, 2+j1*2+0);
      tri::Allocator<MeshType>::AddFace(faceM, 2+ j0*2 + 1 , 2+ j1*2+1, 2+j1*2+0);
    }

    for(size_t j=0;j<2*vn;++j)
      faceM.face[j].SetS();

    if(smoothFlag)
    {
      faceM.face.EnableFFAdjacency();
      tri::UpdateTopology<MeshType>::FaceFace(faceM);
      tri::UpdateFlags<MeshType>::FaceBorderFromFF(faceM);
      tri::Refine(faceM, MidPoint<MeshType>(&faceM),0,true);
      tri::Refine(faceM, MidPoint<MeshType>(&faceM),0,true);
      tri::UpdateSelection<MeshType>::VertexFromFaceStrict(faceM);
      tri::Smooth<MeshType>::VertexCoordLaplacian(faceM,2,true,true);
    }

    tri::Append<MeshType,MeshType>::Mesh(mOut,faceM);

  } // end main loop for each face;
}


template <class MeshType>
void BuildCylinderEdgeShell(MeshType &mIn, MeshType &mOut, float radius=0, int slices=16, int stacks=1 )
{
  if(radius==0) radius = mIn.bbox.Diag()/100.0f;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  std::vector<PEdge> edgeVec;
  tri::UpdateTopology<MeshType>::FillUniqueEdgeVector(mIn,edgeVec,false);
  for(size_t i=0;i<edgeVec.size();++i)
  {
    MeshType mCyl;
    tri::OrientedCylinder(mCyl,edgeVec[i].v[0]->P(),edgeVec[i].v[1]->P(),radius,true,slices,stacks);
    tri::Append<MeshType,MeshType>::Mesh(mOut,mCyl);
  }
}

template <class MeshType>
void BuildSphereVertexShell(MeshType &mIn, MeshType &mOut, float radius=0, int recDiv=2 )
{
  if(radius==0) radius = mIn.bbox.Diag()/100.0f;
  for(size_t i=0;i<mIn.vert.size();++i)
  {
    MeshType mSph;
    tri::Sphere(mSph,recDiv);
    tri::UpdatePosition<MeshType>::Scale(mSph,radius);
    tri::UpdatePosition<MeshType>::Translate(mSph,mIn.vert[i].P());
    tri::Append<MeshType,MeshType>::Mesh(mOut,mSph);
  }
}

template <class MeshType>
void BuildCylinderVertexShell(MeshType &mIn, MeshType &mOut, float radius=0, float height=0, int slices=16, int stacks=1 )
{
  typedef typename MeshType::CoordType CoordType;
  if(radius==0) radius = mIn.bbox.Diag()/100.0f;
  if(height==0) height = mIn.bbox.Diag()/200.0f;
  for(size_t i=0;i<mIn.vert.size();++i)
  {
    CoordType p = mIn.vert[i].P();
    CoordType n = mIn.vert[i].N().Normalize();

    MeshType mCyl;
    tri::OrientedCylinder(mCyl,p-n*height,p+n*height,radius,true,slices,stacks);
    tri::Append<MeshType,MeshType>::Mesh(mOut,mCyl);
  }
}


template <class MeshType>
void GenerateCameraMesh(MeshType &in){
    typedef typename MeshType::CoordType MV;
    MV vv[52]={
        MV(-0.000122145 , -0.2 ,0.35),
        MV(0.000122145 , -0.2 ,-0.35),MV(-0.000122145 , 0.2 ,0.35),MV(0.000122145 , 0.2 ,-0.35),MV(0.999878 , -0.2 ,0.350349),MV(1.00012 , -0.2 ,-0.349651),MV(0.999878 , 0.2 ,0.350349),MV(1.00012 , 0.2 ,-0.349651),MV(1.28255 , 0.1 ,0.754205),MV(1.16539 , 0.1 ,1.03705),MV(0.88255 , 0.1 ,1.15421),
        MV(0.599707 , 0.1 ,1.03705),MV(0.48255 , 0.1 ,0.754205),MV(0.599707 , 0.1 ,0.471362),MV(0.88255 , 0.1 ,0.354205),MV(1.16539 , 0.1 ,0.471362),MV(1.28255 , -0.1 ,0.754205),MV(1.16539 , -0.1 ,1.03705),MV(0.88255 , -0.1 ,1.15421),MV(0.599707 , -0.1 ,1.03705),MV(0.48255 , -0.1 ,0.754205),
        MV(0.599707 , -0.1 ,0.471362),MV(1.16539 , -0.1 ,0.471362),MV(0.88255 , -0.1 ,0.354205),MV(3.49164e-005 , 0 ,-0.1),MV(1.74582e-005 , -0.0866025 ,-0.05),MV(-1.74582e-005 , -0.0866025 ,0.05),MV(-3.49164e-005 , 8.74228e-009 ,0.1),MV(-1.74582e-005 , 0.0866025 ,0.05),MV(1.74582e-005 , 0.0866025 ,-0.05),MV(-0.399913 , 1.99408e-022 ,-0.25014),
        MV(-0.399956 , -0.216506 ,-0.12514),MV(-0.400044 , -0.216506 ,0.12486),MV(-0.400087 , 2.18557e-008 ,0.24986),MV(-0.400044 , 0.216506 ,0.12486),MV(-0.399956 , 0.216506 ,-0.12514),MV(0.479764 , 0.1 ,0.754205),MV(0.362606 , 0.1 ,1.03705),MV(0.0797637 , 0.1 ,1.15421),MV(-0.203079 , 0.1 ,1.03705),MV(-0.320236 , 0.1 ,0.754205),
        MV(-0.203079 , 0.1 ,0.471362),MV(0.0797637 , 0.1 ,0.354205),MV(0.362606 , 0.1 ,0.471362),MV(0.479764 , -0.1 ,0.754205),MV(0.362606 , -0.1 ,1.03705),MV(0.0797637 , -0.1 ,1.15421),MV(-0.203079 , -0.1 ,1.03705),MV(-0.320236 , -0.1 ,0.754205),MV(0.0797637 , -0.1 ,0.354205),MV(0.362606 , -0.1 ,0.471362),
        MV(-0.203079 , -0.1 ,0.471362),	};
    int ff[88][3]={
        {0,2,3},
        {3,1,0},{4,5,7},{7,6,4},{0,1,5},{5,4,0},{1,3,7},{7,5,1},{3,2,6},{6,7,3},{2,0,4},
        {4,6,2},{10,9,8},{10,12,11},{10,13,12},{10,14,13},{10,15,14},{10,8,15},{8,17,16},{8,9,17},{9,18,17},
        {9,10,18},{10,19,18},{10,11,19},{11,20,19},{11,12,20},{12,21,20},{12,13,21},{13,23,21},{13,14,23},{14,22,23},
        {14,15,22},{15,16,22},{15,8,16},{23,16,17},{23,17,18},{23,18,19},{23,19,20},{23,20,21},{23,22,16},{25,27,26},
        {25,28,27},{25,29,28},{25,24,29},{24,31,30},{24,25,31},{25,32,31},{25,26,32},{26,33,32},{26,27,33},{27,34,33},
        {27,28,34},{28,35,34},{28,29,35},{29,30,35},{29,24,30},{35,30,31},{35,31,32},{35,32,33},{35,33,34},{42,37,36},
        {42,38,37},{42,39,38},{42,40,39},{42,41,40},{42,36,43},{36,45,44},{36,37,45},{37,46,45},{37,38,46},{38,47,46},
        {38,39,47},{39,48,47},{39,40,48},{40,51,48},{40,41,51},{41,49,51},{41,42,49},{42,50,49},{42,43,50},{43,44,50},
        {43,36,44},{51,44,45},{51,45,46},{51,46,47},{51,47,48},{51,49,50},{51,50,44},
    };

     in.Clear();
     Allocator<MeshType>::AddVertices(in,52);
     Allocator<MeshType>::AddFaces(in,88);

    in.vn=52;in.fn=88;
    int i,j;
    for(i=0;i<in.vn;i++)
                in.vert[i].P()=vv[i];;

    std::vector<typename MeshType::VertexPointer> index(in.vn);

    typename MeshType::VertexIterator vi;
    for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;
    for(j=0;j<in.fn;++j)
    {
        in.face[j].V(0)=index[ff[j][0]];
        in.face[j].V(1)=index[ff[j][1]];
        in.face[j].V(2)=index[ff[j][2]];
    }
}

template <class MeshType>
void OrientedRect(MeshType &square, float width, float height, Point3f c, Point3f dir=Point3f(0,0,0), float angleDeg=0,Point3f preRotTra = Point3f(0,0,0))
{
  float zeros[4]={0,0,0,0};
  square.Clear();
  Matrix44f rotM;
  tri::Grid(square,2,2,width,height,zeros);
  tri::UpdatePosition<MeshType>::Translate(square,Point3f(-width/2.0f,-height/2.0f,0.0f));
  if(angleDeg!=0){
    tri::UpdatePosition<MeshType>::Translate(square,preRotTra);
    rotM.SetRotateDeg(angleDeg,dir);
    tri::UpdatePosition<MeshType>::Matrix(square,rotM);
  }
  tri::UpdatePosition<MeshType>::Translate(square,c);
  tri::UpdateBounding<MeshType>::Box(square);
}

template <class MeshType>
void OrientedSquare(MeshType &square, float width, Point3f c, Point3f dir=Point3f(0,0,0), float angleDeg=0,Point3f preRotTra = Point3f(0,0,0))
{
  OrientedRect(square,width,width,c,dir,angleDeg,preRotTra);
}



//@}

} // End Namespace TriMesh
} // End Namespace vcg
#endif
