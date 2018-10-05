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

// This sample require gl.
#ifndef GLU_VERSIONS
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
  #include <windows.h>
#endif
#include <GL/gl.h>
#endif
#endif

#include<vcg/complex/complex.h>
#include<vcg/complex/append.h>

// input output
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

// topology computation
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include<vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/position.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/stat.h>

#include <vcg/complex/algorithms/intersection.h>
#include <vcg/complex/algorithms/refine.h>
#include <wrap/gl/glu_tessellator_cap.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
                                        Use<MyEdge>     ::AsEdgeType,
                                        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf,vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags >{};
class MyEdge    : public Edge<MyUsedTypes, edge::VertexRef,edge::BitFlags,edge::EEAdj>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};



void CapHole(MyMesh &m, MyMesh &capMesh, bool reverseFlag)
{
  capMesh.Clear();
  std::vector< std::vector<Point3f> > outlines;
  std::vector<Point3f> outline;

  tri::Allocator<MyMesh>::CompactVertexVector(m);
  tri::Allocator<MyMesh>::CompactFaceVector(m);
  tri::UpdateFlags<MyMesh>::FaceClearV(m);
  tri::UpdateFlags<MyMesh>::VertexClearV(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  int nv=0;

  for(size_t i=0;i<m.face.size();i++)
  {
    for (int j=0;j<3;j++)
    if (!m.face[i].IsV() && face::IsBorder(m.face[i],j))
    {
      MyFace* startB=&(m.face[i]);
      vcg::face::Pos<MyFace> p(startB,j);
      assert(p.IsBorder());
      do
      {
        assert(p.IsManifold());
        p.F()->SetV();
        outline.push_back(p.V()->P());
        p.NextB();
        nv++;
      }
      while(!p.F()->IsV());
      if (reverseFlag)
        std::reverse(outline.begin(),outline.end());

      outlines.push_back(outline);
      outline.clear();
    }
  }
  if (nv<2) return;
  MyMesh::VertexIterator vi=vcg::tri::Allocator<MyMesh>::AddVertices(capMesh,nv);
  for (size_t i=0;i<outlines.size();i++)
  {
    for(size_t j=0;j<outlines[i].size();++j,++vi)
      (&*vi)->P()=outlines[i][j];
  }

  std::vector<int> indices;
  glu_tesselator::tesselate(outlines, indices);

  std::vector<Point3f> points;
  glu_tesselator::unroll(outlines, points);
  MyMesh::FaceIterator fi=tri::Allocator<MyMesh>::AddFaces(capMesh,nv-2);
  for (size_t i=0; i<indices.size(); i+=3,++fi)
  {
    (*&fi)->V(0)=&capMesh.vert[ indices[i+0] ];
    (*&fi)->V(1)=&capMesh.vert[ indices[i+1] ];
    (*&fi)->V(2)=&capMesh.vert[ indices[i+2] ];
  }
  tri::Clean<MyMesh>::RemoveDuplicateVertex(capMesh);
  tri::UpdateBounding<MyMesh>::Box(capMesh);
}




bool SplitMesh(MyMesh &m,             /// The mesh that has to be splitted. It is NOT changed
               MyMesh &A, MyMesh &B,  /// The two resulting pieces, correct only if true is returned
               Plane3f plane)
{
  tri::Append<MyMesh,MyMesh>::Mesh(A,m);
  tri::UpdateQuality<MyMesh>::VertexFromPlane(A, plane);
  tri::QualityMidPointFunctor<MyMesh> slicingfunc(0.0f);
  tri::QualityEdgePredicate<MyMesh> slicingpred(0.0f);
  tri::UpdateTopology<MyMesh>::FaceFace(A);
  // The Actual Slicing
  tri::RefineE<MyMesh, tri::QualityMidPointFunctor<MyMesh>, tri::QualityEdgePredicate<MyMesh> > (A, slicingfunc, slicingpred, false);

  tri::Append<MyMesh,MyMesh>::Mesh(B,A);

  tri::UpdateSelection<MyMesh>::VertexFromQualityRange(A,-std::numeric_limits<float>::max(),0);
  tri::UpdateSelection<MyMesh>::FaceFromVertexStrict(A);
  for(MyMesh::FaceIterator fi=A.face.begin();fi!=A.face.end();++fi)
      if(!(*fi).IsD() && (*fi).IsS() ) tri::Allocator<MyMesh>::DeleteFace(A,*fi);
  tri::Clean<MyMesh>::RemoveUnreferencedVertex(A);

  tri::UpdateSelection<MyMesh>::VertexFromQualityRange(B,0,std::numeric_limits<float>::max());
  tri::UpdateSelection<MyMesh>::FaceFromVertexStrict(B);
  for(MyMesh::FaceIterator fi=B.face.begin();fi!=B.face.end();++fi)
      if(!(*fi).IsD() && (*fi).IsS() ) tri::Allocator<MyMesh>::DeleteFace(B,*fi);
  tri::Clean<MyMesh>::RemoveUnreferencedVertex(B);

  tri::UpdateTopology<MyMesh>::FaceFace(m);

  MyMesh Cap;
  CapHole(A,Cap,0);
  tri::Append<MyMesh,MyMesh>::Mesh(A,Cap);

  CapHole(B,Cap,0);
  tri::Append<MyMesh,MyMesh>::Mesh(B,Cap);

  tri::Clean<MyMesh>::RemoveDuplicateVertex(A);
  tri::Clean<MyMesh>::RemoveDuplicateVertex(B);
  return true;
}

void GetRandPlane(Box3f &bb, Plane3f &plane)
{
  Point3f planeCenter = bb.Center();
  Point3f planeDir = Point3f(-0.5f+float(rand())/RAND_MAX,-0.5f+float(rand())/RAND_MAX,-0.5f+float(rand())/RAND_MAX);
  planeDir.Normalize();

  plane.Init(planeCenter+planeDir*0.3f*bb.Diag()*float(rand())/RAND_MAX,planeDir);
}

int main( int argc, char **argv )
{
  if(argc<2)
  {
    printf("Usage trimesh_base <meshfilename.ply>\n");
    return -1;
  }

  MyMesh m, // The loaded mesh
         em, // the 2D polyline representing the section
         slice, // the planar mesh resulting from the triangulation of the above
         sliced; // the 3D mesh resulting by the actual slicing of m into two capped sub pieces

  if(tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
  {
    printf("Error reading file  %s\n",argv[1]);
    exit(0);
  }
  tri::UpdateBounding<MyMesh>::Box(m);
  printf("Input mesh  vn:%i fn:%i\n",m.VN(),m.FN());
  srand(time(0));

  Plane3f slicingPlane;
  GetRandPlane(m.bbox,slicingPlane);
  printf("slicing dir %5.2f %5.2f %5.2f\n",slicingPlane.Direction()[0],slicingPlane.Direction()[1],slicingPlane.Direction()[2]);
  vcg::IntersectionPlaneMesh<MyMesh, MyMesh, float>(m, slicingPlane, em );
  tri::Clean<MyMesh>::RemoveDuplicateVertex(em);
  vcg::tri::CapEdgeMesh(em,slice);
  printf("Slice  mesh has %i vert and %i faces\n", slice.VN(), slice.FN() );

  MyMesh A,B;
  SplitMesh(m,A,B,slicingPlane);
  tri::UpdatePosition<MyMesh>::Translate(A, slicingPlane.Direction()*m.bbox.Diag()/80.0);
  tri::UpdatePosition<MyMesh>::Translate(B,-slicingPlane.Direction()*m.bbox.Diag()/80.0);
  tri::Append<MyMesh,MyMesh>::Mesh(sliced,A);
  tri::Append<MyMesh,MyMesh>::Mesh(sliced,B);
  printf("Sliced mesh has %i vert and %i faces\n", sliced.VN(), sliced.FN() );

  tri::io::ExporterPLY<MyMesh>::Save(slice,"slice.ply",false);
  tri::io::ExporterPLY<MyMesh>::Save(sliced,"sliced.ply",false);

  return 0;
}
