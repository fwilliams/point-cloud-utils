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
#ifndef VCG_TRI_CONVEX_HULL_H
#define VCG_TRI_CONVEX_HULL_H

#include <queue>
#include <unordered_map>
#include <algorithm>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/clean.h>

namespace vcg
{

namespace tri
{

template <class InputMesh, class CHMesh>
class ConvexHull
{

public:

  typedef typename InputMesh::ScalarType		ScalarType;
  typedef typename InputMesh::CoordType		CoordType;
  typedef typename InputMesh::VertexPointer	InputVertexPointer;
  typedef typename InputMesh::VertexIterator	InputVertexIterator;
  typedef typename CHMesh::VertexIterator		CHVertexIterator;
  typedef typename CHMesh::VertexPointer		CHVertexPointer;
  typedef typename CHMesh::FaceIterator		CHFaceIterator;
  typedef typename CHMesh::FacePointer		CHFacePointer;

private:

  typedef std::pair<InputVertexPointer, ScalarType> Pair;


  // Initialize the convex hull with the biggest tetraedron created using the vertices of the input mesh
  static void InitConvexHull(InputMesh& mesh, CHMesh& convexHull)
  {
  typename CHMesh:: template PerVertexAttributeHandle<size_t> indexInputVertex = Allocator<CHMesh>::template GetPerVertexAttribute<size_t>(convexHull, std::string("indexInput"));
  InputVertexPointer v[3];
    //Find the 6 points with min/max coordinate values
    InputVertexIterator vi = mesh.vert.begin();
    std::vector<InputVertexPointer> minMax(6, &(*vi));
    for (; vi != mesh.vert.end(); vi++)
    {
      if ((*vi).P().X() < (*minMax[0]).P().X())
        minMax[0] = &(*vi);
      if ((*vi).P().Y() < (*minMax[1]).P().Y())
        minMax[1] = &(*vi);
      if ((*vi).P().Z() < (*minMax[2]).P().Z())
        minMax[2] = &(*vi);
      if ((*vi).P().X() > (*minMax[3]).P().X())
        minMax[3] = &(*vi);
      if ((*vi).P().Y() > (*minMax[4]).P().Y())
        minMax[4] = &(*vi);
      if ((*vi).P().Z() > (*minMax[5]).P().Z())
        minMax[5] = &(*vi);
    }
    //Find the farthest two points
    ScalarType maxDist = 0;
    for (int i = 0; i < 6; i++)
    {
      for (int j = i + 1; j < 6; j++)
      {
        float dist = (minMax[i]->P() - minMax[j]->P()).SquaredNorm();
        if (dist > maxDist)
        {
          maxDist = dist;
          v[0] = minMax[i];
          v[1] = minMax[j];
        }
      }
    }
    //Find the third point to create the base of the tetrahedron
    vcg::Line3<ScalarType> line(v[0]->P(), (v[0]->P() - v[1]->P()));
    maxDist = 0;
    for (vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++)
    {
      ScalarType dist = vcg::Distance(line, (*vi).P());
      if (dist > maxDist)
      {
        maxDist = dist;
        v[2] = &(*vi);
      }
    }
    //Create face in the convex hull
    CHVertexIterator chVi = vcg::tri::Allocator<CHMesh>::AddVertices(convexHull, 3);
    for (int i = 0; i < 3; i++)
    {
      (*chVi).P().Import(v[i]->P());
    v[i]->SetV();
      indexInputVertex[chVi] = vcg::tri::Index(mesh, v[i]);
      chVi++;
    }
    CHFaceIterator fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, 0, 1, 2);
    (*fi).N() = vcg::NormalizedTriangleNormal(*fi);

    //Find the fourth point to create the tetrahedron
    InputVertexPointer v4;
    float distance = 0;
    float absDist = -1;
    for (vi = mesh.vert.begin(); vi != mesh.vert.end(); vi++)
    {
      float tempDist = ((*vi).P() - (*fi).P(0)).dot((*fi).N());
      if (fabs(tempDist) > absDist)
      {
        distance = tempDist;
        v4 = &(*vi);
        absDist = fabs(distance);
      }
    }

    //Flip the previous face if the fourth point is above the face
    if (distance > 0)
    {
      (*fi).N() = -(*fi).N();
      CHVertexPointer tempV = (*fi).V(1);
      (*fi).V(1) = (*fi).V(2);
      (*fi).V(2) = tempV;
    }

    //Create the other 3 faces of the tetrahedron
    chVi = vcg::tri::Allocator<CHMesh>::AddVertices(convexHull, 1);
    (*chVi).P().Import(v4->P());
    indexInputVertex[chVi] = vcg::tri::Index(mesh, v4);
    v4->SetV();
    fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, &convexHull.vert[3], convexHull.face[0].V0(1), convexHull.face[0].V0(0));
    (*fi).N() = vcg::NormalizedTriangleNormal(*fi);
    fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, &convexHull.vert[3], convexHull.face[0].V1(1), convexHull.face[0].V1(0));
    (*fi).N() = vcg::NormalizedTriangleNormal(*fi);
    fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, &convexHull.vert[3], convexHull.face[0].V2(1), convexHull.face[0].V2(0));
    (*fi).N() = vcg::NormalizedTriangleNormal(*fi);
    vcg::tri::UpdateTopology<CHMesh>::FaceFace(convexHull);
  }


public:


  /**
    Return the convex hull of the input mesh using the Quickhull algorithm.
    For each vertex of the convex hull the algorithm stores the vertex index
    of the original mesh in attribute "indexInput".

    "The quickhull algorithm for convex hulls" by C. Bradford Barber et al.
    ACM Transactions on Mathematical Software, Volume 22 Issue 4, Dec. 1996
  */
  static bool ComputeConvexHull(InputMesh& mesh, CHMesh& convexHull)
  {
    vcg::tri::RequireFFAdjacency(convexHull);
    vcg::tri::RequirePerFaceNormal(convexHull);
    vcg::tri::Allocator<InputMesh>::CompactVertexVector(mesh);
    typename CHMesh:: template PerVertexAttributeHandle<size_t> indexInputVertex = Allocator<CHMesh>::template GetPerVertexAttribute<size_t>(convexHull, std::string("indexInput"));
    if (mesh.vert.size() < 4)
      return false;
    vcg::tri::UpdateFlags<InputMesh>::VertexClearV(mesh);
    InitConvexHull(mesh, convexHull);

    //Build list of visible vertices for each convex hull face and find the furthest vertex for each face
    std::vector<std::vector<InputVertexPointer>> listVertexPerFace(convexHull.face.size());
    std::vector<Pair> furthestVexterPerFace(convexHull.face.size(), std::make_pair((InputVertexPointer)NULL, 0.0f));
    for (size_t i = 0; i < mesh.vert.size(); i++)
    {
    if (!mesh.vert[i].IsV())
    {
      for (size_t j = 0; j < convexHull.face.size(); j++)
      {
        ScalarType dist = (mesh.vert[i].P() - convexHull.face[j].P(0)).dot(convexHull.face[j].N());
        if (dist > 0)
        {
          listVertexPerFace[j].push_back(&mesh.vert[i]);
          if (dist > furthestVexterPerFace[j].second)
          {
            furthestVexterPerFace[j].second = dist;
            furthestVexterPerFace[j].first = &mesh.vert[i];
          }
        }
      }
    }
    }

    for (size_t i = 0; i < listVertexPerFace.size(); i++)
    {
      if (listVertexPerFace[i].size() > 0)
      {
        //Find faces to remove and face on the border where to connect the new fan faces
        InputVertexPointer vertex = furthestVexterPerFace[i].first;
        std::queue<int> queue;
        std::vector<int> visFace;
        std::vector<int> borderFace;
        visFace.push_back(i);
        queue.push(i);
        while (queue.size() > 0)
        {
          CHFacePointer fp = &convexHull.face[queue.front()];
          queue.pop();
          fp->SetV();
          for (int ii = 0; ii < 3; ii++)
          {
            CHFacePointer nextF = fp->FFp(ii);
            if (!nextF->IsV())
            {
              int indexF = vcg::tri::Index(convexHull, nextF);
              ScalarType dist = (vertex->P() - nextF->P(0)).dot(nextF->N());
              if (dist < 0)
              {
                borderFace.push_back(indexF);
                fp->SetB(ii);
                nextF->SetB(fp->FFi(ii));
              }
              else
              {
                visFace.push_back(indexF);
                queue.push(indexF);
              }
            }
          }
        }
        if (borderFace.size() > 0)
        {
          CHVertexIterator vi = vcg::tri::Allocator<CHMesh>::AddVertices(convexHull, 1);
          (*vi).P().Import((*vertex).P());
      vertex->SetV();
          indexInputVertex[vi] = vcg::tri::Index(mesh, vertex);
        }

        //Add a new face for each border
        std::unordered_map< CHVertexPointer, std::pair<int, char> > fanMap;
        for (size_t jj = 0; jj < borderFace.size(); jj++)
        {
          int indexFace = borderFace[jj];
          CHFacePointer f = &convexHull.face[indexFace];
          for (int j = 0; j < 3; j++)
          {
            if (f->IsB(j))
            {
              f->ClearB(j);
              //Add new face
              CHFaceIterator fi = vcg::tri::Allocator<CHMesh>::AddFace(convexHull, &convexHull.vert.back(), f->V1(j), f->V0(j));
              (*fi).N() = vcg::NormalizedTriangleNormal(*fi);
              f = &convexHull.face[indexFace];
              int newFace = vcg::tri::Index(convexHull, *fi);
              //Update convex hull FF topology
              CHVertexPointer vp[] = { f->V1(j), f->V0(j) };
              for (int ii = 0; ii < 2; ii++)
              {
                int indexE = ii * 2;
                typename std::unordered_map< CHVertexPointer, std::pair<int, char> >::iterator vIter = fanMap.find(vp[ii]);
                if (vIter != fanMap.end())
                {
                  CHFacePointer f2 = &convexHull.face[(*vIter).second.first];
                  char edgeIndex = (*vIter).second.second;
                  f2->FFp(edgeIndex) = &convexHull.face.back();
                  f2->FFi(edgeIndex) = indexE;
                  fi->FFp(indexE) = f2;
                  fi->FFi(indexE) = edgeIndex;
                }
                else
                {
                  fanMap[vp[ii]] = std::make_pair(newFace, indexE);
                }
              }
              //Build the visibility list for the new face
              std::vector<InputVertexPointer> tempVect;
              int indices[2] = { indexFace, int(vcg::tri::Index(convexHull, f->FFp(j)) )};
              std::vector<InputVertexPointer> vertexToTest(listVertexPerFace[indices[0]].size() + listVertexPerFace[indices[1]].size());
              typename std::vector<InputVertexPointer>::iterator tempIt = std::set_union(listVertexPerFace[indices[0]].begin(), listVertexPerFace[indices[0]].end(), listVertexPerFace[indices[1]].begin(), listVertexPerFace[indices[1]].end(), vertexToTest.begin());
              vertexToTest.resize(tempIt - vertexToTest.begin());

              Pair newInfo = std::make_pair((InputVertexPointer)NULL , 0.0f);
              for (size_t ii = 0; ii < vertexToTest.size(); ii++)
              {
                if (!(*vertexToTest[ii]).IsV())
                {
                  float dist = ((*vertexToTest[ii]).P() - (*fi).P(0)).dot((*fi).N());
                  if (dist > 0)
                  {
                    tempVect.push_back(vertexToTest[ii]);
                    if (dist > newInfo.second)
                    {
                      newInfo.second = dist;
                      newInfo.first = vertexToTest[ii];
                    }
                  }
                }
              }
              listVertexPerFace.push_back(tempVect);
              furthestVexterPerFace.push_back(newInfo);
              //Update topology of the new face
              CHFacePointer ffp = f->FFp(j);
              int ffi = f->FFi(j);
              ffp->FFp(ffi) = ffp;
              ffp->FFi(ffi) = ffi;
              f->FFp(j) = &convexHull.face.back();
              f->FFi(j) = 1;
              fi->FFp(1) = f;
              fi->FFi(1) = j;
            }
          }
        }
        //Delete the faces inside the updated convex hull
        for (size_t j = 0; j < visFace.size(); j++)
        {
          if (!convexHull.face[visFace[j]].IsD())
          {
            std::vector<InputVertexPointer> emptyVec;
            vcg::tri::Allocator<CHMesh>::DeleteFace(convexHull, convexHull.face[visFace[j]]);
            listVertexPerFace[visFace[j]].swap(emptyVec);
          }
        }
      }
    }

  tri::UpdateTopology<CHMesh>::ClearFaceFace(convexHull);
    vcg::tri::Allocator<CHMesh>::CompactFaceVector(convexHull);
    vcg::tri::Clean<CHMesh>::RemoveUnreferencedVertex(convexHull);
    return true;
  }
  /**
  * @brief ComputePointVisibility
  * Select the <b>visible points</b> in a point cloud, as viewed from a given viewpoint.
  * It uses the Qhull implementation of che convex hull in the vcglibrary
  * The algorithm used (Katz, Tal and Basri 2007) determines visibility without
  * reconstructing a surface or estimating normals.
  * A point is considered visible if its transformed point lies on the convex hull
  * of a trasformed points cloud from the original mesh points.
  *
  * @param m         The point cloud
  * @param visible   The mesh that will contain the visible hull
  * @param viewpoint
  * @param logR      Bounds the radius of the sphere used to select visible points.
  *   It is used to adjust the radius of the sphere (calculated as distance between
  *   the center and the farthest point from it) according to the following equation:
  *       radius = radius * pow(10,threshold);
  *   As the radius increases more points are marked as visible.
  *   Use a big threshold for dense point clouds, a small one for sparse clouds.
  */

 static void ComputePointVisibility(InputMesh& m, CHMesh& visible, CoordType viewpoint, ScalarType logR=2)
 {
   visible.Clear();
   tri::RequireCompactness(m);
   InputMesh flipM;

   printf("Input mesh m %i %i\n",m.vn,m.fn);

   tri::Allocator<InputMesh>::AddVertices(flipM,m.vn);
   ScalarType maxDist=0;
   InputVertexIterator ci=flipM.vert.begin();
   for(InputVertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
   {
     ci->P()=vi->P()-viewpoint;
     maxDist = std::max(maxDist,Norm(ci->P()));
     ++ci;
   }
   ScalarType R = maxDist*pow(10,logR);
   printf("Using R = %f logR = %f maxdist=%f \n",R,logR,maxDist);
   for(InputVertexIterator vi=flipM.vert.begin();vi!=flipM.vert.end();++vi)
   {
     ScalarType d = Norm(vi->P());
     vi->P() = vi->P() + vi->P()*ScalarType(2.0*(R - d)/d);
   }

   tri::Allocator<InputMesh>::AddVertex(flipM,CoordType(0,0,0));
   assert(m.vn+1 == flipM.vn);

   ComputeConvexHull(flipM,visible);
   assert(flipM.vert[m.vn].P()==Point3f(0,0,0));
   int vpInd=-1; // Index of the viewpoint in the ConvexHull mesh
   int selCnt=0;
   typename CHMesh:: template PerVertexAttributeHandle<size_t> indexInputVertex = Allocator<InputMesh>::template GetPerVertexAttribute<size_t>(visible, std::string("indexInput"));
   for(int i=0;i<visible.vn;++i)
   {
     size_t ind = indexInputVertex[i];
     if(ind==m.vn) vpInd = i;
     else
     {
       visible.vert[i].P() = m.vert[ind].P();
       m.vert[ind].SetS();
       m.vert[ind].C() = Color4b::LightBlue;
       selCnt++;
     }
   }
   printf("Selected %i visible points\n",selCnt);

   assert(vpInd != -1);
   // Final pass delete all the faces of the convex hull incident in the viewpoint
   for(int i=0;i<visible.fn;++i)
   {
     if( (Index(visible,visible.face[i].V(0)) == vpInd) ||
         (Index(visible,visible.face[i].V(1)) == vpInd) ||
         (Index(visible,visible.face[i].V(2)) == vpInd) )
       tri::Allocator<CHMesh>::DeleteFace(visible,visible.face[i]);
   }

   tri::Allocator<CHMesh>::CompactEveryVector(visible);
   tri::Clean<CHMesh>::FlipMesh(visible);

 }
};

} // end namespace tri

} // end namespace vcg

#endif //VCG_TRI_CONVEX_HULL_H
