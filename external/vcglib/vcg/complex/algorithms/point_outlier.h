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
#ifndef VCG_TRI_OUTLIERS__H
#define VCG_TRI_OUTLIERS__H

#include <vcg/space/index/kdtree/kdtree.h>


namespace vcg
{

namespace tri
{

template <class MeshType>
class OutlierRemoval
{
    public:

    typedef typename MeshType::ScalarType					ScalarType;
    typedef typename vcg::KdTree<ScalarType>				KdTreeType;
    typedef typename vcg::KdTree<ScalarType>::PriorityQueue	PriorityQueue;


    /**
      Compute an outlier probability value for each vertex of the mesh using the approch
      in the paper "LoOP: Local Outlier Probabilities". The outlier probability is stored in the
      vertex attribute "outlierScore". It use the input kdtree to find the kNearest of each vertex.

      "LoOP: local outlier probabilities" by 	Hans-Peter Kriegel et al.
      Proceedings of the 18th ACM conference on Information and knowledge management
    */
    static void ComputeLoOPScore(MeshType& mesh, KdTreeType& kdTree, int kNearest)
    {
      vcg::tri::RequireCompactness(mesh);
      typename MeshType::template PerVertexAttributeHandle<ScalarType> outlierScore = tri::Allocator<MeshType>:: template GetPerVertexAttribute<ScalarType>(mesh, std::string("outlierScore"));
      typename MeshType::template PerVertexAttributeHandle<ScalarType> sigma =        tri::Allocator<MeshType>:: template GetPerVertexAttribute<ScalarType>(mesh, std::string("sigma"));
      typename MeshType::template PerVertexAttributeHandle<ScalarType> plof =         tri::Allocator<MeshType>:: template GetPerVertexAttribute<ScalarType>(mesh, std::string("plof"));

#pragma omp parallel for schedule(dynamic, 10)
      for (size_t i = 0; i < mesh.vert.size(); i++)
      {
        PriorityQueue queue;
        kdTree.doQueryK(mesh.vert[i].cP(), kNearest, queue);
        ScalarType sum = 0;
        for (int j = 0; j < queue.getNofElements(); j++)
          sum += queue.getWeight(j);
        sum /= (queue.getNofElements());
        sigma[i] = sqrt(sum);
      }

      float mean = 0;
#pragma omp parallel for reduction(+: mean) schedule(dynamic, 10)
      for (size_t i = 0; i < mesh.vert.size(); i++)
      {
        PriorityQueue queue;
        kdTree.doQueryK(mesh.vert[i].cP(), kNearest, queue);
        ScalarType sum = 0;
        for (int j = 0; j < queue.getNofElements(); j++)
          sum += sigma[queue.getIndex(j)];
        sum /= (queue.getNofElements());
        plof[i] = sigma[i] / sum  - 1.0f;
        mean += plof[i] * plof[i];
      }

      mean /= mesh.vert.size();
      mean = sqrt(mean);

#pragma omp parallel for schedule(dynamic, 10)
      for (size_t i = 0; i < mesh.vert.size(); i++)
      {
        ScalarType value = plof[i] / (mean * sqrt(2.0f));
        double dem = 1.0 + 0.278393 * value;
        dem += 0.230389 * value * value;
        dem += 0.000972 * value * value * value;
        dem += 0.078108 * value * value * value * value;
        ScalarType op = max(0.0, 1.0 - 1.0 / dem);
        outlierScore[i] = op;
      }

      tri::Allocator<MeshType>::DeletePerVertexAttribute(mesh, std::string("sigma"));
      tri::Allocator<MeshType>::DeletePerVertexAttribute(mesh, std::string("plof"));
    };

    /**
    Select all the vertex of the mesh with an outlier probability above the input threshold [0.0, 1.0].
    */
    static int SelectLoOPOutliers(MeshType& mesh, KdTreeType& kdTree, int kNearest, float threshold)
    {
      ComputeLoOPScore(mesh, kdTree, kNearest);
      int count = 0;
      typename MeshType:: template PerVertexAttributeHandle<ScalarType> outlierScore = tri::Allocator<MeshType>::template GetPerVertexAttribute<ScalarType>(mesh, std::string("outlierScore"));
      for (int i = 0; i < mesh.vert.size(); i++)
      {
        if (outlierScore[i] > threshold)
        {
          mesh.vert[i].SetS();
          count++;
        }
      }
      return count;
    }



    /**
    Delete all the vertex of the mesh with an outlier probability above the input threshold [0.0, 1.0].
    */
    static int DeleteLoOPOutliers(MeshType& m, KdTreeType& kdTree, int kNearest, float threshold)
    {
      SelectLoOPOutliers(m,kdTree,kNearest,threshold);
      int ovn = m.vn;

      for(typename MeshType::VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
          if((*vi).IsS() ) tri::Allocator<MeshType>::DeleteVertex(m,*vi);
      tri::Allocator<MeshType>::CompactVertexVector(m);
      tri::Allocator<MeshType>::DeletePerVertexAttribute(m, std::string("outlierScore"));
      return m.vn - ovn;
    }
};

} // end namespace tri

} // end namespace vcg

#endif // VCG_TRI_OUTLIERS_H
