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

#ifndef KDTREE_FACE_H
#define KDTREE_FACE_H

#include <vector>
#include <vcg/space/distance3.h>

namespace vcg {

  /**
  * This class allows to create a Kd-Tree thought to perform the neighbour query using the mesh faces (closest search).
  * The class implemetantion is thread-safe.
  */
  template<class MeshType>
  class KdTreeFace
  {
  public:

    typedef typename MeshType::ScalarType Scalar;
    typedef typename MeshType::CoordType VectorType;
    typedef typename MeshType::BoxType AxisAlignedBoxType;
    typedef typename MeshType::FacePointer FacePointer;

    class Node
    {
    public:
      Scalar splitValue;
      unsigned int firstChildId : 24;
      unsigned int dim : 2;
      unsigned int leaf : 1;
      AxisAlignedBoxType aabb;
      std::vector<FacePointer> list;
    };
    typedef std::vector<Node> NodeListPointer;

  public:

    KdTreeFace(MeshType& mesh, unsigned int maxObjPerCell = 64, unsigned int maxDepth = 64) : epsilon(std::numeric_limits<Scalar>::epsilon())
    {
      targetCellSize = maxObjPerCell;
      targetMaxDepth = maxDepth;
      mNodes.resize(1);
      Node& node = mNodes.back();
      node.leaf = 0;
      node.aabb = mesh.bbox;
      node.aabb.Offset(VectorType(epsilon, epsilon, epsilon));
      for (int i = 0; i < mesh.face.size(); i++)
        node.list.push_back(&mesh.face[i]);
      numLevel = createTree(0, 1);
    };

    ~KdTreeFace()
    {

    };


    template <class ObjectMarker> FacePointer doQueryClosest(const VectorType& queryPoint, VectorType& narestPoint, Scalar& dist, ObjectMarker& marker, Scalar maxDist = std::numeric_limits<Scalar>::max())
    {
      if (maxDist < std::numeric_limits<Scalar>::max() && !mNodes[0].aabb.IsIn(queryPoint) && vcg::PointFilledBoxDistance(queryPoint, mNodes[0].aabb) >= maxDist)
      {
        dist = maxDist;
        return NULL;
      }
      std::vector<QueryNode> mNodeStack(numLevel + 1);
      mNodeStack[0].nodeId = 0;
      mNodeStack[0].sq = 0.f;
      unsigned int count = 1;

      Scalar minDist = maxDist;
      VectorType p;
      FacePointer face = NULL;
      while (count)
      {
        QueryNode& qnode = mNodeStack[count - 1];
        Node& node = mNodes[qnode.nodeId];

        if (qnode.sq < minDist)
        {
          if (node.leaf)
          {
            --count; // pop
            for (int i = 0; i < node.list.size(); i++)
            {
              if (!marker.IsMarked(node.list[i]))
              {
                marker.Mark(node.list[i]);
                Scalar tempDist = minDist;
                VectorType tempP;
                if (vcg::face::PointDistanceBase(*node.list[i], queryPoint, tempDist, tempP))
                {
                  if (tempDist < minDist)
                  {
                    minDist = tempDist;
                    p = tempP;
                    face = node.list[i];
                  }
                }
              }
            }
          }
          else
          {
            // replace the stack top by the farthest and push the closest
            float new_off = queryPoint[node.dim] - node.splitValue;
            float abs_off = abs(new_off);
            if (abs_off < minDist)
            {
              if (new_off < 0.)
              {
                mNodeStack[count].nodeId = node.firstChildId;
                qnode.nodeId = node.firstChildId + 1;
                new_off = vcg::PointFilledBoxDistance(queryPoint, mNodes[node.firstChildId + 1].aabb);
              }
              else
              {
                mNodeStack[count].nodeId = node.firstChildId + 1;
                qnode.nodeId = node.firstChildId;
                new_off = vcg::PointFilledBoxDistance(queryPoint, mNodes[node.firstChildId].aabb);
              }
              mNodeStack[count].sq = qnode.sq;
              qnode.sq = new_off;
              ++count;
            }
            else
            {
              if (new_off < 0.)
                qnode.nodeId = node.firstChildId;
              else
                qnode.nodeId = node.firstChildId + 1;
            }
          }
        }
        else
        {
          // pop
          --count;
        }
      }
      dist = minDist;
      narestPoint = p;
      return face;
    }

  protected:

    // element of the stack
    struct QueryNode
    {
      QueryNode() {}
      QueryNode(unsigned int id) : nodeId(id) {}
      unsigned int nodeId;  // id of the next node
      Scalar sq;            // distance to the next node
    };


    int createTree(unsigned int nodeId, unsigned int level)
    {
      Node& node = mNodes[nodeId];
      VectorType diag = node.aabb.max - node.aabb.min;
      unsigned int dim;
      if (diag.X() > diag.Y())
        dim = diag.X() > diag.Z() ? 0 : 2;
      else
        dim = diag.Y() > diag.Z() ? 1 : 2;

      node.splitValue = Scalar(0.5*(node.aabb.max[dim] + node.aabb.min[dim]));
      node.dim = dim;

      AxisAlignedBoxType leftBox, rightBox;
      leftBox.Add(node.aabb.min);
      rightBox.Add(node.aabb.max);
      if (node.dim == 0)
      {
        leftBox.Add(VectorType(node.splitValue, node.aabb.max[1], node.aabb.max[2]));
        rightBox.Add(VectorType(node.splitValue, node.aabb.min[1], node.aabb.min[2]));
      }
      else if (node.dim == 1)
      {
        leftBox.Add(VectorType(node.aabb.max[0], node.splitValue, node.aabb.max[2]));
        rightBox.Add(VectorType(node.aabb.min[0], node.splitValue, node.aabb.min[2]));
      }
      else if (node.dim == 2)
      {
        leftBox.Add(VectorType(node.aabb.max[0], node.aabb.max[1], node.splitValue));
        rightBox.Add(VectorType(node.aabb.min[0], node.aabb.min[1], node.splitValue));
      }
      leftBox.Offset(VectorType(epsilon, epsilon, epsilon));
      rightBox.Offset(VectorType(epsilon, epsilon, epsilon));

      node.firstChildId = mNodes.size();
      int firstChildId = node.firstChildId;
      mNodes.resize(mNodes.size() + 2);
      Node& parent = mNodes[nodeId];
      Node& leftChild = mNodes[firstChildId];
      Node& rightChild = mNodes[firstChildId + 1];
      leftChild.aabb.SetNull();
      rightChild.aabb.SetNull();

      for (int i = 0; i < parent.list.size(); i++)
      {
        unsigned int state = 0;
        FacePointer fp = parent.list[i];
        for (int j = 0; j < 3; j++)
        {
          if (fp->P(j)[dim] < parent.splitValue)
            state |= (1 << 0);
          else if (fp->P(j)[dim] > parent.splitValue)
            state |= (1 << 1);
          else
          {
            state |= (1 << 0);
            state |= (1 << 1);
          }
        }
        if (state & (1 << 0))
        {
          leftChild.list.push_back(fp);
          leftChild.aabb.Add(fp->P(0));
          leftChild.aabb.Add(fp->P(1));
          leftChild.aabb.Add(fp->P(2));
        }
        if (state & (1 << 1))
        {
          rightChild.list.push_back(fp);
          rightChild.aabb.Add(fp->P(0));
          rightChild.aabb.Add(fp->P(1));
          rightChild.aabb.Add(fp->P(2));
        }
      }
      parent.list.clear();
      leftChild.aabb.Intersect(leftBox);
      rightChild.aabb.Intersect(rightBox);

      int leftLevel, rightLevel;
      {
        if (leftChild.list.size() <= targetCellSize || level >= targetMaxDepth)
        {
          leftChild.leaf = 1;
          leftLevel = level;
        }
        else
        {
          leftChild.leaf = 0;
          leftLevel = createTree(firstChildId, level + 1);
        }
      }

      {
        Node& rightChild = mNodes[firstChildId + 1];
        if (rightChild.list.size() <= targetCellSize || level >= targetMaxDepth)
        {
          rightChild.leaf = 1;
          rightLevel = level;
        }
        else
        {
          rightChild.leaf = 0;
          rightLevel = createTree(firstChildId + 1, level + 1);
        }
      }
      if (leftLevel > rightLevel)
        return leftLevel;
      return rightLevel;
    };


  protected:

    NodeListPointer mNodes; //kd-tree nodes
    unsigned int numLevel;
    const Scalar epsilon;
    unsigned int targetCellSize;
    unsigned int targetMaxDepth;
  };

}
#endif
