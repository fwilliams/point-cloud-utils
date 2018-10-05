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

#ifndef KDTREE_VCG_H
#define KDTREE_VCG_H

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/index/kdtree/priorityqueue.h>

#include <vector>
#include <limits>
#include <iostream>
#include <cstdint>

namespace vcg {

  template<typename _DataType>
  class ConstDataWrapper
  {
  public:
    typedef _DataType DataType;
    inline ConstDataWrapper()
      : mpData(0), mStride(0), mSize(0)
    {}
    inline ConstDataWrapper(const DataType* pData, int size, int64_t stride = sizeof(DataType))
      : mpData(reinterpret_cast<const unsigned char*>(pData)), mStride(stride), mSize(size)
    {}
    inline const DataType& operator[] (int i) const
    {
      return *reinterpret_cast<const DataType*>(mpData + i*mStride);
    }
    inline size_t size() const { return mSize; }
  protected:
    const unsigned char* mpData;
    int64_t mStride;
    size_t mSize;
  };

  template<class StdVectorType>
  class VectorConstDataWrapper :public  ConstDataWrapper<typename StdVectorType::value_type>
  {
  public:
    inline VectorConstDataWrapper(StdVectorType &vec) :
      ConstDataWrapper<typename StdVectorType::value_type>(&(vec[0]), vec.size(), sizeof(typename StdVectorType::value_type))
    {}
  };

  template<class MeshType>
  class VertexConstDataWrapper :public  ConstDataWrapper<typename MeshType::CoordType>
  {
  public:
    inline VertexConstDataWrapper(MeshType &m) :
      ConstDataWrapper<typename MeshType::CoordType>(&(m.vert[0].P()), m.vert.size(), sizeof(typename MeshType::VertexType))
    {}
  };

  /**
  * This class allows to create a Kd-Tree thought to perform the neighbour query (radius search, knn-nearest serach and closest search).
  * The class implemetantion is thread-safe.
  */
  template<typename _Scalar>
  class KdTree
  {
  public:

    typedef _Scalar Scalar;
    typedef vcg::Point3<Scalar> VectorType;
    typedef vcg::Box3<Scalar> AxisAlignedBoxType;

    typedef HeapMaxPriorityQueue<int, Scalar> PriorityQueue;

    struct Node
    {
      union {
        //standard node
        struct {
          Scalar splitValue;
          unsigned int firstChildId : 24;
          unsigned int dim : 2;
          unsigned int leaf : 1;
        };
        //leaf
        struct {
          unsigned int start;
          unsigned short size;
        };
      };
    };
    typedef std::vector<Node> NodeList;

    // return the protected members which store the nodes and the points list
    inline const NodeList& _getNodes(void) { return mNodes; }
    inline const std::vector<VectorType>& _getPoints(void) { return mPoints; }
    inline unsigned int _getNumLevel(void) { return numLevel; }
    inline const AxisAlignedBoxType& _getAABBox(void) { return mAABB; }

  public:

    KdTree(const ConstDataWrapper<VectorType>& points, unsigned int nofPointsPerCell = 16, unsigned int maxDepth = 64, bool balanced = false);

    ~KdTree();

    void doQueryK(const VectorType& queryPoint, int k, PriorityQueue& mNeighborQueue);

    void doQueryDist(const VectorType& queryPoint, float dist, std::vector<unsigned int>& points, std::vector<Scalar>& sqrareDists);

    void doQueryClosest(const VectorType& queryPoint, unsigned int& index, Scalar& dist);

  protected:

    // element of the stack
    struct QueryNode
    {
      QueryNode() {}
      QueryNode(unsigned int id) : nodeId(id) {}
      unsigned int nodeId;  // id of the next node
      Scalar sq;            // squared distance to the next node
    };

    // used to build the tree: split the subset [start..end[ according to dim and splitValue,
    // and returns the index of the first element of the second subset
    unsigned int split(int start, int end, unsigned int dim, float splitValue);

    int createTree(unsigned int nodeId, unsigned int start, unsigned int end, unsigned int level);

  protected:

    AxisAlignedBoxType mAABB; //BoundingBox
    NodeList mNodes; //kd-tree nodes
    std::vector<VectorType> mPoints; //points read from the input DataWrapper
    std::vector<unsigned int> mIndices; //points indices
    unsigned int targetCellSize; //min number of point in a leaf
    unsigned int targetMaxDepth; //max tree depth
    unsigned int numLevel; //actual tree depth
    bool isBalanced; //true if the tree is balanced
  };


  template<typename Scalar>
  KdTree<Scalar>::KdTree(const ConstDataWrapper<VectorType>& points, unsigned int nofPointsPerCell, unsigned int maxDepth, bool balanced)
    : mPoints(points.size()), mIndices(points.size())
  {
    // compute the AABB of the input
    mPoints[0] = points[0];
    mAABB.Set(mPoints[0]);
    for (unsigned int i = 1; i < mPoints.size(); ++i)
    {
      mPoints[i] = points[i];
      mIndices[i] = i;
      mAABB.Add(mPoints[i]);
    }

    targetMaxDepth = maxDepth;
    targetCellSize = nofPointsPerCell;
    isBalanced = balanced;
    //mNodes.reserve(4 * mPoints.size() / nofPointsPerCell);
    //first node inserted (no leaf). The others are made by the createTree function (recursively)
    mNodes.resize(1);
    mNodes.back().leaf = 0;
    numLevel = createTree(0, 0, mPoints.size(), 1);
  }

  template<typename Scalar>
  KdTree<Scalar>::~KdTree()
  {
  }


  /** Performs the kNN query.
  *
  * This algorithm uses the simple distance to the split plane to prune nodes.
  * A more elaborated approach consists to track the closest corner of the cell
  * relatively to the current query point. This strategy allows to save about 5%
  * of the leaves. However, in practice the slight overhead due to this tracking
  * reduces the overall performance.
  *
  * This algorithm also use a simple stack while a priority queue using the squared
  * distances to the cells as a priority values allows to save about 10% of the leaves.
  * But, again, priority queue insertions and deletions are quite involved, and therefore
  * a simple stack is by far much faster.
  *
  * The result of the query, the k-nearest neighbors, are stored into the stack mNeighborQueue, where the
  * topmost element [0] is NOT the nearest but the farthest!! (they are not sorted but arranged into a heap).
  */
  template<typename Scalar>
  void KdTree<Scalar>::doQueryK(const VectorType& queryPoint, int k, PriorityQueue& mNeighborQueue)
  {
    mNeighborQueue.setMaxSize(k);
    mNeighborQueue.init();

    std::vector<QueryNode> mNodeStack(numLevel + 1);
    mNodeStack[0].nodeId = 0;
    mNodeStack[0].sq = 0.f;
    unsigned int count = 1;

    while (count)
    {
      //we select the last node (AABB) inserted in the stack
      QueryNode& qnode = mNodeStack[count - 1];

      //while going down the tree qnode.nodeId is the nearest sub-tree, otherwise,
      //in backtracking, qnode.nodeId is the other sub-tree that will be visited iff
      //the actual nearest node is further than the split distance.
      Node& node = mNodes[qnode.nodeId];

      //if the distance is less than the top of the max-heap, it could be one of the k-nearest neighbours
      if (mNeighborQueue.getNofElements() < k || qnode.sq < mNeighborQueue.getTopWeight())
      {
        //when we arrive to a leaf
        if (node.leaf)
        {
          --count; //pop of the leaf

          //end is the index of the last element of the leaf in mPoints
          unsigned int end = node.start + node.size;
          //adding the element of the leaf to the heap
          for (unsigned int i = node.start; i < end; ++i)
            mNeighborQueue.insert(mIndices[i], vcg::SquaredNorm(queryPoint - mPoints[i]));
        }
        //otherwise, if we're not on a leaf
        else
        {
          // the new offset is the distance between the searched point and the actual split coordinate
          float new_off = queryPoint[node.dim] - node.splitValue;

          //left sub-tree
          if (new_off < 0.)
          {
            mNodeStack[count].nodeId = node.firstChildId;
            //in the father's nodeId we save the index of the other sub-tree (for backtracking)
            qnode.nodeId = node.firstChildId + 1;
          }
          //right sub-tree (same as above)
          else
          {
            mNodeStack[count].nodeId = node.firstChildId + 1;
            qnode.nodeId = node.firstChildId;
          }
          //distance is inherited from the father (while descending the tree it's equal to 0)
          mNodeStack[count].sq = qnode.sq;
          //distance of the father is the squared distance from the split plane
          qnode.sq = new_off*new_off;
          ++count;
        }
      }
      else
      {
        // pop
        --count;
      }
    }
  }


  /** Performs the distance query.
  *
  * The result of the query, all the points within the distance dist form the query point, is the vector of the indeces
  * and the vector of the squared distances from the query point.
  */
  template<typename Scalar>
  void KdTree<Scalar>::doQueryDist(const VectorType& queryPoint, float dist, std::vector<unsigned int>& points, std::vector<Scalar>& sqrareDists)
  {
    std::vector<QueryNode> mNodeStack(numLevel + 1);
    mNodeStack[0].nodeId = 0;
    mNodeStack[0].sq = 0.f;
    unsigned int count = 1;

    float sqrareDist = dist*dist;
    while (count)
    {
      QueryNode& qnode = mNodeStack[count - 1];
      Node   & node = mNodes[qnode.nodeId];

      if (qnode.sq < sqrareDist)
      {
        if (node.leaf)
        {
          --count; // pop
          unsigned int end = node.start + node.size;
          for (unsigned int i = node.start; i < end; ++i)
          {
            float pointSquareDist = vcg::SquaredNorm(queryPoint - mPoints[i]);
            if (pointSquareDist < sqrareDist)
            {
              points.push_back(mIndices[i]);
              sqrareDists.push_back(pointSquareDist);
            }
          }
        }
        else
        {
          // replace the stack top by the farthest and push the closest
          float new_off = queryPoint[node.dim] - node.splitValue;
          if (new_off < 0.)
          {
            mNodeStack[count].nodeId = node.firstChildId;
            qnode.nodeId = node.firstChildId + 1;
          }
          else
          {
            mNodeStack[count].nodeId = node.firstChildId + 1;
            qnode.nodeId = node.firstChildId;
          }
          mNodeStack[count].sq = qnode.sq;
          qnode.sq = new_off*new_off;
          ++count;
        }
      }
      else
      {
        // pop
        --count;
      }
    }
  }


  /** Searchs the closest point.
  *
  * The result of the query, the closest point to the query point, is the index of the point and
  * and the squared distance from the query point.
  */
  template<typename Scalar>
  void KdTree<Scalar>::doQueryClosest(const VectorType& queryPoint, unsigned int& index, Scalar& dist)
  {
    std::vector<QueryNode> mNodeStack(numLevel + 1);
    mNodeStack[0].nodeId = 0;
    mNodeStack[0].sq = 0.f;
    unsigned int count = 1;

    int minIndex = mIndices.size() / 2;
    Scalar minDist = vcg::SquaredNorm(queryPoint - mPoints[minIndex]);
    minIndex = mIndices[minIndex];

    while (count)
    {
      QueryNode& qnode = mNodeStack[count - 1];
      Node   & node = mNodes[qnode.nodeId];

      if (qnode.sq < minDist)
      {
        if (node.leaf)
        {
          --count; // pop
          unsigned int end = node.start + node.size;
          for (unsigned int i = node.start; i < end; ++i)
          {
            float pointSquareDist = vcg::SquaredNorm(queryPoint - mPoints[i]);
            if (pointSquareDist < minDist)
            {
              minDist = pointSquareDist;
              minIndex = mIndices[i];
            }
          }
        }
        else
        {
          // replace the stack top by the farthest and push the closest
          float new_off = queryPoint[node.dim] - node.splitValue;
          if (new_off < 0.)
          {
            mNodeStack[count].nodeId = node.firstChildId;
            qnode.nodeId = node.firstChildId + 1;
          }
          else
          {
            mNodeStack[count].nodeId = node.firstChildId + 1;
            qnode.nodeId = node.firstChildId;
          }
          mNodeStack[count].sq = qnode.sq;
          qnode.sq = new_off*new_off;
          ++count;
        }
      }
      else
      {
        // pop
        --count;
      }
    }
    index = minIndex;
    dist = minDist;
  }



  /**
  * Split the subarray between start and end in two part, one with the elements less than splitValue,
  * the other with the elements greater or equal than splitValue. The elements are compared
  * using the "dim" coordinate [0 = x, 1 = y, 2 = z].
  */
  template<typename Scalar>
  unsigned int KdTree<Scalar>::split(int start, int end, unsigned int dim, float splitValue)
  {
    int l(start), r(end - 1);
    for (; l < r; ++l, --r)
    {
      while (l < end && mPoints[l][dim] < splitValue)
        l++;
      while (r >= start && mPoints[r][dim] >= splitValue)
        r--;
      if (l > r)
        break;
      std::swap(mPoints[l], mPoints[r]);
      std::swap(mIndices[l], mIndices[r]);
    }
    //returns the index of the first element on the second part
    return (mPoints[l][dim] < splitValue ? l + 1 : l);
  }

  /** recursively builds the kdtree
  *
  *  The heuristic is the following:
  *   - if the number of points in the node is lower than targetCellsize then make a leaf
  *   - else compute the AABB of the points of the node and split it at the middle of
  *     the largest AABB dimension.
  *
  *  This strategy might look not optimal because it does not explicitly prune empty space,
  *  unlike more advanced SAH-like techniques used for RT. On the other hand it leads to a shorter tree,
  *  faster to traverse and our experience shown that in the special case of kNN queries,
  *  this strategy is indeed more efficient (and much faster to build). Moreover, for volume data
  *  (e.g., fluid simulation) pruning the empty space is useless.
  *
  *  Actually, storing at each node the exact AABB (we therefore have a binary BVH) allows
  *  to prune only about 10% of the leaves, but the overhead of this pruning (ball/ABBB intersection)
  *  is more expensive than the gain it provides and the memory consumption is x4 higher !
  */
  template<typename Scalar>
  int KdTree<Scalar>::createTree(unsigned int nodeId, unsigned int start, unsigned int end, unsigned int level)
  {
    //select the first node
    Node& node = mNodes[nodeId];
    AxisAlignedBoxType aabb;

    //putting all the points in the bounding box
    aabb.Set(mPoints[start]);
    for (unsigned int i = start + 1; i < end; ++i)
      aabb.Add(mPoints[i]);

    //bounding box diagonal
    VectorType diag = aabb.max - aabb.min;

    //the split "dim" is the dimension of the box with the biggest value
    unsigned int dim;
    if (diag.X() > diag.Y())
      dim = diag.X() > diag.Z() ? 0 : 2;
    else
      dim = diag.Y() > diag.Z() ? 1 : 2;

    node.dim = dim;
    if (isBalanced) //we divide the points using the median value along the "dim" dimension
    {
      std::vector<Scalar> tempVector;
      for (unsigned int i = start + 1; i < end; ++i)
        tempVector.push_back(mPoints[i][dim]);
      std::sort(tempVector.begin(), tempVector.end());
      node.splitValue = (tempVector[tempVector.size() / 2.0] + tempVector[tempVector.size() / 2.0 + 1]) / 2.0;
    }
    else //we divide the bounding box in 2 partitions, considering the average of the "dim" dimension
      node.splitValue = Scalar(0.5*(aabb.max[dim] + aabb.min[dim]));

    //midId is the index of the first element in the second partition
    unsigned int midId = split(start, end, dim, node.splitValue);

    node.firstChildId = mNodes.size();
    mNodes.resize(mNodes.size() + 2);
    bool flag = (midId == start) || (midId == end);
    int leftLevel, rightLevel;
    {
      // left child
      unsigned int childId = mNodes[nodeId].firstChildId;
      Node& child = mNodes[childId];
      if (flag || (midId - start) <= targetCellSize || level >= targetMaxDepth)
      {
        child.leaf = 1;
        child.start = start;
        child.size = midId - start;
        leftLevel = level;
      }
      else
      {
        child.leaf = 0;
        leftLevel = createTree(childId, start, midId, level + 1);
      }
    }

    {
      // right child
      unsigned int childId = mNodes[nodeId].firstChildId + 1;
      Node& child = mNodes[childId];
      if (flag || (end - midId) <= targetCellSize || level >= targetMaxDepth)
      {
        child.leaf = 1;
        child.start = midId;
        child.size = end - midId;
        rightLevel = level;
      }
      else
      {
        child.leaf = 0;
        rightLevel = createTree(childId, midId, end, level + 1);
      }
    }
    if (leftLevel > rightLevel)
      return leftLevel;
    return rightLevel;
  }

}

#endif
