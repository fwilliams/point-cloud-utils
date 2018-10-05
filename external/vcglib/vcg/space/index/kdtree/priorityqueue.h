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

#ifndef _PRIORITYQUEUE_H_
#define _PRIORITYQUEUE_H_

#include <algorithm>

namespace vcg {

  /** Implements a bounded-size max priority queue using a heap
  */
  template <typename Index, typename Weight>
  class HeapMaxPriorityQueue
  {
    struct Element
    {
      Weight weight;
      Index index;
    };


    struct
    {
      bool operator()(const Element& a, const Element& b) const
      {
        return a.weight < b.weight;
      }
    } lessElement;


    struct
    {
      bool operator()(const Element& a, const Element& b) const
      {
        return a.weight > b.weight;
      }
    } greaterElement;



  public:

    HeapMaxPriorityQueue(void)
    {
      mElements = 0;
      mMaxSize = 0;
    }

    ~HeapMaxPriorityQueue()
    {
      if (mElements)
        delete[] mElements;
    }


    inline void setMaxSize(int maxSize)
    {
      if (mMaxSize!=maxSize)
      {
        mMaxSize = maxSize;
        delete[] mElements;
        mElements = new Element[mMaxSize];
        mpOffsetedElements = (mElements-1);
      }
      init();
    }

    inline void init() { mCount = 0; }

    inline bool isFull() const { return mCount == mMaxSize; }

    /** returns number of elements inserted in queue
    */
    inline int getNofElements() const { return mCount; }

    inline Weight getWeight(int i) const { return mElements[i].weight; }
    inline Index  getIndex(int i) const { return mElements[i].index; }

    inline Weight getTopWeight() const { return mElements[0].weight; }

    inline void insert(Index index, Weight weight)
    {
      if (mCount==mMaxSize)
      {
        if (weight<mElements[0].weight)
        {
          int j, k;
          j = 1;
          k = 2;
          while (k <= mMaxSize)
          {
            Element* z = &(mpOffsetedElements[k]);
            if ((k < mMaxSize) && (z->weight < mpOffsetedElements[k+1].weight))
              z = &(mpOffsetedElements[++k]);

            if(weight >= z->weight)
              break;
            mpOffsetedElements[j] = *z;
            j = k;
            k = 2 * j;
          }
          mpOffsetedElements[j].weight = weight;
          mpOffsetedElements[j].index = index;
        }
      }
      else
      {
        int i, j;
        i = ++mCount;
        while (i >= 2)
        {
          j = i >> 1;
          Element& y = mpOffsetedElements[j];
          if(weight <= y.weight)
            break;
          mpOffsetedElements[i] = y;
          i = j;
        }
        mpOffsetedElements[i].index = index;
        mpOffsetedElements[i].weight = weight;
      }
    }

    inline void sort(bool ascending = true)
    {
      if (ascending)
        std::sort(mElements, mElements + mCount, lessElement);
      else
        std::sort(mElements, mElements + mCount, greaterElement);
    }

  protected:

    int mCount;
    int mMaxSize;
    Element* mElements;
    Element* mpOffsetedElements;
  };

}
#endif
