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

#ifndef __VCGLIB_UGRID_2D
#define __VCGLIB_UGRID_2D

#include <vector>
#include <algorithm>
#include <stdio.h>

#include <vcg/space/box2.h>
#include <vcg/space/line2.h>
#include <vcg/space/index/grid_util2d.h>
#include <vcg/space/index/grid_closest2d.h>

namespace vcg {

template < class OBJTYPE, class FLT=float >
class GridStaticPtr2D: public BasicGrid2D<FLT>, SpatialIndex2D<OBJTYPE,FLT>
{
public:
    typedef OBJTYPE ObjType;
    typedef ObjType* ObjPtr;
    typedef typename ObjType::ScalarType ScalarType;
    typedef Point2<ScalarType> CoordType;
    typedef Box2<ScalarType> Box2x;
    typedef GridStaticPtr2D<OBJTYPE,FLT> MyGridType;
    typedef typename std::vector<OBJTYPE*>::iterator CellIterator;
    std::vector<std::vector<std::vector<OBJTYPE*> > >  data;

private:



    void Add( const int x,
              const int y,
              ObjType *elem)
    {
        assert((x>=0)&&(x<(int)data.size()));
        assert((y>=0)&&(y<(int)data[x].size()));
        data[x][y].push_back(elem);
    }

    template <class OBJITER>
    inline void Set(const OBJITER & _oBegin,
                     const OBJITER & _oEnd,
                     const Box2x &_bbox,
                     Point2i _siz)
    {
        this->bbox=_bbox;
        this->siz=_siz;

        // find voxel size starting from the provided bbox and grid size.

        this->dim  = this->bbox.max - this->bbox.min;
        this->voxel.X() = (this->dim.X()/(ScalarType)this->siz.X());
        this->voxel.Y() = (this->dim.Y()/(ScalarType)this->siz.Y());

        ///allocate space
        data.resize(this->siz.X());
        for (size_t x=0;x<data.size();x++)
            data[x].resize(this->siz.Y());


        OBJITER IteObj;
        for (IteObj=_oBegin;IteObj!=_oEnd;IteObj++)
        {
            Box2<ScalarType> Box2D;
            (*IteObj).GetBBox(Box2D);

            //get index of intersected cells
            Point2i minIndex=this->GridP(Box2D.min);
            Point2i maxIndex=this->GridP(Box2D.max);

            for (int x=minIndex.X();x<=maxIndex.X();x++)
                for (int y=minIndex.Y();y<=maxIndex.Y();y++)
                    Add(x,y,&(*IteObj));
        }
    }

    template <class OBJITER>
    inline void Set(const OBJITER & _oBegin,
                    const OBJITER & _oEnd,
                    const Box2x &_bbox)
    {
      int _size=(int)std::distance<OBJITER>(_oBegin,_oEnd);

      Point2<FLT> _dim = _bbox.max - _bbox.min;

      Point2i _siz;

      BestDim2D( _size, _dim, _siz );

      Set(_oBegin,_oEnd,_bbox,_siz);
    }

public:

    void Grid( const int x, const int y,
               CellIterator & first,
               CellIterator & last )
    {
      first = data[x][y].begin();
      last  = data[x][y].end();
    }

    void Grid( const Point2i pi,
               CellIterator & first,
               CellIterator & last )
    {
        return Grid(pi.X(),pi.Y(),first,last);
    }

    template <class OBJITER>
    inline void Set(const OBJITER & _oBegin,
                    const OBJITER & _oEnd)
    {
      Box2x bbox,ibbox;

      OBJITER IteObj;
      for (IteObj=_oBegin;IteObj!=_oEnd;IteObj++)
      {
        (*IteObj).GetBBox(ibbox);
        bbox.Add(ibbox);
      }
      
      ScalarType diag=bbox.Diag();
      bbox.Offset(diag*0.01);
      Set<OBJITER>(_oBegin,_oEnd,bbox);
    }

    template <class OBJPOINTDISTFUNCTOR, class OBJMARKER>
    ObjPtr  GetClosest(OBJPOINTDISTFUNCTOR & _getPointDistance,
                       OBJMARKER & _marker,
                       const typename OBJPOINTDISTFUNCTOR::QueryType & _p,
                       const ScalarType & _maxDist,
                       ScalarType & _minDist,
                       CoordType & _closestPt)
    {
        return (vcg::GridClosest2D<MyGridType,OBJPOINTDISTFUNCTOR,OBJMARKER>(*this,_getPointDistance,_marker, _p,_maxDist,_minDist,_closestPt));
    }

    template <class OBJMARKER, class OBJPTRCONTAINER>
    unsigned int GetInBox(OBJMARKER & _marker,
                          const vcg::Box2<ScalarType> _bbox,
                          OBJPTRCONTAINER & _objectPtrs)
    {
        return(vcg::GridGetInBox2D<MyGridType,OBJMARKER,OBJPTRCONTAINER>
               (*this,_marker,_bbox,_objectPtrs));
    }


    template <class OBJRAYISECTFUNCTOR, class OBJMARKER>
    ObjPtr DoRay(OBJRAYISECTFUNCTOR & _rayIntersector, OBJMARKER & _marker,
                 const Ray2<ScalarType> & _ray, const ScalarType & _maxDist,
                 ScalarType & _t)
    {
        return(vcg::GridDoRay2D<MyGridType,OBJRAYISECTFUNCTOR,OBJMARKER>(*this,_rayIntersector,_marker,_ray,_maxDist,_t));
    }


}; //end class GridStaticPtr_2D

} // end namespace

#endif
