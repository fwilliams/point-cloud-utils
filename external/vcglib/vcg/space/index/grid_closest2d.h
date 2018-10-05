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

#ifndef __VCGLIB_GRID_CLOSEST_2D
#define __VCGLIB_GRID_CLOSEST_2D

#include <vcg/space/index/space_iterators2d.h>
#include <vcg/space/intersection2.h>

namespace vcg{


template <class SPATIALINDEXING,class OBJMARKER, class OBJPTRCONTAINER>
unsigned int GridGetInBox2D(SPATIALINDEXING &_Si,
                            OBJMARKER & _marker,
                            const vcg::Box2<typename SPATIALINDEXING::ScalarType> &_bbox,
                            OBJPTRCONTAINER & _objectPtrs,
                            bool update_global_mark=true)
{
    typename SPATIALINDEXING::CellIterator first,last,l;
    //_objectPtrs.clear();
    vcg::Box2i ibbox;
    Box2i Si_ibox(Point2i(0,0),_Si.siz-Point2i(1,1));
    _Si.BoxToIBox(_bbox, ibbox);
    ibbox.Intersect(Si_ibox);

    if (update_global_mark)
        _marker.UnMarkAll();

    if (ibbox.IsNull())
        return 0;
    else
    {
        int ix,iy;
        for (ix=ibbox.min[0]; ix<=ibbox.max[0]; ix++)
            for (iy=ibbox.min[1]; iy<=ibbox.max[1]; iy++)
            {
                _Si.Grid( ix, iy, first, last );
                for(l=first;l!=last;++l)
                    if (!(**l).IsD())
                    {
                        typename SPATIALINDEXING::ObjPtr elem=&(**l);
                        vcg::Box2<typename SPATIALINDEXING::ScalarType> box_elem;
												elem->GetBBox(box_elem);
                        if (update_global_mark)
                        {
                            if(( ! _marker.IsMarked(elem))&&(box_elem.Collide(_bbox)))
                            {
                                _objectPtrs.push_back(elem);
                                _marker.Mark(elem);
                            }
                        }else
                        {
                            if(box_elem.Collide(_bbox))
                             _objectPtrs.push_back(elem);
                        }
                    }
            }
        return (static_cast<unsigned int>(_objectPtrs.size()));
    }
}

template <class SPATIAL_INDEX, class OBJPOINTDISTFUNCTOR, class OBJMARKER>
typename SPATIAL_INDEX::ObjPtr  GridClosest2D(SPATIAL_INDEX &Si,
                                              OBJPOINTDISTFUNCTOR _getPointDistance,
                                              OBJMARKER & _marker,
                                              const typename OBJPOINTDISTFUNCTOR::QueryType  & _p_obj,
                                              const typename SPATIAL_INDEX::ScalarType & _maxDist,
                                              typename SPATIAL_INDEX::ScalarType & _minDist,
                                              typename SPATIAL_INDEX:: CoordType &_closestPt)
{
    typedef typename SPATIAL_INDEX::ObjPtr ObjPtr;
    typedef SPATIAL_INDEX SpatialIndex;
    typedef typename SPATIAL_INDEX::CoordType CoordType;
    typedef typename SPATIAL_INDEX::ScalarType ScalarType;
    typedef typename SPATIAL_INDEX::Box2x Box2x;

    Point2<ScalarType> _p = OBJPOINTDISTFUNCTOR::Pos(_p_obj);

    // Initialize min_dist with _maxDist to exploit early rejection test.
    _minDist = _maxDist;

    ObjPtr winner=NULL;
    _marker.UnMarkAll();
    ScalarType newradius = Si.voxel.Norm();
    ScalarType radius;
    Box2i iboxdone,iboxtodo;
    CoordType t_res;
    typename SPATIAL_INDEX::CellIterator first,last,l;
    if(Si.bbox.IsInEx(_p))
    {
        Point2i _ip;
        Si.PToIP(_p,_ip);
        Si.Grid( _ip[0],_ip[1], first, last );
        for(l=first;l!=last;++l)
        {
            ObjPtr elem=&(**l);
            if (!elem->IsD())
            {
                if (_getPointDistance((**l), _p_obj,_minDist, t_res))
                {
                    winner=elem;
                    _closestPt=t_res;
                    newradius=_minDist; //
                }
                _marker.Mark(elem);
            }
        }
        iboxdone=Box2i(_ip,_ip);
    }

    int ix,iy;
    Box2i ibox(Point2i(0,0),Si.siz-Point2i(1,1));

    do
    {
        radius=newradius;
        Box2x boxtodo=Box2x(_p,radius);
        //boxtodo.Intersect(Si.bbox);
        Si.BoxToIBox(boxtodo, iboxtodo);
        iboxtodo.Intersect(ibox);
        if(!boxtodo.IsNull())
        {
            for (ix=iboxtodo.min[0]; ix<=iboxtodo.max[0]; ix++)
                for (iy=iboxtodo.min[1]; iy<=iboxtodo.max[1]; iy++)
                    if(ix<iboxdone.min[0] || ix>iboxdone.max[0] ||  // this test is to avoid to re-process already analyzed cells.
                            iy<iboxdone.min[1] || iy>iboxdone.max[1] )
                    {
                        Si.Grid( ix, iy, first, last );
                        for(l=first;l!=last;++l) if (!(**l).IsD())
                        {
                            ObjPtr elem=&(**l);
                            if (!elem->IsD())
                            {
                                if( ! _marker.IsMarked(elem))
                                {
                                    if (_getPointDistance((**l), _p_obj, _minDist, t_res))
                                    {
                                        winner=elem;
                                        _closestPt=t_res;
                                    };
                                    _marker.Mark(elem);
                                }
                            }
                        }
                    }
        }
        if(!winner) newradius=radius+Si.voxel.Norm();
        else newradius = _minDist;
        iboxdone=iboxtodo;
    }
    while (_minDist>radius);

    return winner;
}

template <class SPATIALINDEXING,class OBJRAYISECTFUNCTOR, class OBJMARKER>
typename SPATIALINDEXING::ObjPtr GridDoRay2D(SPATIALINDEXING &_Si,
                                           OBJRAYISECTFUNCTOR &_rayIntersector,
                                           OBJMARKER &_marker,
                                           const Ray2<typename SPATIALINDEXING::ScalarType> & _ray,
                                           const typename SPATIALINDEXING::ScalarType & _maxDist,
                                           typename SPATIALINDEXING::ScalarType & _t)
{
    typedef vcg::RayIterator2D<SPATIALINDEXING,OBJRAYISECTFUNCTOR,OBJMARKER> RayIteratorType;
    RayIteratorType RayIte=RayIteratorType(_Si,_rayIntersector,_maxDist,_marker);
    RayIte.Init(_ray);

    if (!RayIte.End())
    {
        _t=RayIte.Dist();
        return(&(*RayIte));
    }
    return NULL;
}


}//end namespace vcg
#endif

