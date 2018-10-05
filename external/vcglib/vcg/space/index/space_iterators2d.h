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

#ifndef __VCGLIB_SPATIAL_ITERATORS_2D
#define __VCGLIB_SPATIAL_ITERATORS_2D

#include <vector>
#include <vcg/space/intersection2.h>
#include <vcg/space/point2.h>
#include <vcg/space/box2.h>
#include <vcg/space/ray2.h>
#include <vcg/math/base.h>
#include <algorithm>
#include <float.h>
#include <limits>



namespace vcg{
template <class Spatial_Idexing,class INTFUNCTOR,class TMARKER>
class RayIterator2D
{
public:
    typedef typename Spatial_Idexing::ScalarType ScalarType;
    typedef typename vcg::Ray2<ScalarType> RayType;
    typedef typename Spatial_Idexing::Box2x IndexingBoxType;
protected:
    typedef typename Spatial_Idexing::ObjType ObjType;
    typedef typename vcg::Point2<ScalarType>  CoordType;
    typedef typename Spatial_Idexing::CellIterator CellIterator;
    ScalarType max_dist;


    bool _controlEnd()
    {
        return (currBox.Collide(Si.bbox));
    }


    void _NextCell()
    {
        currBox.min+=step;
        currBox.max+=step;
        dist+=step.Norm();
        end=!_controlEnd();
    }

    //refresh current cell intersection ,
    // return true if there is at lest 1 intersection
    bool Refresh()
    {
        std::vector<ObjType*> objectPtrs;
        GridGetInBox2D(Si,tm,currBox,objectPtrs,false);
        //printf(" size %d \n",objectPtrs.size());
        for(size_t i=0;i<objectPtrs.size();i++)
        {
            ObjType* elem=objectPtrs[i];
            if (elem->IsD())continue;
            if (tm.IsMarked(elem))continue;
            //tm.Mark(elem);

            ScalarType t;
            CoordType Int;
            if((int_funct((*elem),r,t))&&
               (t<=max_dist))
            {
                Int=r.Origin()+r.Direction()*t;
                Elems.push_back(Entry_Type(elem,t,Int));
            }
        }
        if (Elems.size()==0) return false;
        //then control if there are more than 1 element
        std::sort(Elems.begin(),Elems.end());
        CurrentElem=Elems.rbegin();

        return(Dist()<dist);
    }

public:


    //contructor
    RayIterator2D(Spatial_Idexing &_Si,
                  INTFUNCTOR &_int_funct,
                  const ScalarType &_max_dist,
                  TMARKER  &_tm)
        :Si(_Si),int_funct(_int_funct),tm(_tm)
    {
        max_dist=_max_dist;
    };



    void Init(const RayType _r)
    {
        r=_r;
        r.Normalize();

        //initialization
        end=false;
        tm.UnMarkAll();
        Elems.clear();

        CoordType start;

        //control if intersect the bounding box of the grid
        if (Si.bbox.IsIn(r.Origin()))
            start=r.Origin();
        else
            if (!(vcg::RayBoxIntersection<ScalarType>(r,Si.bbox,start)))
            {
                end=true;
                return;
            }


        stepsize=Si.voxel.Norm()*2;
        step=r.Direction()*stepsize;

        //create initial BB, inflate in case the direction is orthogonal to one axis
        currBox.SetNull();
        currBox.Add(start);
        currBox.Add(start+step);
        ScalarType diag=currBox.Diag();
        currBox.Offset(diag*0.01);
        dist=currBox.Diag();
        end=!_controlEnd();

        while ((!End()) && (!Refresh()))
            _NextCell();

        fflush(stdout);
    }

    bool End()
    {return end;}


    ObjType &operator *(){return *((*CurrentElem).elem);}

    CoordType IntPoint()
    {return ((*CurrentElem).intersection);}

    ScalarType Dist()
    {
        if (Elems.size()>0)
            return ((*CurrentElem).dist);
        else
            return ((ScalarType)FLT_MAX);
    }

    void operator ++()
    {
         if (!Elems.empty()) Elems.pop_back();

         CurrentElem = Elems.rbegin();

         if (Dist()>dist)
         {
             if (!End())
             {
                 _NextCell();
                 while ((!End()) && (!Refresh()))
                     _NextCell();
             }
         }
     }

protected:

    ///structure that mantain for the current cell pre-calculated data
    struct Entry_Type
    {
    public:

        Entry_Type(ObjType* _elem,ScalarType _dist,CoordType _intersection)
        {
            elem=_elem;
            dist=_dist;
            intersection=_intersection;
        }

        Entry_Type(const Entry_Type &e)
        {
            elem=e.elem;
            dist=e.dist;
            intersection=e.intersection;
        }

        inline bool operator <  ( const Entry_Type & l ) const{return (dist > l.dist); }
        ObjType* elem;
        ScalarType dist;
        CoordType intersection;
    };

    RayType r;							//ray to find intersections
    Spatial_Idexing &Si;	  //reference to spatial index algorithm
    bool end;								//true if the scan is terminated
    INTFUNCTOR &int_funct;
    TMARKER &tm;

    std::vector<Entry_Type> Elems;					//element loaded from curren cell
    typedef typename std::vector<Entry_Type>::reverse_iterator ElemIterator;
    ElemIterator CurrentElem;	//iterator to current element

    vcg::Box2<ScalarType> currBox;
    CoordType step;
    ScalarType stepsize;
    ScalarType dist;

};

}

#endif
