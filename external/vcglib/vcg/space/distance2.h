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
/****************************************************************************/

#ifndef __VCG_DISTANCE2
#define __VCG_DISTANCE2

#include <limits>
#include <vcg/space/segment2.h>
#include <vcg/space/intersection2.h>

namespace vcg {

	/*
	* Computes the minimum distance between a two segments
	* @param[in] S0				The input segment0
	* @param[in] S1				The input segment1
	* return the distance between the two segments
	*/
	template<class ScalarType>
	ScalarType Segment2DSegment2DDistance(const vcg::Segment2<ScalarType> &S0,
									const vcg::Segment2<ScalarType> &S1,
                                    vcg::Point2<ScalarType> &p_clos0,
                                    vcg::Point2<ScalarType> &p_clos1)
	{
        //first test if they intersect
        vcg::Point2<ScalarType> IntPoint;
        if (vcg:: SegmentSegmentIntersection(S0,S1,IntPoint))
        {
            p_clos0=IntPoint;
            p_clos1=IntPoint;
            return 0;
        }
		vcg::Point2<ScalarType> Pclos0=ClosestPoint(S0,S1.P0());
		vcg::Point2<ScalarType> Pclos1=ClosestPoint(S0,S1.P1());
		vcg::Point2<ScalarType> Pclos2=ClosestPoint(S1,S0.P0());
		vcg::Point2<ScalarType> Pclos3=ClosestPoint(S1,S0.P1());
		ScalarType d0=(Pclos0-S1.P0()).Norm();
		ScalarType d1=(Pclos1-S1.P1()).Norm();
		ScalarType d2=(Pclos2-S0.P0()).Norm();
		ScalarType d3=(Pclos3-S0.P1()).Norm();

        //then return the minimuim distance
        if ((d0<d1)&&(d0<d2)&&(d0<d3))
        {
            p_clos0=Pclos0;
            p_clos1=S1.P0();
            return d0;
        }
        if ((d1<d0)&&(d1<d2)&&(d1<d3))
        {
            p_clos0=Pclos1;
            p_clos1=S1.P1();
            return d1;
        }
        if ((d2<d0)&&(d2<d1)&&(d2<d3))
        {
            p_clos0=S0.P0();
            p_clos1=Pclos2;
            return d2;
        }
        else
        {
            p_clos0=S0.P1();
            p_clos1=Pclos3;
            return d3;
        }
	}


}///end namespace vcg

#endif
