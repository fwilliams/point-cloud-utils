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
/****************************************************************************
History

$Log: not supported by cvs2svn $
Revision 1.6  2007/05/08 12:11:58  pietroni
added circle-line intersection


****************************************************************************/



#ifndef __VCGLIB_INTERSECTION_2
#define __VCGLIB_INTERSECTION_2
#include <vcg/space/line2.h>
#include <vcg/space/ray2.h>
#include <vcg/space/segment2.h>
#include <vcg/space/point2.h>
#include <vcg/space/triangle2.h>
#include <vcg/space/box2.h>
#include <vector>



namespace vcg {
	/** \addtogroup space */
	/*@{*/
	/** 
	Function computing the intersection between couple of geometric primitives in
	2 dimension
	*/

	/// return true if the algle is convex (right rotation)
	template<class SCALAR_TYPE>
	inline bool Convex(const Point2<SCALAR_TYPE> & p0,const Point2<SCALAR_TYPE> & p1,const Point2<SCALAR_TYPE> & p2)
	{
		const SCALAR_TYPE EPS= SCALAR_TYPE(1e-8);
		return (((p0-p1)^(p2-p1))<=EPS);
	}

	///return if exist the intersection point
	///between 2 lines in a 2d plane
	template<class SCALAR_TYPE>
	inline bool LineLineIntersection(const vcg::Line2<SCALAR_TYPE> & l0,
		const vcg::Line2<SCALAR_TYPE> & l1,
		Point2<SCALAR_TYPE> &p)
	{
		const SCALAR_TYPE Eps= SCALAR_TYPE(1e-8);
		///first line
		SCALAR_TYPE x1=l0.Origin().X();
		SCALAR_TYPE y1=l0.Origin().Y(); 
		SCALAR_TYPE x2=x1+l0.Direction().X();
		SCALAR_TYPE y2=y1+l0.Direction().Y(); 

		///second line
		SCALAR_TYPE x3=l1.Origin().X();
		SCALAR_TYPE y3=l1.Origin().Y(); 
		SCALAR_TYPE x4=x3+l1.Direction().X();
		SCALAR_TYPE y4=y3+l1.Direction().Y(); 

		///then  find intersection

		///denominator
		SCALAR_TYPE den=((x1-x2)*(y3-y4))-((y1-y2)*(x3-x4));
		if (fabs(den)<Eps)
			return false;

		SCALAR_TYPE d0=(x1*y2)-(y1*x2);
		SCALAR_TYPE d1=(x3*y4)-(y3*x4);
		SCALAR_TYPE numx=(d0*(x3-x4))-(d1*(x1-x2));
		SCALAR_TYPE numy=(d0*(y3-y4))-(d1*(y1-y2));

		p.X()=numx/den;
		p.Y()=numy/den;
		return true;
	}

	///return if exist the intersection point
	///between 2 lines in a 2d plane
	template<class SCALAR_TYPE>
	inline bool RayLineIntersection(const vcg::Line2<SCALAR_TYPE> & l,
		const vcg::Ray2<SCALAR_TYPE> & r,
		Point2<SCALAR_TYPE> &p)
	{
		///construct line from ray
		vcg::Line2<SCALAR_TYPE> l_test;
		l_test.Set(r.Origin(),r.Direction());
		if (!LineLineIntersection(l,l_test,p))
			return false;
		Point2<SCALAR_TYPE> dir=p-r.Origin();
		dir.Normalize();
		return (dir*r.Direction()>0);
	}


	/// interseciton between point and triangle
	template<class SCALAR_TYPE>
	inline bool RaySegmentIntersection(const vcg::Ray2<SCALAR_TYPE> & r,
		const vcg::Segment2<SCALAR_TYPE> &seg,
		Point2<SCALAR_TYPE> &p_inters)
	{
		///first compute intersection between lines
		vcg::Line2<SCALAR_TYPE> line2;
		line2.SetOrigin(seg.P0());
		vcg::Point2<SCALAR_TYPE> dir=seg.P1()-seg.P0();
		dir.Normalize();
		line2.SetDirection(dir);
		if(!RayLineIntersection<SCALAR_TYPE>(line2,r,p_inters))
			return false;
		///then test if intersection point is nearest 
		///to both extremes then length of the segment 
		SCALAR_TYPE d0=(seg.P1()-p_inters).Norm();
		SCALAR_TYPE d1=(seg.P0()-p_inters).Norm();
		SCALAR_TYPE length=(seg.P0()-seg.P1()).Norm();
		return ((d0<length)&&(d1<length));
	}

    /// interseciton between point and triangle
    template<class SCALAR_TYPE>
    inline bool RayBoxIntersection(const vcg::Ray2<SCALAR_TYPE> & r,
                                   const vcg::Box2<SCALAR_TYPE> &bbox,
                                    Point2<SCALAR_TYPE> &p_inters)
    {
        ///first create the 4 segments
        vcg::Segment2<SCALAR_TYPE> S[4];
        for (int i=0;i<4;i++)
            S[i]=vcg::Segment2<SCALAR_TYPE>(bbox.P(i),bbox.P((i+1)%4));

        SCALAR_TYPE mind=std::numeric_limits<SCALAR_TYPE>::max();
        bool found=false;
        for (int i=0;i<4;i++)
        {
             Point2<SCALAR_TYPE> p_inters_test;
            if (!RaySegmentIntersection(r,S[i],p_inters_test))continue;
            SCALAR_TYPE Norm=(p_inters_test-r.Origin()).Norm();
            if (Norm<mind)
            {
                mind=Norm;
                p_inters=p_inters_test;
                found=true;
            }
        }
        return found;
    }

	/// interseciton between point and triangle
	template<class SCALAR_TYPE>
	inline bool LineSegmentIntersection(const vcg::Line2<SCALAR_TYPE> & line,
		const vcg::Segment2<SCALAR_TYPE> &seg,
		Point2<SCALAR_TYPE> &p_inters)
	{
		///first compute intersection between lines
		vcg::Line2<SCALAR_TYPE> line2;
		line2.SetOrigin(seg.P0());
		vcg::Point2<SCALAR_TYPE> dir=seg.P1()-seg.P0();
		dir.Normalize();
		line2.SetDirection(dir);
		if(!LineLineIntersection(line,line2,p_inters))
			return false;
		///then test if intersection point is nearest 
		///to both extremes then length of the segment 
		SCALAR_TYPE d0=(seg.P1()-p_inters).Norm();
		SCALAR_TYPE d1=(seg.P0()-p_inters).Norm();
		SCALAR_TYPE length=(seg.P0()-seg.P1()).Norm();
		return ((d0<length)&&(d1<length));
	}

	/// interseciton between two segments
	template<class SCALAR_TYPE>
	inline bool SegmentSegmentIntersection(const vcg::Segment2<SCALAR_TYPE> &seg0,
		const vcg::Segment2<SCALAR_TYPE> &seg1,
		Point2<SCALAR_TYPE> &p_inters)
	{
		const SCALAR_TYPE Eps= SCALAR_TYPE(1e-8);
		SCALAR_TYPE lambda0,lambda1;
        const Point2<SCALAR_TYPE> & p0 = seg0.P0();
        const Point2<SCALAR_TYPE> & p1 = seg0.P1();
        const Point2<SCALAR_TYPE> & p2 = seg1.P0();
        const Point2<SCALAR_TYPE> & p3 = seg1.P1();

		SCALAR_TYPE a = (p1-p0)[0];
		SCALAR_TYPE b = (p2-p3)[0];
		SCALAR_TYPE c = (p1-p0)[1];
		SCALAR_TYPE d = (p2-p3)[1];

		SCALAR_TYPE e = (p2-p0)[0];
        SCALAR_TYPE f = (p2-p0)[1];

		SCALAR_TYPE det = a*d-b*c;

		lambda0 = (d*e-b*f)/det;
		lambda1 = (-c*e+a*f)/det;
        if (fabs(det)<Eps)
			return false;// they are parallell
		
		if (!(lambda0 >= 0.0 && lambda0 <= 1.0 && lambda1 >= 0.0 && lambda1 <= 1.0))
			return false;
        p_inters = p0*(1-lambda0)+p1*lambda0;
		return true;
	}
	/// interseciton between point and triangle
	template<class SCALAR_TYPE>
	inline bool IsInsideTrianglePoint( const Triangle2<SCALAR_TYPE> & t,const Point2<SCALAR_TYPE> & p)
	{
		Point2<SCALAR_TYPE> p0=t.P0(0);
		Point2<SCALAR_TYPE> p1=t.P0(1);
		Point2<SCALAR_TYPE> p2=t.P0(2);

		///first test with bounding box
		vcg::Box2<SCALAR_TYPE> b2d;
		b2d.Add(p0);
		b2d.Add(p1);
		b2d.Add(p2);
		if (!b2d.IsIn(p))
			return false;

		///then text convex
		if (!Convex(p0,p1,p2))
			std::swap<Point2<SCALAR_TYPE> >(p1,p2);
		return((Convex(p,p0,p1))&&(Convex(p,p1,p2))&&(Convex(p,p2,p0)));
		//return((Convex(p,p0,p1))&&(Convex(p,p1,p2))&&(Convex(p,p2,p0)));
	}

	template<class ScalarType>
	bool TriangleTriangleIntersect2D(const vcg::Triangle2<ScalarType> &tr0,
		const vcg::Triangle2<ScalarType> &tr1)
	{
		///test BBox Intersection
		vcg::Box2<ScalarType> bbtr0;
		bbtr0.Add(tr0.P(0));
		bbtr0.Add(tr0.P(1));
		bbtr0.Add(tr0.P(2));
		vcg::Box2<ScalarType> bbtr1;
		bbtr1.Add(tr1.P(0));
		bbtr1.Add(tr1.P(1));
		bbtr1.Add(tr1.P(2));
		if (!bbtr0.Collide(bbtr1)) return false;
		///test vertex in face
		for (int i=0;i<3;i++)
		{
			bool inside0=vcg::IsInsideTrianglePoint(tr0,tr1.P(i));
			bool inside1=vcg::IsInsideTrianglePoint(tr1,tr0.P(i));
			if (inside0 || inside1) return true;
		}
		///test segment 
		///to segment intersection
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++)
			{
				if (i>j) continue;
				vcg::Segment2<ScalarType> seg0=vcg::Segment2<ScalarType>(tr0.P(i),tr0.P((i+1)%3));
				vcg::Segment2<ScalarType> seg1=vcg::Segment2<ScalarType>(tr1.P(j),tr1.P((j+1)%3));
				vcg::Point2<ScalarType> p_inters;
				bool intersect=SegmentSegmentIntersection(seg0,seg1,p_inters);
				if (intersect) return true;
			}
		}
		return false;
	}
		
	template <class ScalarType>
	bool PointInsidePolygon(vcg::Point2<ScalarType> p,
							const std::vector<vcg::Segment2<ScalarType> > &polygon)
	{
		int n=polygon.size();
		vcg::Box2<ScalarType> BB;
		for (int i=0;i<n;i++)
		{
			BB.Add(polygon[i].P0());
			BB.Add(polygon[i].P1());
		}
		if (!BB.IsIn(p))return false;
		ScalarType size=BB.Diag();
		///take 4 directions
		int inside_test=0;
		for (int dir=0;dir<4;dir++)
		{
			int intersection=0;
			vcg::Ray2<ScalarType> r;
			vcg::Point2<ScalarType> direct=vcg::Point2<ScalarType>(0,0);
			switch (dir) 
			{
				case 0 : direct.X()=1;break;
				case 1 : direct.Y()=1;break;
				case 2 : direct.X()=-1; break;
				default :direct.Y()=-1;
			}			
			r.SetOrigin(p);
			r.SetDirection(direct);
			for (int i=0;i<n;i++) 
			{
				Point2<ScalarType> p_inters;
				if (vcg::RaySegmentIntersection(r,polygon[i],p_inters))intersection++;
			}
			if ((intersection%2)==1)
				inside_test++;
		}
		return(inside_test>2);
	}

	//intersection between a circle and a line
	template<class ScalarType>
	inline bool CircleLineIntersection(const vcg::Line2<ScalarType> & line,
		const vcg::Point2<ScalarType> &center,
		const ScalarType &radius,
		vcg::Point2<ScalarType> &p0,
		vcg::Point2<ScalarType> &p1)
	{
		///translate with origin on the center
		ScalarType x1,x2,y1,y2;
		x1=line.Origin().X()-center.X();
		y1=line.Origin().Y()-center.Y();
		x2=x1+line.Direction().X();
		y2=y1+line.Direction().Y();

		ScalarType dx,dy,dr,D,delta,sign;
		dx=x2-x1;
		dy=y2-y1;
		dr=sqrt(dx*dx+dy*dy);
		D=x1*y2-x2*y1;
		delta=radius*radius*dr*dr-D*D;
		if (dy>=0)
			sign=1;
		else
			sign=-1;

		if (delta<0.000001)
			return false;///no intersection
		else
		{
			p0.X()=(D*dy+sign*dx*sqrt(delta))/dr*dr;
			p0.Y()=(-D*dx+fabs(dy)*sqrt(delta))/dr*dr;
			p1.X()=(D*dy-sign*dx*sqrt(delta))/dr*dr;
			p1.Y()=(-D*dx-fabs(dy)*sqrt(delta))/dr*dr;
			p0+=center;
			p1+=center;
			return true;
		}
	}


    // Ray-Segment Functor
    class RaySegmentIntersectionFunctor {
    public:

        template <class SEGMENTTYPE, class SCALARTYPE>
        inline bool operator () (const SEGMENTTYPE & S,
                                 const Ray2<SCALARTYPE> & ray,
                                 SCALARTYPE & t)
        {
            typedef SCALARTYPE ScalarType;
            typedef vcg::Point2<ScalarType> CoordType;

            CoordType inters_test;
            bool bret = RaySegmentIntersection(ray,S, inters_test);
            if (bret)
                t=(inters_test-ray.Origin()).Norm();
            return (bret);
        }
    };

	/*@}*/
} // end namespace
#endif
