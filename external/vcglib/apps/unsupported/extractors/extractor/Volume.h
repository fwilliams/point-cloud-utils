#ifndef __VCGTEST_VOLUME
#define __VCGTEST_VOLUME

#include "ImplicitSphere.h"
#include "SphereUnion.h"
#include "SphereDifference.h"

class Volume
{
public:
	Volume()
	{
		ImplicitSphere     s1(vcg::Point3f(-5.0,  0.0,  0.0), 10.0);
		ImplicitSphere     s2(vcg::Point3f( 5.0,  5.0,  3.0), 7.0);
		ImplicitSphere     s3(vcg::Point3f( 1.0,  0.0,  10.0), 6.0);
		SphereUnion        sphere_union(s1, s2);
		SphereDifference	 sphere_difference(sphere_union, s3);
		_sphere_diff = sphere_difference;		
	}

	float V(const int pi, const int pj, const int pk)
	{
		return _sphere_diff.V(pi, pj, pk);
	}

	void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
	{
		vcg::Point3f p, n;
		float d;
		if (_sphere_diff.DirectedDistance(p1, p2, p, n, d))
		{
			v->P() = p;
			v->N() = n;
		}
		else
		{
			float f1 = V(p1.X(), p1.Y(), p1.Z());
			float f2 = V(p2.X(), p2.Y(), p2.Z());
			float u = (float) f1/(f1-f2);
			v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
			v->P().Y() = (float) p1.Y();
			v->P().Z() = (float) p1.Z();
			
		}
	}
	void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
	{
		vcg::Point3f p, n;
		float d;
		if (_sphere_diff.DirectedDistance(p1, p2, p, n, d))
		{
			v->P() = p;
			v->N() = n;
		}
		else
		{
			float f1 = V(p1.X(), p1.Y(), p1.Z());
			float f2 = V(p2.X(), p2.Y(), p2.Z());
			float u = (float) f1/(f1-f2);
			v->P().X() = (float) p1.X();
			v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
			v->P().Z() = (float) p1.Z();
			
		}
	}
	void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
	{
		vcg::Point3f p, n;
		float d;
		if (_sphere_diff.DirectedDistance(p1, p2, p, n, d))
		{
			v->P() = p;
			v->N() = n;
		}
		else
		{
			float f1 = V(p1.X(), p1.Y(), p1.Z());
			float f2 = V(p2.X(), p2.Y(), p2.Z());
			float u = (float) f1/(f1-f2);
			v->P().X() = (float) p1.X();
			v->P().Y() = (float) p1.Y();
			v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
			
		}
	}

private:
  SphereDifference _sphere_diff;
};

#endif // __VCGTEST_VOLUME