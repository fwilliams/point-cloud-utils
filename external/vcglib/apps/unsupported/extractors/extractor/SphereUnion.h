#ifndef __VCGTEST_SPHEREUNION
#define __VCGTEST_SPHEREUNION

class SphereUnion
{
public:
	SphereUnion()
	{};

	SphereUnion(const ImplicitSphere &sphere1, const ImplicitSphere &sphere2)
	{
		_sphere1 = sphere1; 
		_sphere2 = sphere2;
	};

	SphereUnion(const SphereUnion &sphere_union)
	{
		_sphere1 = sphere_union._sphere1;
		_sphere2 = sphere_union._sphere2;
	}

	SphereUnion& operator=(const SphereUnion &sphere_union)
	{
		if (this != &sphere_union)
		{
			_sphere1 = sphere_union._sphere1;
			_sphere2 = sphere_union._sphere2;
		}
		return *this;
	}

	bool operator!=(const SphereUnion &sphere_union)
	{
		bool comp1 = _sphere1 != sphere_union._sphere1;
		bool comp2 = _sphere2 != sphere_union._sphere2;
		return (comp1 && comp2);
	}

	float V(int x, int y, int z)
	{
		return vcg::math::Min<float>(_sphere1.V(x, y, z), _sphere2.V(x, y, z));
	};

	bool DirectedDistance(const vcg::Point3i &p1, const vcg::Point3i &p2, vcg::Point3f &v, vcg::Point3f &n, float &d)
	{
		vcg::Point3f v1, n1;
		vcg::Point3f v2, n2;
		float				 d1, d2;

		bool ok1 = _sphere1.DirectedDistance(p1, p2, v1, n1, d1);
		bool ok2 = _sphere2.DirectedDistance(p1, p2, v2, n2, d2);

		if (ok1 && ok2)
		{ 
			if (d1 < d2) 
				ok2 = false; 
			else 
				ok1 = false; 
		}
		
		if (ok1)
		{
			v = v1;
			n = n1;
			d = d1;
			return true;
		}
		else if (ok2)
		{
			v = v2;
			n = n2;
			d = d2;
			return true;
		}
		else
			return false;
	};

private:
	ImplicitSphere _sphere1;
	ImplicitSphere _sphere2;
};

#endif // __VCGTEST_SPHEREUNION