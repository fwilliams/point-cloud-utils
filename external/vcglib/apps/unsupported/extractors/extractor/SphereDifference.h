#ifndef __VCGTEST_SPHEREDIFFERENCE
#define __VCGTEST_SPHEREDIFFERENCE

class SphereDifference
{
public:
	SphereDifference() 
	{}

	SphereDifference( const SphereDifference &sphere_difference)
	{
		_union	= sphere_difference._union;
		_sphere = sphere_difference._sphere;
	}

	SphereDifference( const SphereUnion &sphere_union,  const ImplicitSphere &sphere)
	{
		_union = sphere_union;
		_sphere = sphere;
	}

	float V(int x, int y, int z) 
	{
		return vcg::math::Max<float>(_union.V(x, y, z), -_sphere.V(x, y, z));
	}

	inline bool DirectedDistance(const vcg::Point3i p1, const vcg::Point3i p2, vcg::Point3f &p, vcg::Point3f &n, float &d)
	{
		vcg::Point3f v1, n1;
		vcg::Point3f v2, n2;
		float				 d1, d2;

		bool ok1 = _union.DirectedDistance(p1, p2, v1, n1, d1);
		bool ok2 = _sphere.DirectedDistance(p1, p2, v2, n2, d2);
		d2 = -d2;

		if (ok1 && ok2) 
		{ 
			if (d1 > d2) 
				ok2 = false; 
			else 
				ok1 = false; 
		}
    
    if (ok1)  
		{ 
			p = v1; 
			n = n1; 
			d = d1; 
			return true; 
		}
    else if (ok2)  
		{ 
			p = v2; 
			n = n2; 
			d = d2; 
			return true; 
		}
    
    return false;
 	}

	private:
	SphereUnion			_union;
	ImplicitSphere	_sphere;

};

#endif // __VCGTEST_SPHEREDIFFERENCE