#ifndef __VCGTEST_IMPLICITSPHERE
#define __VCGTEST_IMPLICITSPHERE

class ImplicitSphere
{
public:
	ImplicitSphere()
	{
		_center.SetZero();
		_radius = _sqr_radius = 0.0;
	};

	ImplicitSphere(vcg::Point3f &center, float radius)
	{
		_center = center;
		_radius = radius;
		_sqr_radius = radius*radius;
	};

	ImplicitSphere(const ImplicitSphere &sphere)
	{
		_center = sphere._center;
		_radius = sphere._radius;
		_sqr_radius = sphere._sqr_radius;
	};

	ImplicitSphere& operator=(const ImplicitSphere &sphere)
	{
		if (this != &sphere)
		{
			_center = sphere._center;
			_radius = sphere._radius;
			_sqr_radius = sphere._sqr_radius;
		}
		return *this;
	};
		
	~ImplicitSphere() 
	{};

	bool operator!=(const ImplicitSphere &sphere)
	{
		return (sphere._center!=_center && sphere._radius!=_radius);
	};


	float V(int x, int y, int z) const
	{
		vcg::Point3f point((float) x, (float) y, (float) z);
		return (_center-point).Norm() - _radius;
	};

	bool DirectedDistance(const vcg::Point3i &p1, const vcg::Point3i &p2, vcg::Point3f &v, vcg::Point3f &n, float &dist)
	{
		vcg::Point3f orig, dir;
		orig.X() = (float) p1.X();				
		orig.Y() = (float) p1.Y();				
		orig.Z() = (float) p1.Z();
		dir.X()	 = (float) p2.X()-p1.X();	
		dir.Y()  = (float) p2.Y()-p1.Y();	
		dir.Z()	 = (float) p2.Z()-p1.Z();

		double a = dir.SquaredNorm();
		double b = 2.0*(dir*(orig - _center));
		double c = (orig - _center).SquaredNorm() - _radius*_radius;
		double d = b*b - 4.0*a*c;
		
		if (d >= 0)
		{
			d = sqrt(d);

			double t1 = (-b-d) / (2.0*a);
			double t2 = (-b+d) / (2.0*a);
			double t  = 1.00001;
			if (t1 >= 0.0 && t1 < t) t = t1;
			if (t2 >= 0.0 && t2 < t) t = t2;

			if (t != 1.00001)
			{
				v = (vcg::Point3f) (orig + dir*((float)t));
				n = (vcg::Point3f) ((v - _center) / _radius);
				dist = (float) ((dir*n) < 0.0) ? dir.Norm()*t : -dir.Norm()*t;
				return true;
			}
		}
		return false;
	};

private:
	vcg::Point3f	_center;
	float					_radius;
	float					_sqr_radius;
};

#endif // __VCGTEST_IMPLICITSPHERE