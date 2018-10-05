#ifndef __EXTRS_IMPLICIT
#define __EXTRS_IMPLICIT

#include <vcg/space/point3.h>

class Implicit
{
public:
	Implicit() {};

	virtual ~Implicit() {};

	virtual float V(int x, int y, int z) const = 0;

	virtual vcg::Point3f N(float x, float y, float z) const = 0;
};

#endif // __EXTRS_IMPLICIT