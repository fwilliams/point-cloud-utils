#ifndef __VCG_SIMPLE_VOLUME
#define __VCG_SIMPLE_VOLUME
#include<vector>
namespace vcg 
{

template <class VOX_TYPE>
class SimpleVolume
{
public:
  typedef VOX_TYPE VoxelType;
  std::vector<VoxelType> Vol;
  
  Point3i sz;   /// Dimensioni griglia come numero di celle per lato
  
  const Point3i &ISize() {return sz;};   /// Dimensioni griglia come numero di celle per lato
 	
  void Init(Point3i _sz)
  {
    sz=_sz;
    Vol.resize(sz[0]*sz[1]*sz[2]);
  }

  float Val(const int &x,const int &y,const int &z) const {
      return cV(x,y,z).V(); 
    //else return numeric_limits<float>::quiet_NaN( ); 
  }

  float &Val(const int &x,const int &y,const int &z) {
      return V(x,y,z).V(); 
    //else return numeric_limits<float>::quiet_NaN( ); 
  }

	VOX_TYPE &V(const int &x,const int &y,const int &z) {
		return Vol[x+y*sz[0]+z*sz[0]*sz[1]]; 
	}

	const VOX_TYPE &cV(const int &x,const int &y,const int &z) const {
		return Vol[x+y*sz[0]+z*sz[0]*sz[1]]; 
	}


typedef enum { XAxis=0,YAxis=1,ZAxis=2} VolumeAxis;

template < class VertexPointerType,  VolumeAxis AxisVal >
  void GetIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{
			float f1 = Val(p1.X(), p1.Y(), p1.Z())-thr;
			float f2 = Val(p2.X(), p2.Y(), p2.Z())-thr;
			float u = (float) f1/(f1-f2);
			if(AxisVal==XAxis) v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
                    else v->P().X() = (float) p1.X();
			if(AxisVal==YAxis) v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
			              else v->P().Y() = (float) p1.Y();
			if(AxisVal==ZAxis) v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
			              else v->P().Z() = (float) p1.Z();
}

template < class VertexPointerType >
  void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,XAxis>(p1,p2,v,thr); }

template < class VertexPointerType >
  void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,YAxis>(p1,p2,v,thr); }

template < class VertexPointerType >
  void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
{ GetIntercept<VertexPointerType,ZAxis>(p1,p2,v,thr); }
};
template <class VolumeType>
class RawVolumeImporter
{
public:
  enum DataType
{
		// Funzioni superiori
  UNDEF=0,
  BYTE=1,
  SHORT=2,
  FLOAT=3
};

static bool Open(const char *filename, VolumeType &V, Point3i sz, DataType d)
{
return true;
}
};

class SimpleVoxel
{
private:
  float _v;
public:
  float &V() {return _v;};
  float V() const {return _v;};
};
} // end namespace 
#endif // __VCG_SIMPLE_VOLUME