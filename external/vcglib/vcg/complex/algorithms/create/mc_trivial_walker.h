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
#ifndef __VCG_TRIVIAL_WALKER
#define __VCG_TRIVIAL_WALKER

#include<vcg/space/index/grid_util.h>

namespace vcg {

// Very simple volume class.
// just an example of the interface that the trivial walker expects

template <class VOX_TYPE>
class SimpleVolume : public BasicGrid<typename VOX_TYPE::ScalarType>
{
public:
  typedef VOX_TYPE VoxelType;
  typedef typename VoxelType::ScalarType ScalarType;
  typedef typename BasicGrid<typename VOX_TYPE::ScalarType>::Box3x Box3x;

  const Point3i &ISize() {return this->siz;}   /// Dimensioni griglia come numero di celle per lato

  ScalarType Val(const int &x,const int &y,const int &z) const {
    return cV(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

  ScalarType &Val(const int &x,const int &y,const int &z) {
    return V(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

  VOX_TYPE &V(const int &x,const int &y,const int &z) {
    return Vol[x+y*this->siz[0]+z*this->siz[0]*this->siz[1]];
  }

  VOX_TYPE &V(const Point3i &pi) {
    return Vol[ pi[0] + pi[1]*this->siz[0] + pi[2]*this->siz[0]*this->siz[1] ];
  }

  const VOX_TYPE &cV(const int &x,const int &y,const int &z) const {
    return Vol[x+y*this->siz[0]+z*this->siz[0]*this->siz[1]];
  }

  bool ValidCell(const Point3i & /*p0*/, const Point3i & /*p1*/) const { return true;}

  template < class VertexPointerType >
  void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,XAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,YAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,ZAxis>(p1,p2,v,thr); }

  /// The following members/methods are just for this particular case.
  /// The above one are the one required by the marching cube interface.

  std::vector<VoxelType> Vol;

  typedef enum { XAxis=0,YAxis=1,ZAxis=2} VolumeAxis;

  template < class VertexPointerType,  VolumeAxis AxisVal >
  void GetIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  {
    float f1 = V(p1).V()-thr;
    float f2 = V(p2).V()-thr;
    float u = (float) f1/(f1-f2);
    if(AxisVal==XAxis) v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
    else v->P().X() = (float) p1.X();
    if(AxisVal==YAxis) v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
    else v->P().Y() = (float) p1.Y();
    if(AxisVal==ZAxis) v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
    else v->P().Z() = (float) p1.Z();
    this->IPfToPf(v->P(),v->P());
    if(VoxelType::HasNormal()) v->N().Import(V(p1).N()*(1-u) + V(p2).N()*u);
  }



  void Init(Point3i _sz, Box3x bb)
  {
    this->siz=_sz;
    this->bbox = bb;
    Vol.resize(this->siz[0]*this->siz[1]*this->siz[2]);
    this->ComputeDimAndVoxel();
  }
};

template <class _ScalarType=float>
class SimpleVoxel
{
private:
  _ScalarType _v;
public:
  typedef _ScalarType ScalarType;
  ScalarType &V() {return _v;}
  ScalarType V() const {return _v;}
  static bool HasNormal() {return false;}
  vcg::Point3<ScalarType> N() const {return Point3<ScalarType>(0,0,0);}
  vcg::Point3<ScalarType> &N()  { static Point3<ScalarType> _p(0,0,0); return _p;}
};

template <class _ScalarType=float>
class SimpleVoxelWithNormal
{
private:
  _ScalarType _v;
  vcg::Point3<_ScalarType> _n;
public:
  typedef _ScalarType ScalarType;
  ScalarType &V() {return _v;}
  ScalarType V() const {return _v;}
  vcg::Point3<ScalarType> &N() {return _n;}
  vcg::Point3<ScalarType> N() const {return _n;}
  static bool HasNormal() {return true;}

};


namespace tri {


// La classe Walker implementa la politica di visita del volume; conoscendo l'ordine di visita del volume
// Ë conveniente che il Walker stesso si faccia carico del caching dei dati utilizzati durante l'esecuzione
// degli algoritmi MarchingCubes ed ExtendedMarchingCubes, in particolare il calcolo del volume ai vertici
// delle celle e delle intersezioni della superficie con le celle. In questo esempio il volume da processare
// viene suddiviso in fette; in questo modo se il volume ha dimensione h*l*w (rispettivamente altezza,
// larghezza e profondit‡), lo spazio richiesto per il caching dei vertici gi‡ allocati passa da O(h*l*w)
// a O(h*l).

template <class MeshType, class VolumeType>
class TrivialWalker
{
private:
    typedef int VertexIndex;
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::VertexPointer VertexPointer;
    public:

  // SetExtractionBox set the portion of the volume to be traversed
  void SetExtractionBox(Box3i subbox)
    {
        _bbox = subbox;
        _slice_dimension = _bbox.DimX()*_bbox.DimZ();

        _x_cs = new VertexIndex[ _slice_dimension ];
        _y_cs = new VertexIndex[ _slice_dimension ];
        _z_cs = new VertexIndex[ _slice_dimension ];
        _x_ns = new VertexIndex[ _slice_dimension ];
        _z_ns = new VertexIndex[ _slice_dimension ];
    }
   
    TrivialWalker()
    { 
      _bbox.SetNull();
      _slice_dimension=0;
    }

    template<class EXTRACTOR_TYPE>
  void BuildMesh(MeshType &mesh, VolumeType &volume, EXTRACTOR_TYPE &extractor, const float threshold, vcg::CallBackPos * cb=0)
  {
    if(_bbox.IsNull() || _slice_dimension==0)
      SetExtractionBox(Box3i(Point3i(0,0,0),volume.ISize()));
    _volume = &volume;
    _mesh		= &mesh;
    _mesh->Clear();
    _thr=threshold;
    Begin();
    extractor.Initialize();
    for (int j=_bbox.min.Y(); j<(_bbox.max.Y()-1)-1; j+=1)
    {
      if(cb && ((j%10)==0) ) 	cb(j*_bbox.DimY()/100.0,"Marching volume");
      for (int i=_bbox.min.X(); i<(_bbox.max.X()-1)-1; i+=1)
      {
        for (int k=_bbox.min.Z(); k<(_bbox.max.Z()-1)-1; k+=1)
        {
          Point3i p1(i,j,k);
          Point3i p2(i+1,j+1,k+1);
          if(volume.ValidCell(p1,p2))
            extractor.ProcessCell(p1, p2);
        }
      }
      NextYSlice();
    }
    extractor.Finalize();
    _volume = NULL;
    _mesh		= NULL;
  }

    float V(int pi, int pj, int pk)
    {
    return _volume->Val(pi, pj, pk)-_thr;
    }

    bool Exist(const vcg::Point3i &p0, const vcg::Point3i &p1, VertexPointer &v)
    {
        int pos = p0.X()+p0.Z()*_bbox.DimX();
        int vidx;

        if (p0.X()!=p1.X()) // punti allineati lungo l'asse X
            vidx = (p0.Y()==_current_slice) ? _x_cs[pos] : _x_ns[pos];
        else if (p0.Y()!=p1.Y()) // punti allineati lungo l'asse Y
            vidx = _y_cs[pos];
        else if (p0.Z()!=p1.Z()) // punti allineati lungo l'asse Z
            vidx = (p0.Y()==_current_slice)? _z_cs[pos] : _z_ns[pos];
        else
            assert(false);

        v = (vidx!=-1)? &_mesh->vert[vidx] : NULL;
        return v!=NULL;
    }

    void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
    {
        int i = p1.X() - _bbox.min.X();
        int z = p1.Z() - _bbox.min.Z();
        VertexIndex index = i+z*_bbox.DimX();
        VertexIndex pos=-1;
        if (p1.Y()==_current_slice)
        {
            if ((pos=_x_cs[index])==-1)
            {
                _x_cs[index] = (VertexIndex) _mesh->vert.size();
                pos = _x_cs[index];
                Allocator<MeshType>::AddVertices( *_mesh, 1 );
                v = &_mesh->vert[pos];
                _volume->GetXIntercept(p1, p2, v, _thr);
                return;
            }
        }
        if (p1.Y()==_current_slice+1)
        {
            if ((pos=_x_ns[index])==-1)
            {
                _x_ns[index] = (VertexIndex) _mesh->vert.size();
                pos = _x_ns[index];
                Allocator<MeshType>::AddVertices( *_mesh, 1 );
                v = &_mesh->vert[pos];
                _volume->GetXIntercept(p1, p2, v,_thr);
                return;
            }
        }
    assert(pos >=0 && size_t(pos)< _mesh->vert.size());
        v = &_mesh->vert[pos];
    }
    void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
    {
        int i = p1.X() - _bbox.min.X();
        int z = p1.Z() - _bbox.min.Z();
        VertexIndex index = i+z*_bbox.DimX();
        VertexIndex pos;
        if ((pos=_y_cs[index])==-1)
        {
            _y_cs[index] = (VertexIndex) _mesh->vert.size();
            pos = _y_cs[index];
            Allocator<MeshType>::AddVertices( *_mesh, 1);
            v = &_mesh->vert[ pos ];
            _volume->GetYIntercept(p1, p2, v,_thr);
        }
        v = &_mesh->vert[pos];
    }
    void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v)
    {
        int i = p1.X() - _bbox.min.X();
        int z = p1.Z() - _bbox.min.Z();
        VertexIndex index = i+z*_bbox.DimX();
        VertexIndex pos;
        if (p1.Y()==_current_slice)
        {
            if ((pos=_z_cs[index])==-1)
            {
                _z_cs[index] = (VertexIndex) _mesh->vert.size();
                pos = _z_cs[index];
                Allocator<MeshType>::AddVertices( *_mesh, 1 );
                v = &_mesh->vert[pos];
                _volume->GetZIntercept(p1, p2, v,_thr);
                return;
            }
        }
        if (p1.Y()==_current_slice+1)
        {
            if ((pos=_z_ns[index])==-1)
            {
                _z_ns[index] = (VertexIndex) _mesh->vert.size();
                pos = _z_ns[index];
                Allocator<MeshType>::AddVertices( *_mesh, 1 );
                v = &_mesh->vert[pos];
                _volume->GetZIntercept(p1, p2, v,_thr);
                return;
            }
        }
        v = &_mesh->vert[pos];
    }

protected:
    Box3i		_bbox;

    int _slice_dimension;
    int	_current_slice;

    VertexIndex *_x_cs; // indici dell'intersezioni della superficie lungo gli Xedge della fetta corrente
    VertexIndex	*_y_cs; // indici dell'intersezioni della superficie lungo gli Yedge della fetta corrente
    VertexIndex *_z_cs; // indici dell'intersezioni della superficie lungo gli Zedge della fetta corrente
    VertexIndex *_x_ns; // indici dell'intersezioni della superficie lungo gli Xedge della prossima fetta
    VertexIndex *_z_ns; // indici dell'intersezioni della superficie lungo gli Zedge della prossima fetta

    MeshType		*_mesh;
    VolumeType	*_volume;

  float _thr;
    void NextYSlice()
    {
        memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
        memset(_y_cs,	-1, _slice_dimension*sizeof(VertexIndex));
        memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));

        std::swap(_x_cs, _x_ns);
        std::swap(_z_cs, _z_ns);

        _current_slice += 1;
    }

    void Begin()
    {
        _current_slice = _bbox.min.Y();

        memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
        memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
        memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
        memset(_x_ns, -1, _slice_dimension*sizeof(VertexIndex));
        memset(_z_ns, -1, _slice_dimension*sizeof(VertexIndex));

    }
};
} // end namespace tri
} // end namespace vcg
#endif // __VCGTEST_WALKER
