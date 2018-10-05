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
#ifndef __VOXEL_H__
#define __VOXEL_H__

namespace vcg {

// Stato di un voxel

// B() dice se ci sono dati in uno stadio usabile.
// Cnt() dice quanti ce ne sono stati sommati (per la normalizzazione)

// b==false cnt==0 totalmente non inzializzato (Zero)
// b==false cnt >0 da normalizzare
// b==true  cnt==0 gia' normalizzato
// b==true  cnt >0 Errore!!!


/**
 *
 */

template<class SCALAR_TYPE=float>
class Voxel
{
public:
  typedef SCALAR_TYPE scalar;

  Voxel(SCALAR_TYPE vv, bool bb, Point3<scalar> nn, short _cnt) {v=vv;b=bb;n=nn;cnt=_cnt;}
  Voxel(SCALAR_TYPE vv, Point3<scalar> nn, scalar qq) {v=vv;b=true;n=nn;cnt=0;q=qq;}

  const Point3<scalar> &N() const  	{  return n; 		}

  void SetN(const Point3<scalar> &nn) 	{ n=nn;		}
  const scalar &V() const 		{ 	 return v; 		}

  void SetV(const scalar &vv) 		{		 v=vv;		}

  const scalar &Q() const		{		 return q;		}

  void SetQ(const scalar &qq) 		{		 q=qq;		}


  bool B() const {return b;};
  void SetB(bool val) {b=val;}
  int Cnt() const {return cnt;}
  void SetCnt(int val) {cnt=val;}
  inline void Blend( Voxel const & vx, scalar w)
  {
    float w1=1.0-w;
    v=v*w1+vx.v*w;
    q=q*w1+vx.q*w;
    n=n*w1+vx.n*w;
    //return *this;
  }

  inline Voxel & operator += ( Voxel const & vx)
  {
    assert(!b);
    if(cnt==0)
    {
      v=vx.v;
      q=vx.q;
      n=vx.n;
      cnt=1;
      b=false;
    }
    else
    {
      v+=vx.v;
      q+=vx.q;
      n+=vx.n;
      ++cnt;
    }
    return *this;
  }

  inline bool Normalize(int thr)
  {
    assert(cnt>0);
    assert(!B());
    if(cnt<thr)
    {
      (*this) = Zero();
      return false;
    }
    v/=cnt;
    q/=cnt;
    n/=cnt;
    cnt=0;
    b=true;
    return true;
  }

  static const Voxel &Zero() {
    static Voxel tt(0,false,Point3f(0,0,0),0);
    return tt;
  }
  void Merge(const Voxel &VOX)
  {
    v=(v*q+VOX.Q()*VOX.v)/(q+VOX.Q());
    n=(n*q+VOX.n*VOX.Q())/(q+VOX.Q());
    q=q+VOX.Q();
  }

  void Set(const Voxel &VOX)
  {
    v=VOX.v;
    n=VOX.n;
    q=VOX.q;
  }

protected:
  bool b;
  short cnt;
  scalar v;
  scalar q;
  Point3<SCALAR_TYPE> n;

};


class Voxelfc :public Voxel<float>
{
public:

  Voxelfc(float vv, bool bb, Point3f nn, short _cnt) :Voxel<float>(vv,bb,nn,_cnt){}
  Voxelfc(float vv, Point3f nn, scalar qq) :Voxel<float>(vv,nn,qq) {}
  Voxelfc(float vv, Point3f nn, scalar qq,Color4b cc) :Voxel<float>(vv,nn,qq)
  {
    c[0]=cc[0];
    c[1]=cc[1];
    c[2]=cc[2];
  }

  inline bool Normalize(int thr)
  {
    if(cnt>=thr)  c/=cnt;
    return Voxel<float>::Normalize(thr);
  }

  static const Voxelfc &Zero() {
    static Voxelfc tt(0,false,Point3f(0,0,0),0);
    return tt;
  }

  void Merge(const Voxelfc &VOX)
  {
    c=( c*q + VOX.C()*VOX.Q() )/(q+VOX.Q());
    Voxel<float>::Merge(VOX);
  }

  void Set(const Voxelfc &VOX)
  {
    Voxel<float>::Set(VOX);
    c=VOX.c;
  }

  const float &C(const int i) const { return c[i]; }
  const Point3f &C() const  	{  return c; 		}
  void SetC(const Point3f &cc) 	{ c=cc;		}
  Color4b C4b() const
  {
    static Color4b cc;
    cc=Color4b(c[0],c[1],c[2],255);
    return cc;
  }
  inline void Blend( Voxelfc const & vx, scalar w)
  {
    float w1=1.0-w;
    v=v*w1+vx.v*w;
    q=q*w1+vx.q*w;
    n=n*w1+vx.n*w;
    c=c*w1+vx.c*w;
    //return *this;
  }

  inline Voxelfc & operator += ( Voxelfc const & vx)
  {
    Voxel<float>::operator +=(vx);
    if(cnt==1)	c =vx.c;
    else c+=vx.c;
    return *this;
  }

private:
  Point3f c;
};
}
#endif
