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

$Log: gen_normal.h,v $
****************************************************************************/

#ifndef __VCG_GEN_NORMAL
#define __VCG_GEN_NORMAL

#include <algorithm>

namespace vcg {

template <class ScalarType>
class GenNormal
{
public:
typedef Point3<ScalarType> Point3x;

static void  Random(int vn, std::vector<Point3<ScalarType > > &NN)
{
  NN.clear();
  while(NN.size()<vn)
  {
    Point3x pp(((float)rand())/RAND_MAX,
           ((float)rand())/RAND_MAX,
           ((float)rand())/RAND_MAX);
    pp=pp*2.0-Point3x(1,1,1);
    if(pp.SquaredNorm()<1)
    {
      Normalize(pp);
      NN.push_back(pp);
    }
  }
}


static Point3x FibonacciPt(int i, int n)
{
  const ScalarType Phi =  ScalarType(std::sqrt(ScalarType(5))*0.5 + 0.5);
  const ScalarType phi = 2.0*M_PI* (i/Phi - floor(i/Phi));
  ScalarType cosTheta = 1.0 - (2*i + 1.0)/ScalarType(n);
    ScalarType sinTheta = 1 - cosTheta*cosTheta;
    sinTheta = std::sqrt(std::min(ScalarType(1),std::max(ScalarType(0),sinTheta)));
    return Point3x(
      cos(phi)*sinTheta,
      sin(phi)*sinTheta,
      cosTheta);
}

// Implementation of the Spherical Fibonacci Point Sets
// according to the description of
// Spherical Fibonacci Mapping
// Benjamin Keinert, Matthias Innmann, Michael Sanger, Marc Stamminger
// TOG 2015
static void Fibonacci(int n, std::vector<Point3x > &NN)
{
  NN.resize(n);
  for(int i=0;i<n;++i)
    NN[i]=FibonacciPt(i,n);
}

static void UniformCone(int vn, std::vector<Point3<ScalarType > > &NN, ScalarType AngleRad, Point3x dir=Point3x(0,1,0))
{
  std::vector<Point3<ScalarType > > NNT;
  NN.clear();
  // per prima cosa si calcola il volume della spherical cap di angolo AngleRad
  ScalarType Height= 1.0 - cos(AngleRad); // height is measured from top...
  // Surface is the one of the tangent cylinder
  ScalarType CapArea = 2.0*M_PI*Height;
  ScalarType Ratio = CapArea / (4.0*M_PI );

  printf("----------AngleRad %f Angledeg %f ratio %f vn %i vn2 %i \n",AngleRad,math::ToDeg(AngleRad),Ratio,vn,int(vn/Ratio));
  Fibonacci(vn/Ratio,NNT);
  printf("asked %i got %i (expecting %i instead of %i)\n", int(vn/Ratio), NNT.size(), int(NNT.size()*Ratio), vn);
  typename std::vector<Point3<ScalarType> >::iterator vi;

  ScalarType cosAngle = cos(AngleRad);
  for(vi=NNT.begin();vi!=NNT.end();++vi)
  {
    if(dir.dot(*vi) >= cosAngle) NN.push_back(*vi);
  }
 }

// This is an Implementation of the Dave Rusin’s Disco Ball algorithm
// You can spread the points as follows:
// Put  N+1  points on the meridian from north to south poles, equally spaced.
// If you swing this meridian around the sphere, you'll sweep out the entire
// surface; in the process, each of the points will sweep out a circle. You
// can show that the  ith  point will sweep out a circle of radius sin(pi i/N).
// If you space points equally far apart on this circle, keeping the
// displacement roughly the same as on that original meridian, you'll be
// able to fit about  2N sin(pi i/N) points here. This process will put points
// pretty evenly spaced on the sphere; the number of such points is about
//     2+ 2N*Sum(i=1 to N-1) sin(pi i/N).
// The closed form of this summation
//     2.0 - ( (2.0*N * sin (M_PI/N))/(cos(M_PI/N) - 1.0));
static void DiscoBall(int vn, std::vector<Point3<ScalarType > > &NN)
{
  // Guess the right N
  ScalarType N=0;

  for(N=1;N<vn;++N)
  {
    ScalarType expectedPoints = 2.0 - ( (2.0*N * sin (M_PI/N))/(cos(M_PI/N) - 1.0));
    if(expectedPoints >= vn) break;
  }

  ScalarType VerticalAngle = M_PI / N;
  NN.push_back(Point3<ScalarType>(0,0,1.0));
  for (int i =1; i<N; ++i)
  {
    // Z is the north/south axis
    ScalarType HorizRadius = sin(i*VerticalAngle);
    ScalarType CircleLength = 2.0 * M_PI * HorizRadius;
    ScalarType Z = cos(i*VerticalAngle);
    ScalarType PointNumPerCircle = floor( CircleLength / VerticalAngle);
    ScalarType HorizontalAngle = 2.0*M_PI/PointNumPerCircle;
    for(ScalarType j=0;j<PointNumPerCircle;++j)
    {
      ScalarType X = cos(j*HorizontalAngle)*HorizRadius;
      ScalarType Y = sin(j*HorizontalAngle)*HorizRadius;
      NN.push_back(Point3<ScalarType>(X,Y,Z));
    }
  }
  NN.push_back(Point3<ScalarType>(0,0,-1.0));
}

static void RecursiveOctahedron(int vn, std::vector<Point3<ScalarType > > &NN)
{
  OctaLevel pp;

  int ll=10;
  while(pow(4.0f,ll)+2>vn) ll--;

  pp.Init(ll);
  sort(pp.v.begin(),pp.v.end());
  int newsize = unique(pp.v.begin(),pp.v.end())-pp.v.begin();
  pp.v.resize(newsize);

  NN=pp.v;
  //Perturb(NN);
 }

static void Perturb(std::vector<Point3<ScalarType > > &NN)
{
  float width=0.2f/sqrt(float(NN.size()));

  typename std::vector<Point3<ScalarType> >::iterator vi;
  for(vi=NN.begin(); vi!=NN.end();++vi)
  {
    Point3x pp(((float)rand())/RAND_MAX,
           ((float)rand())/RAND_MAX,
           ((float)rand())/RAND_MAX);
    pp=pp*2.0-Point3x(1,1,1);
    pp*=width;
    (*vi)+=pp;
    (*vi).Normalize();
  }

}

/*
Trova la normale piu vicina a quella data.
Assume che tutte normale in ingresso sia normalizzata;
*/
static int BestMatchingNormal(const Point3x &n, std::vector<Point3x> &nv)
{
    int ret=-1;
    ScalarType bestang=-1;
    ScalarType cosang;
    typename std::vector<Point3x>::iterator ni;
    for(ni=nv.begin();ni!=nv.end();++ni)
        {
            cosang=(*ni).dot(n);
            if(cosang>bestang) {
                bestang=cosang;
                ret=ni-nv.begin();
            }
        }
  assert(ret>=0 && ret <int(nv.size()));
    return ret;
}


private :
class OctaLevel
  {
  public:
    std::vector<Point3x> v;
    int level;
    int sz;
    int sz2;

    Point3x &Val(int i, int j) {

      assert(i>=-sz2 && i<=sz2);
      assert(j>=-sz2 && j<=sz2);
      return v[i+sz2 +(j+sz2)*sz];
    }
/*
 *  Only the first quadrant is generated and replicated onto the other ones.
 *
 *   o         lev == 1
 *   | \       sz2 = 2^lev = 2
 *   o - o     sz = 5 (eg. all the points lie in a 5x5 squre)
 *   | \ | \
 *   o - o - o
 *
 *     |
 *     V
 *
 *   o
 *   | \       lev == 1
 *   o - o     sz2 = 4
 *   | \ | \   sz = 9 (eg. all the points lie in a 9x9 squre)
 *   o - o - o
 *   | \ | \ | \
 *   o - o - o - o
 *   | \ | \ | \ | \
 *   o - o - o - o - o
 *
 *
 */
    void Init(int lev)
    {
      sz2=pow(2.0f,lev);
      sz=sz2*2+1;
      v.resize(sz*sz,Point3x(0,0,0));
      if(lev==0)
      {
        Val( 0,0)=Point3x( 0, 0, 1);
        Val( 1,0)=Point3x( 1, 0, 0);
        Val( 0,1)=Point3x( 0, 1, 0);
      }
      else
      {
        OctaLevel tmp;
        tmp.Init(lev-1);
        int i,j;
        for(i=0;i<=sz2;++i)
          for(j=0;j<=(sz2-i);++j)
          {
            if((i%2)==0 && (j%2)==0)
              Val(i,j)=tmp.Val(i/2,j/2);
            if((i%2)!=0 && (j%2)==0)
                Val(i,j)=(tmp.Val((i-1)/2,j/2)+tmp.Val((i+1)/2,j/2))/2.0;
            if((i%2)==0 && (j%2)!=0)
                Val(i,j)=(tmp.Val(i/2,(j-1)/2)+tmp.Val(i/2,(j+1)/2))/2.0;
            if((i%2)!=0 && (j%2)!=0)
                 Val(i,j)=(tmp.Val((i-1)/2,(j+1)/2)+tmp.Val((i+1)/2,(j-1)/2))/2.0;

            Val( sz2-j, sz2-i)[0] = Val(i,j)[0]; Val( sz2-j, sz2-i)[1] = Val(i,j)[1]; Val( sz2-j, sz2-i)[2] = -Val(i,j)[2];
            Val(-sz2+j, sz2-i)[0] =-Val(i,j)[0]; Val(-sz2+j, sz2-i)[1] = Val(i,j)[1]; Val(-sz2+j, sz2-i)[2] = -Val(i,j)[2];
            Val( sz2-j,-sz2+i)[0] = Val(i,j)[0]; Val( sz2-j,-sz2+i)[1] =-Val(i,j)[1]; Val( sz2-j,-sz2+i)[2] = -Val(i,j)[2];
            Val(-sz2+j,-sz2+i)[0] =-Val(i,j)[0]; Val(-sz2+j,-sz2+i)[1] =-Val(i,j)[1]; Val(-sz2+j,-sz2+i)[2] = -Val(i,j)[2];

            Val(-i,-j)[0] = -Val(i,j)[0]; Val(-i,-j)[1] = -Val(i,j)[1]; Val(-i,-j)[2] = Val(i,j)[2];
            Val( i,-j)[0] =  Val(i,j)[0]; Val( i,-j)[1] = -Val(i,j)[1]; Val( i,-j)[2] = Val(i,j)[2];
            Val(-i, j)[0] = -Val(i,j)[0]; Val(-i, j)[1] =  Val(i,j)[1]; Val(-i, j)[2] = Val(i,j)[2];
          }

        typename std::vector<Point3<ScalarType> >::iterator vi;
        for(vi=v.begin(); vi!=v.end();++vi)
            (*vi).Normalize();
      }
    }
  };
};
}
#endif
