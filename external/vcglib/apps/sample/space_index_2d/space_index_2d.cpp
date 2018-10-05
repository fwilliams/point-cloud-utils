/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2009                                           \/)\/    *
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
#include <stdio.h>
#include <time.h>
#include <vcg/space/distance2.h>
#include<vcg/space/segment2.h>
#include<vcg/space/index/grid_static_ptr2d.h>
#include<vcg/space/index/grid_closest2d.h>
#include<vcg/space/intersection2.h>

typedef double MyScalarType;
typedef vcg::Point2<MyScalarType> MyCoordType;
typedef vcg::Ray2<MyScalarType> MyRayType;

//**BASIC SEGMENT CLASS
class MySegmentType:public vcg::Segment2<MyScalarType>
{
public:
	int mark;
    bool IsD(){return false;}

    MySegmentType(const vcg::Point2<MyScalarType> &_P0,
                  const vcg::Point2<MyScalarType> &_P1)
    {
        P0()=_P0;
        P1()=_P1;
        mark=0;
    }

    int &TMark(){return mark;}

	MySegmentType(){}

    MySegmentType(const MySegmentType &s1):vcg::Segment2<MyScalarType>(s1)
	{
		P0()=s1.P0();
		P1()=s1.P1();
		mark=s1.mark;
	}
};


//**ALLOCATED SEGMENTS**//
std::vector<MySegmentType> AllocatedSeg;

//**GENERATION OF RANDOM SEGMENTS

vcg::Point2<MyScalarType> RandomPoint(MyScalarType SpaceSize=100)
{
    int dimension=RAND_MAX;
    int X=rand();
    int Y=rand();
    vcg::Point2<MyScalarType> P0=vcg::Point2<MyScalarType>((MyScalarType)X/dimension,(MyScalarType)Y/dimension);
    P0*=SpaceSize;
    return P0;
}


void RandomSeg(vcg::Point2<MyScalarType> &P0,
                vcg::Point2<MyScalarType> &P1,
                MyScalarType SpaceSize=100,
                MyScalarType maxdim=1)
{

    P0=RandomPoint(SpaceSize);
    vcg::Point2<MyScalarType> D=RandomPoint(SpaceSize);
    D.Normalize();
    D*=maxdim;
    P1=P0+D;
}

void InitRandom(int num,
                MyScalarType SpaceSize=100,
                MyScalarType maxdim=1)
{
    AllocatedSeg.clear();
    AllocatedSeg.resize(num);
    srand(clock());
    for (int i=0;i<num;i++)
    {
        vcg::Point2<MyScalarType> P0,P1;
        RandomSeg(P0,P1,SpaceSize,maxdim);
        AllocatedSeg[i]=MySegmentType(P0,P1);
        //AllocatedSeg[i].deleted=false;
    }

}


//**MARKER CLASSES**//
class MyMarker
{
	
public:
	int mark;

	MyMarker(){mark=0;}

    void UnMarkAll(){mark++;}

	bool IsMarked(MySegmentType* obj)
    {
        int markObj=obj->TMark();
        return(markObj==mark);
    }

	void Mark(MySegmentType* obj)
    {obj->TMark()=mark;}

};

//**GRID-RELATED STUFF**//
MyMarker mf;
vcg::GridStaticPtr2D<MySegmentType,MyScalarType> Grid2D;

//**QUERIES
MySegmentType * GetClosestSegment(MyCoordType & _p,
                                  MyCoordType &_closestPt)
{
    vcg::PointSegment2DEPFunctor<MyScalarType> PDistFunct;

    MyScalarType _minDist;
    MyScalarType _maxDist=std::numeric_limits<MyScalarType>::max();
    return (Grid2D.GetClosest(PDistFunct,mf,_p,_maxDist,_minDist,_closestPt));
}

void GetInBoxSegments(vcg::Box2<MyScalarType> bbox,std::vector<MySegmentType*> &result)
{
    Grid2D.GetInBox(mf,bbox,result);
}

MySegmentType * DoRay(MyRayType & _r,
                      MyCoordType &_closestPt)
{
    MyRayType _ray1=_r;
    _ray1.Normalize();
    typedef vcg::RaySegmentIntersectionFunctor SintFunct;
    SintFunct rs;
    MyScalarType _maxDist=std::numeric_limits<MyScalarType>::max();
    MyScalarType _t;
    MySegmentType *seg=Grid2D.DoRay(rs,mf,_ray1,_maxDist,_t);

    if (seg==NULL)return NULL;
        _closestPt=_ray1.Origin()+_ray1.Direction()*_t;
    return seg;
}

//**BRUTE FORCE QUERIES
void GetInBoxSegmentsBruteF( vcg::Box2<MyScalarType> bbox,
                             std::vector<MySegmentType*> &result)
{
  vcg::Box2<MyScalarType> ibbox;
  for (size_t i=0;i<AllocatedSeg.size();i++)
  {
    AllocatedSeg[i].GetBBox(ibbox);
    if (!ibbox.Collide(bbox)) continue;
    result.push_back(&AllocatedSeg[i]);
  }
}

MySegmentType* GetClosesestSegmentBruteF(MyCoordType & _p,
                                        MyCoordType &_closestPt)
{
    MyScalarType _minDist=std::numeric_limits<MyScalarType>::max();
    MySegmentType *ret=NULL;
    for (size_t i=0;i<AllocatedSeg.size();i++)
    {
        vcg::Point2<MyScalarType> test;
        test=vcg::ClosestPoint(AllocatedSeg[i],_p);
        MyScalarType currD=(test-_p).Norm();
        if (currD<_minDist)
        {
            _closestPt=test;
            _minDist=currD;
            ret=&AllocatedSeg[i];
        }
    }
    return ret;
}

MySegmentType * DoRayBruteF(MyRayType & _r,
                            MyCoordType &_closestPt)
{
    MyScalarType _minDist=std::numeric_limits<MyScalarType>::max();
    MySegmentType *ret=NULL;
    for (size_t i=0;i<AllocatedSeg.size();i++)
    {
        vcg::Point2<MyScalarType> test;
        bool inters=vcg::RaySegmentIntersection(_r,AllocatedSeg[i],test);
        if (!inters)continue;
        MyScalarType currD=(test-_r.Origin()).Norm();
        if (currD<_minDist)
        {
            _closestPt=test;
            _minDist=currD;
            ret=&AllocatedSeg[i];
        }
    }
    return ret;
}


void TestBox(int num_test=100000,
             MyScalarType SpaceSize=100)
{
    int numWrong=0;

    for (int i=0;i<num_test;i++)
    {
        vcg::Point2<MyScalarType> P0=RandomPoint(SpaceSize);
        vcg::Point2<MyScalarType> P1=RandomPoint(SpaceSize);
        vcg::Box2<MyScalarType> bbox;
        bbox.Add(P0);
        bbox.Add(P1);
        std::vector<MySegmentType*> result0;
        GetInBoxSegments(bbox,result0);

        std::vector<MySegmentType*> result1;
        GetInBoxSegmentsBruteF(bbox,result1);

        std::sort(result0.begin(),result0.end());
        std::sort(result1.begin(),result1.end());

        std::vector<MySegmentType*>::iterator new_end=std::unique(result1.begin(),result1.end());
        int dist=distance(result1.begin(),new_end);
        result1.resize(dist);

        if (result0.size()!=result1.size())numWrong++;

        for (size_t j = 0; j < result0.size(); j++)
            if (result0[j] != result1[j])
            {
                numWrong++;
            }
    }
    printf("WRONG TESTS BBOX %d ON %d \n",numWrong,num_test);
    fflush(stdout);
}

void TestClosest(int num_test=100000,
                  MyScalarType SpaceSize=100)
{

    int numWrong=0;
    for (int i=0;i<num_test;i++)
    {
        vcg::Point2<MyScalarType> P0=RandomPoint(SpaceSize);

        vcg::Point2<MyScalarType> closest0;
        MySegmentType* result0=GetClosestSegment(P0,closest0);

        vcg::Point2<MyScalarType> closest1;
        MySegmentType* result1=GetClosesestSegmentBruteF(P0,closest1);


        if (result0!=result1)
        {
            numWrong++;
            printf("D0 %5.5f  \n",(closest0-P0).Norm());
            printf("D1 %5.5f  \n",(closest1-P0).Norm());
            fflush(stdout);
        }
    }
    printf("WRONG TESTS CLOSEST %d ON %d \n",numWrong,num_test);
    fflush(stdout);
}

void TestRay(int num_test=100000,
             MyScalarType SpaceSize=100)
{
    int numWrong=0;
    int NUll0=0;
    int NUll1=0;
    for (int i=0;i<num_test;i++)
    {
        vcg::Point2<MyScalarType> P0=RandomPoint(SpaceSize);
        vcg::Point2<MyScalarType> P1=RandomPoint(SpaceSize);
        vcg::Point2<MyScalarType> Orig=P0;
        vcg::Point2<MyScalarType> Dir=P1-P0;
        Dir.Normalize();

        MyRayType r(Orig,Dir);

        vcg::Point2<MyScalarType> closest0;
        MySegmentType* result0=DoRay(r,closest0);

        vcg::Point2<MyScalarType> closest1;
        MySegmentType* result1=DoRayBruteF(r,closest1);


        if (result0!=result1)
        {
            numWrong++;
//            printf("D0 %5.5f  \n",(closest0-P0).Norm());
//            printf("D1 %5.5f  \n",(closest1-P0).Norm());
//            fflush(stdout);
        }
        if (result0==NULL) NUll0++;
        if (result1==NULL) NUll1++;
    }
    printf("WRONG TESTS DORAY %d ON %d \n",numWrong,num_test);
    printf("NULL0 %d \n",NUll0);
    printf("NULL1 %d \n",NUll1);
    fflush(stdout);
}

int main( int argc, char **argv )
{
  (void) argc;
  (void) argv;
  int num_sample=20000;
  int t0=clock();

  printf("** Random Initialization ** \n");
  fflush(stdout);
  InitRandom(num_sample,100,0.3);
  int t1=clock();

  ///Initialization performance
  printf("** Time elapsed for initialization of %d sample is %d\n \n",num_sample,t1-t0);
  Grid2D.Set(AllocatedSeg.begin(),AllocatedSeg.end());
  fflush(stdout);

  //Box Query correctness
  TestBox(num_sample);
  TestClosest(num_sample);
  TestRay(num_sample);

  return 0;
}
