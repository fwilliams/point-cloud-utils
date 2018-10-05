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
#include <QtOpenGL/QtOpenGL>
#include<vcg/space/box2.h>
#include<vcg/space/box3.h>
#include<vcg/math/random_generator.h>
#include<wrap/qt/col_qt_convert.h>
#include <vcg/space/rect_packer.h>
#include <vcg/space/outline2_packer.h>
#include <vcg/space/rasterized_outline2_packer.h>
#include <vcg/complex/algorithms/outline_support.h>
#include <wrap/qt/Outline2ToQImage.h>
#include <wrap/qt/outline2_rasterizer.h>

#include <time.h>

using namespace vcg;
using namespace std;

static void buildRandRectSet(int rectNum, vector<Box2f> &rectVec)
{
  math::MarsenneTwisterRNG rnd;
  float exp=3.0f;
  float ratioMin=0.2;
  float ratioMax=0.9;
  float sizeMin=0.1;
  float sizeMax=1.0f;
  rnd.initialize(time(0));
  for(int i=0;i<rectNum;++i)
  {
    Box2f bb;
    float ratio=ratioMin+(ratioMax-ratioMin)*rnd.generate01();
    float size= sizeMin+(sizeMax-sizeMin)*pow((float)rnd.generate01(),exp);
    bb.min=Point2f(-size*ratio,-size);
    bb.max=Point2f( size*ratio, size);
    rectVec.push_back(bb);
  }
}

int main( int argc, char **argv )
{
  vector<Similarity2f> trVec;
  vector<Similarity2f> trPolyVec;
  vector< vector<Point2f> > outline2Vec;
  vector< vector<Point2f> > multiPolySet;
  Point2f finalSize;
  std::vector<Point2f> finalSizeVec;
  const Point2i containerSize(1024,1024);
  Outline2Dumper::Param pp;
  std::vector<int> contInd;

  vector<Box2f> rectVec;
  buildRandRectSet(10,rectVec);
//  RectPacker<float>::Pack(rectVec,containerSize,trVec,finalSize);
  RectPacker<float>::PackMulti(rectVec,containerSize,3,trVec,contInd,finalSizeVec);
  RectPacker<float>::Stat s = RectPacker<float>::stat();
  printf("RectPacker attempt %i time %5.3f %5.3f\n",s.pack_attempt_num,s.pack_total_time,s.pack_attempt_time);

//  Outline2Dumper::rectSetToPolySet(rectVec,polySet);

//  Outline2Dumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,0,polySet,trPolyVec);
//  Outline2Dumper::dumpPolySetPNG("testpolyEq0.png",polySet,trPolyVec,pp);
//  Outline2Dumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,1,polySet,trPolyVec);
//  Outline2Dumper::dumpPolySetPNG("testpolyEq1.png",polySet,trPolyVec,pp);
//  Outline2Dumper::multiRectSetToSinglePolySet(rectVec,trVec,contInd,2,polySet,trPolyVec);
//  Outline2Dumper::dumpPolySetPNG("testpolyEq2.png",polySet,trPolyVec,pp);


//   buildRandPolySet(100,polySet);
//   PolyPacker<float>::PackMultiAsObjectOrientedRect(polySet,containerSize,3,trVec,contInd,finalSizeVec);

//   Outline2Dumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,0,multiPolySet,trPolyVec);
//   Outline2Dumper::dumpPolySetPNG("testpolyEq0.png",multiPolySet,trPolyVec,pp);

//   Outline2Dumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,1,multiPolySet,trPolyVec);
//   Outline2Dumper::dumpPolySetPNG("testpolyEq1.png",multiPolySet,trPolyVec,pp);

//   Outline2Dumper::multiPolySetToSinglePolySet(polySet,trVec,contInd,2,multiPolySet,trPolyVec);
//   Outline2Dumper::dumpPolySetPNG("testpolyEq2.png",multiPolySet,trPolyVec,pp);

   //  Outline2Dumper::dumpPolySetPNG("testpolyOO.png",polySet,trVec,pp);

  vcg::tri::OutlineUtil<float>::BuildRandomOutlineVec(25,outline2Vec);

  PolyPacker<float>::PackAsEqualSquares(outline2Vec,containerSize,trVec,finalSize);
  Outline2Dumper::dumpOutline2VecPNG("testpolyEq.png",outline2Vec,trVec,pp);

  PolyPacker<float>::PackAsAxisAlignedRect(outline2Vec,containerSize,trVec,finalSize);
  Outline2Dumper::dumpOutline2VecPNG("testpolyAA.png",outline2Vec,trVec,pp);

  PolyPacker<float>::PackAsObjectOrientedRect(outline2Vec,containerSize,trVec,finalSize);
  Outline2Dumper::dumpOutline2VecPNG("testpolyOO.png",outline2Vec,trVec,pp);

  RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters packingParam;

  packingParam.costFunction  = RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Parameters::LowestHorizon;
  packingParam.doubleHorizon = true;
  packingParam.cellSize = 4;
  packingParam.rotationNum = 16; //number of rasterizations in 90°

  RasterizedOutline2Packer<float, QtOutline2Rasterizer>::Pack(outline2Vec,containerSize,trVec,packingParam);
  Outline2Dumper::dumpOutline2VecPNG("testpolyRR.png",outline2Vec,trVec,pp);


  return 0;
}
