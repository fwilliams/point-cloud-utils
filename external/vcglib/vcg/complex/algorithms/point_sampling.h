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

The sampling Class has a set of static functions, that you can call to sample the surface of a mesh.
Each function is templated on the mesh and on a Sampler object s.
Each function calls many time the sample object with the sampling point as parameter.

Sampler Classes and Sampling algorithms are independent.
Sampler classes exploits the sample that are generated with various algorithms.
For example, you can compute Hausdorff distance (that is a sampler) using various
sampling strategies (montecarlo, stratified etc).

****************************************************************************/
#ifndef __VCGLIB_POINT_SAMPLING
#define __VCGLIB_POINT_SAMPLING


#include <vcg/math/random_generator.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/space/segment2.h>
#include <vcg/space/index/grid_static_ptr.h>
namespace vcg
{
namespace tri
{
/// \ingroup trimesh
/// \headerfile point_sampling.h vcg/complex/algorithms/point_sampling.h

/**
 \brief A basic sampler class that show the required interface used by the SurfaceSampling class.

 Most of the methods of sampling classes call the AddFace method of this class with the face containing the sample and its barycentric coord.
 Beside being an example of how to write a sampler it provides a simple way to use the various sampling classes.
 For example if you just want to get a vector with positions over the surface You have just to write

     vector<Point3f> myVec;
     TrivialSampler<MyMesh> ts(myVec);
     SurfaceSampling<MyMesh, TrivialSampler<MyMesh> >::Montecarlo(M, ts, SampleNum);

**/

template <class MeshType>
class TrivialSampler
{
public:
  typedef typename MeshType::ScalarType  ScalarType;
  typedef typename MeshType::CoordType  CoordType;
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::EdgeType   EdgeType;
  typedef typename MeshType::FaceType   FaceType;

  void reset()
  {
    sampleVec->clear();
  }

  TrivialSampler()
  {
    sampleVec = new std::vector<CoordType>();
    vectorOwner=true;
  }

  TrivialSampler(std::vector<CoordType> &Vec)
  {
    sampleVec = &Vec;
    vectorOwner=false;
    reset();
  }

  ~TrivialSampler()
  {
    if(vectorOwner) delete sampleVec;
  }

private:
  std::vector<CoordType> *sampleVec;
  bool vectorOwner;
public:
  
  std::vector<CoordType> &SampleVec() 
  {
    return *sampleVec;
  }

  void AddVert(const VertexType &p)
  {
    sampleVec->push_back(p.cP());
  }
  void AddEdge(const EdgeType& e, ScalarType u ) // u==0 -> v(0) u==1 -> v(1);
  {
    sampleVec->push_back(e.cV(0)->cP()*(1.0-u)+e.cV(1)->cP()*u);
  }

  void AddFace(const FaceType &f, const CoordType &p)
  {
    sampleVec->push_back(f.cP(0)*p[0] + f.cP(1)*p[1] +f.cP(2)*p[2] );
  }

  void AddTextureSample(const FaceType &, const CoordType &, const Point2i &, float )
  {
    // Retrieve the color of the sample from the face f using the barycentric coord p
    // and write that color in a texture image at position <tp[0], texHeight-tp[1]>
    // if edgeDist is > 0 then the corrisponding point is affecting face color even if outside the face area (in texture space)
  }
}; // end class TrivialSampler

template <class MeshType>
class TrivialPointerSampler
{
public:
  typedef typename MeshType::CoordType  CoordType;
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::FaceType   FaceType;

  TrivialPointerSampler() {}
  ~TrivialPointerSampler()  {}

  void reset()
  {
    sampleVec.clear();
  }

public:
  std::vector<VertexType *> sampleVec;

  void AddVert(VertexType &p)
  {
    sampleVec.push_back(&p);
  }
  // This sampler should be used only for getting vertex pointers. Meaningless in other case.
  void AddFace(const FaceType &, const CoordType &)   { assert(0); }
  void AddTextureSample(const FaceType &, const CoordType &, const Point2i &, float ) { assert(0); }
}; // end class TrivialSampler


template <class MeshType>
class MeshSampler
{
public:
  typedef typename MeshType::VertexType  VertexType;
  typedef typename MeshType::FaceType    FaceType;
  typedef typename MeshType::CoordType   CoordType;

  MeshSampler(MeshType &_m):m(_m){
    perFaceNormal = false;
  }
  MeshType &m;

  bool perFaceNormal;  // default false; if true the sample normal is the face normal, otherwise it is interpolated

  void reset()
  {
    m.Clear();
  }

  void AddVert(const VertexType &p)
  {
    tri::Allocator<MeshType>::AddVertices(m,1);
    m.vert.back().ImportData(p);
  }

  void AddFace(const FaceType &f, CoordType p)
  {
    tri::Allocator<MeshType>::AddVertices(m,1);
    m.vert.back().P() = f.cP(0)*p[0] + f.cP(1)*p[1] +f.cP(2)*p[2];
    if(perFaceNormal) m.vert.back().N() = f.cN();
       else m.vert.back().N() = f.cV(0)->N()*p[0] + f.cV(1)->N()*p[1] + f.cV(2)->N()*p[2];
    if(tri::HasPerVertexQuality(m) )
       m.vert.back().Q() = f.cV(0)->Q()*p[0] + f.cV(1)->Q()*p[1] + f.cV(2)->Q()*p[2];
  }
}; // end class BaseSampler



/* This sampler is used to perform compute the Hausdorff measuring.
 * It keep internally the spatial indexing structure used to find the closest point
 * and the partial integration results needed to compute the average and rms error values.
 * Averaged values assume that the samples are equi-distributed (e.g. a good unbiased montecarlo sampling of the surface).
 */
template <class MeshType>
class HausdorffSampler
{
  typedef typename MeshType::FaceType    FaceType;
  typedef typename MeshType::VertexType    VertexType;
  typedef typename MeshType::CoordType   CoordType;
  typedef typename MeshType::ScalarType   ScalarType;
  typedef GridStaticPtr<FaceType, ScalarType > MetroMeshFaceGrid;
  typedef GridStaticPtr<VertexType, ScalarType > MetroMeshVertexGrid;

public:

  HausdorffSampler(MeshType* _m, MeshType* _sampleMesh=0, MeshType* _closestMesh=0 ) :markerFunctor(_m)
  {
    m=_m;
    init(_sampleMesh,_closestMesh);
  }

  MeshType *m;           /// the mesh for which we search the closest points.
  MeshType *samplePtMesh;  /// the mesh containing the sample points
  MeshType *closestPtMesh; /// the mesh containing the corresponding closest points that have been found

  MetroMeshVertexGrid   unifGridVert;
  MetroMeshFaceGrid   unifGridFace;

  // Parameters
  double          min_dist;
  double          max_dist;
  double          mean_dist;
  double          RMS_dist;   /// from the wikipedia defintion RMS DIST is sqrt(Sum(distances^2)/n), here we store Sum(distances^2)
  double          volume;
  double          area_S1;
  Histogramf hist;
  // globals parameters driving the samples.
  int             n_total_samples;
  int             n_samples;
  bool useVertexSampling;
  ScalarType dist_upper_bound;  // samples that have a distance beyond this threshold distance are not considered.
  typedef typename tri::FaceTmark<MeshType> MarkerFace;
  MarkerFace markerFunctor;


  float getMeanDist() const { return mean_dist / n_total_samples; }
  float getMinDist() const { return min_dist ; }
  float getMaxDist() const { return max_dist ; }
  float getRMSDist() const { return sqrt(RMS_dist / n_total_samples); }

  void init(MeshType* _sampleMesh=0, MeshType* _closestMesh=0 )
  {
    samplePtMesh =_sampleMesh;
    closestPtMesh = _closestMesh;
    if(m)
    {
      tri::UpdateNormal<MeshType>::PerFaceNormalized(*m);
      if(m->fn==0) useVertexSampling = true;
      else useVertexSampling = false;

      if(useVertexSampling) unifGridVert.Set(m->vert.begin(),m->vert.end());
      else  unifGridFace.Set(m->face.begin(),m->face.end());
      markerFunctor.SetMesh(m);
      hist.SetRange(0.0, m->bbox.Diag()/100.0, 100);
    }
    min_dist = std::numeric_limits<double>::max();
    max_dist = 0;
    mean_dist =0;
    RMS_dist = 0;
    n_total_samples = 0;
  }

  void AddFace(const FaceType &f, CoordType interp)
  {
    CoordType startPt = f.cP(0)*interp[0] + f.cP(1)*interp[1] +f.cP(2)*interp[2]; // point to be sampled
    CoordType startN  = f.cV(0)->cN()*interp[0] + f.cV(1)->cN()*interp[1] +f.cV(2)->cN()*interp[2]; // Normal of the interpolated point
    AddSample(startPt,startN); // point to be sampled);
  }

  void AddVert(VertexType &p)
  {
    p.Q()=AddSample(p.cP(),p.cN());
  }


  float AddSample(const CoordType &startPt,const CoordType &startN)
  {
    // the results
    CoordType       closestPt;
    ScalarType dist = dist_upper_bound;

    // compute distance between startPt and the mesh S2
    FaceType   *nearestF=0;
    VertexType   *nearestV=0;
    vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
    dist=dist_upper_bound;
    if(useVertexSampling)
      nearestV =  tri::GetClosestVertex<MeshType,MetroMeshVertexGrid>(*m,unifGridVert,startPt,dist_upper_bound,dist);
    else
      nearestF =  unifGridFace.GetClosest(PDistFunct,markerFunctor,startPt,dist_upper_bound,dist,closestPt);

    // update distance measures
    if(dist == dist_upper_bound)
      return dist;

    if(dist > max_dist) max_dist = dist;        // L_inf
    if(dist < min_dist) min_dist = dist;        // L_inf

    mean_dist += dist;	        // L_1
    RMS_dist  += dist*dist;     // L_2
    n_total_samples++;

    hist.Add((float)fabs(dist));
    if(samplePtMesh)
    {
      tri::Allocator<MeshType>::AddVertices(*samplePtMesh,1);
      samplePtMesh->vert.back().P() = startPt;
      samplePtMesh->vert.back().Q() = dist;
      samplePtMesh->vert.back().N() = startN;
    }
    if(closestPtMesh)
    {
      tri::Allocator<MeshType>::AddVertices(*closestPtMesh,1);
      closestPtMesh->vert.back().P() = closestPt;
      closestPtMesh->vert.back().Q() = dist;
      closestPtMesh->vert.back().N() = startN;
    }
    return dist;
  }
}; // end class HausdorffSampler



/* This sampler is used to transfer the detail of a mesh onto another one.
 * It keep internally the spatial indexing structure used to find the closest point
 */
template <class MeshType>
class RedetailSampler
{
  typedef typename MeshType::FaceType    FaceType;
  typedef typename MeshType::VertexType    VertexType;
  typedef typename MeshType::CoordType   CoordType;
  typedef typename MeshType::ScalarType   ScalarType;  
  typedef GridStaticPtr<FaceType, ScalarType > MetroMeshGrid;
  typedef GridStaticPtr<VertexType, ScalarType > VertexMeshGrid;

public:

  RedetailSampler():m(0) {}

  MeshType *m;           /// the source mesh for which we search the closest points (e.g. the mesh from which we take colors etc).
  CallBackPos *cb;
  int sampleNum;  // the expected number of samples. Used only for the callback
  int sampleCnt;
  MetroMeshGrid   unifGridFace;
  VertexMeshGrid   unifGridVert;
  bool useVertexSampling;

  // Parameters
  typedef tri::FaceTmark<MeshType> MarkerFace;
  MarkerFace markerFunctor;

  bool coordFlag;
  bool colorFlag;
  bool normalFlag;
  bool qualityFlag;
  bool selectionFlag;
  bool storeDistanceAsQualityFlag;
  float dist_upper_bound;
  void init(MeshType *_m, CallBackPos *_cb=0, int targetSz=0)
  {
    coordFlag=false;
    colorFlag=false;
    qualityFlag=false;
    selectionFlag=false;
    storeDistanceAsQualityFlag=false;
    m=_m;
    tri::UpdateNormal<MeshType>::PerFaceNormalized(*m);
    if(m->fn==0) useVertexSampling = true;
    else useVertexSampling = false;

    if(useVertexSampling) unifGridVert.Set(m->vert.begin(),m->vert.end());
    else  unifGridFace.Set(m->face.begin(),m->face.end());
    markerFunctor.SetMesh(m);
    // sampleNum and sampleCnt are used only for the progress callback.
    cb=_cb;
    sampleNum = targetSz;
    sampleCnt = 0;
  }

  // this function is called for each vertex of the target mesh.
  // and retrieve the closest point on the source mesh.
  void AddVert(VertexType &p)
  {
    assert(m);
    // the results
    CoordType       closestPt,      normf, bestq, ip;
    ScalarType dist = dist_upper_bound;
    const CoordType &startPt= p.cP();
    // compute distance between startPt and the mesh S2
    if(useVertexSampling)
    {
      VertexType   *nearestV=0;
      nearestV =  tri::GetClosestVertex<MeshType,VertexMeshGrid>(*m,unifGridVert,startPt,dist_upper_bound,dist); //(PDistFunct,markerFunctor,startPt,dist_upper_bound,dist,closestPt);
      if(cb) cb(sampleCnt++*100/sampleNum,"Resampling Vertex attributes");
      if(storeDistanceAsQualityFlag)  p.Q() = dist;
      if(dist == dist_upper_bound) return ;

      if(coordFlag) p.P()=nearestV->P();
      if(colorFlag) p.C() = nearestV->C();
      if(normalFlag) p.N() = nearestV->N();
      if(qualityFlag) p.Q()= nearestV->Q();
      if(selectionFlag) if(nearestV->IsS()) p.SetS();
    }
    else
    {
      FaceType   *nearestF=0;
      vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
      dist=dist_upper_bound;
      if(cb) cb(sampleCnt++*100/sampleNum,"Resampling Vertex attributes");
      nearestF =  unifGridFace.GetClosest(PDistFunct,markerFunctor,startPt,dist_upper_bound,dist,closestPt);
      if(dist == dist_upper_bound) return ;

      CoordType interp;
      InterpolationParameters(*nearestF,(*nearestF).cN(),closestPt, interp);
      interp[2]=1.0-interp[1]-interp[0];

      if(coordFlag) p.P()=closestPt;
      if(colorFlag) p.C().lerp(nearestF->V(0)->C(),nearestF->V(1)->C(),nearestF->V(2)->C(),interp);
      if(normalFlag) p.N() = nearestF->V(0)->N()*interp[0] + nearestF->V(1)->N()*interp[1] + nearestF->V(2)->N()*interp[2];
      if(qualityFlag) p.Q()= nearestF->V(0)->Q()*interp[0] + nearestF->V(1)->Q()*interp[1] + nearestF->V(2)->Q()*interp[2];
      if(selectionFlag) if(nearestF->IsS()) p.SetS();
    }
  }
}; // end class RedetailSampler




/**
 \brief Main Class of the Sampling framework.

This class allows you to perform various kind of random/procedural point sampling over a triangulated surface.
The class is templated over the PointSampler object that allows to customize the use of the generated samples.


**/


template <class MeshType, class VertexSampler = TrivialSampler< MeshType> >
class SurfaceSampling
{
  typedef typename MeshType::CoordType       CoordType;
  typedef typename MeshType::BoxType         BoxType;
  typedef typename MeshType::ScalarType      ScalarType;
  typedef typename MeshType::VertexType      VertexType;
  typedef typename MeshType::VertexPointer   VertexPointer;
  typedef typename MeshType::VertexIterator  VertexIterator;
  typedef typename MeshType::EdgeType        EdgeType;
  typedef typename MeshType::EdgeIterator    EdgeIterator;
  typedef typename MeshType::FaceType        FaceType;
  typedef typename MeshType::FacePointer     FacePointer;
  typedef typename MeshType::FaceIterator    FaceIterator;
  typedef typename MeshType::FaceContainer   FaceContainer;

  typedef typename vcg::SpatialHashTable<FaceType, ScalarType> MeshSHT;
  typedef typename vcg::SpatialHashTable<FaceType, ScalarType>::CellIterator MeshSHTIterator;
  typedef typename vcg::SpatialHashTable<VertexType, ScalarType> MontecarloSHT;
  typedef typename vcg::SpatialHashTable<VertexType, ScalarType>::CellIterator MontecarloSHTIterator;
  typedef typename vcg::SpatialHashTable<VertexType, ScalarType> SampleSHT;
  typedef typename vcg::SpatialHashTable<VertexType, ScalarType>::CellIterator SampleSHTIterator;

  typedef typename MeshType::template PerVertexAttributeHandle<float> PerVertexFloatAttribute;

public:

static math::MarsenneTwisterRNG &SamplingRandomGenerator()
{
    static math::MarsenneTwisterRNG rnd;
    return rnd;
}

// Returns an integer random number in the [0,i-1] interval using the improve Marsenne-Twister method.
// this functor is needed for passing it to the std functions.
static unsigned int RandomInt(unsigned int i)
{
    return (SamplingRandomGenerator().generate(i));
}

// Returns a random number in the [0,1) real interval using the improved Marsenne-Twister method.
static double RandomDouble01()
{
    return SamplingRandomGenerator().generate01();
}

#define FAK_LEN 1024
static double LnFac(int n) {
   // Tabled log factorial function. gives natural logarithm of n!

   // define constants
   static const double                 // coefficients in Stirling approximation
      C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
      C1 =  1./12.,
      C3 = -1./360.;
   // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
   // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
   // static variables
   static double fac_table[FAK_LEN];   // table of ln(n!):
   static bool initialized = false;         // remember if fac_table has been initialized


   if (n < FAK_LEN) {
      if (n <= 1) {
         if (n < 0) assert(0);//("Parameter negative in LnFac function");
         return 0;
      }
      if (!initialized) {              // first time. Must initialize table
         // make table of ln(n!)
         double sum = fac_table[0] = 0.;
         for (int i=1; i<FAK_LEN; i++) {
            sum += log(double(i));
            fac_table[i] = sum;
         }
         initialized = true;
      }
      return fac_table[n];
   }
   // not found in table. use Stirling approximation
   double  n1, r;
   n1 = n;  r  = 1. / n1;
   return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}

static int  PoissonRatioUniforms(double L) {
   /*

   This subfunction generates a integer with the poisson
   distribution using the ratio-of-uniforms rejection method (PRUAt).
   This approach is STABLE even for large L (e.g. it does not suffer from the overflow limit of the classical Knuth implementation)
   Execution time does not depend on L, except that it matters whether
   is within the range where ln(n!) is tabulated.

   Reference:

   E. Stadlober
   "The ratio of uniforms approach for generating discrete random variates".
   Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.

   Partially adapted/inspired from some subfunctions of the Agner Fog stocc library ( www.agner.org/random )
   Same licensing scheme.

   */
  // constants

  const double SHAT1 = 2.943035529371538573;    // 8/e
  const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)
  double u;                                          // uniform random
  double lf;                                         // ln(f(x))
  double x;                                          // real sample
  int k;                                         // integer sample

  double   pois_a = L + 0.5;                               // hat center
  int mode = (int)L;                      // mode
  double   pois_g  = log(L);
  double    pois_f0 = mode * pois_g - LnFac(mode);          // value at mode
  double   pois_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;         // hat width
  double   pois_bound = (int)(pois_a + 6.0 * pois_h);  // safety-bound

  while(1) {
      u = RandomDouble01();
      if (u == 0) continue;                           // avoid division by 0
      x = pois_a + pois_h * (RandomDouble01() - 0.5) / u;
      if (x < 0 || x >= pois_bound) continue;         // reject if outside valid range
      k = (int)(x);
      lf = k * pois_g - LnFac(k) - pois_f0;
      if (lf >= u * (4.0 - u) - 3.0) break;           // quick acceptance
      if (u * (u - lf) > 1.0) continue;               // quick rejection
      if (2.0 * log(u) <= lf) break;                  // final acceptance
   }
   return k;
}


/**
  algorithm poisson random number (Knuth):
    init:
         Let L ← e^−λ, k ← 0 and p ← 1.
    do:
         k ← k + 1.
         Generate uniform random number u in [0,1] and let p ← p × u.
    while p > L.
    return k − 1.

  */
static int Poisson(double lambda)
{
  if(lambda>50) return PoissonRatioUniforms(lambda);
  double L = exp(-lambda);
  int k =0;
  double p = 1.0;
  do
  {
    k = k+1;
    p = p*RandomDouble01();
  } while (p>L);

  return k -1;
}


static void AllVertex(MeshType & m, VertexSampler &ps)
{
	AllVertex(m, ps, false);
}

static void AllVertex(MeshType & m, VertexSampler &ps, bool onlySelected)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if(!(*vi).IsD())
			if ((!onlySelected) || ((*vi).IsS()))
			{
				ps.AddVert(*vi);
			}
}

/// Sample the vertices in a weighted way. Each vertex has a probability of being chosen
/// that is proportional to its quality.
/// It assumes that you are asking a number of vertices smaller than nv;
/// Algorithm:
/// 1) normalize quality so that sum q == 1;
/// 2) shuffle vertices.
/// 3) for each vertices choose it if rand > thr;

static void VertexWeighted(MeshType & m, VertexSampler &ps, int sampleNum)
{
    ScalarType qSum = 0;
    VertexIterator vi;
    for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
            if(!(*vi).IsD())
                        qSum += (*vi).Q();

    ScalarType samplePerUnit = sampleNum/qSum;
    ScalarType floatSampleNum =0;
    std::vector<VertexPointer> vertVec;
    FillAndShuffleVertexPointerVector(m,vertVec);

    std::vector<bool> vertUsed(m.vn,false);

    int i=0; int cnt=0;
    while(cnt < sampleNum)
        {
            if(vertUsed[i])
                {
                        floatSampleNum += vertVec[i]->Q() * samplePerUnit;
                        int vertSampleNum   = (int) floatSampleNum;
                        floatSampleNum -= (float) vertSampleNum;

                        // for every sample p_i in T...
                        if(vertSampleNum > 1)
                            {
                                ps.AddVert(*vertVec[i]);
                                cnt++;
                                vertUsed[i]=true;
                            }
                }
            i = (i+1)%m.vn;
        }
}

/// Sample the vertices in a uniform way. Each vertex has a probability of being chosen
/// that is proportional to the area it represent.
static void VertexAreaUniform(MeshType & m, VertexSampler &ps, int sampleNum)
{
    VertexIterator vi;
    for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        if(!(*vi).IsD())
                        (*vi).Q() = 0;

    FaceIterator fi;
    for(fi = m.face.begin(); fi != m.face.end(); ++fi)
        if(!(*fi).IsD())
        {
            ScalarType areaThird = DoubleArea(*fi)/6.0;
            (*fi).V(0)->Q()+=areaThird;
            (*fi).V(1)->Q()+=areaThird;
            (*fi).V(2)->Q()+=areaThird;
        }

    VertexWeighted(m,ps,sampleNum);
}

static void	FillAndShuffleFacePointerVector(MeshType & m, std::vector<FacePointer> &faceVec)
{
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
        if(!(*fi).IsD())	faceVec.push_back(&*fi);

    assert((int)faceVec.size()==m.fn);

    unsigned int (*p_myrandom)(unsigned int) = RandomInt;
    std::random_shuffle(faceVec.begin(),faceVec.end(), p_myrandom);
}
static void	FillAndShuffleVertexPointerVector(MeshType & m, std::vector<VertexPointer> &vertVec)
{
    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
                if(!(*vi).IsD())	vertVec.push_back(&*vi);

    assert((int)vertVec.size()==m.vn);

    unsigned int (*p_myrandom)(unsigned int) = RandomInt;
    std::random_shuffle(vertVec.begin(),vertVec.end(), p_myrandom);
}

/// Sample the vertices in a uniform way. Each vertex has the same probabiltiy of being chosen.
static void VertexUniform(MeshType & m, VertexSampler &ps, int sampleNum, bool onlySelected)
{
	if (sampleNum >= m.vn) {
		AllVertex(m, ps, onlySelected);
		return;
	}

	std::vector<VertexPointer> vertVec;
	FillAndShuffleVertexPointerVector(m, vertVec);

	int added = 0;
	for (int i = 0; ((i < m.vn) && (added < sampleNum)); ++i)
		if (!(*vertVec[i]).IsD())
			if ((!onlySelected) || (*vertVec[i]).IsS())
			{
				ps.AddVert(*vertVec[i]);
				added++;
			}
		
}


static void VertexUniform(MeshType & m, VertexSampler &ps, int sampleNum)
{
	VertexUniform(m, ps, sampleNum, false);
}


/// Perform an uniform sampling over an EdgeMesh.
///
/// It assumes that the mesh is 1-manifold.
/// each connected component is sampled in a independent way.
/// For each component of lenght <L> we place on it floor(L/radius)+1 samples.
/// (if conservative argument is false we place ceil(L/radius)+1 samples)
///
static void EdgeMeshUniform(MeshType &m, VertexSampler &ps, float radius, bool conservative = true)
{
	tri::RequireEEAdjacency(m);
	tri::RequireCompactness(m);
	tri::RequirePerEdgeFlags(m);
	tri::RequirePerVertexFlags(m);
	tri::UpdateTopology<MeshType>::EdgeEdge(m);
	tri::UpdateFlags<MeshType>::EdgeClearV(m);

	for (EdgeIterator ei = m.edge.begin(); ei != m.edge.end(); ++ei)
	{
		if (!ei->IsV())
		{
			edge::Pos<EdgeType> ep(&*ei,0);
			edge::Pos<EdgeType> startep     = ep;
			VertexPointer       startVertex = 0;
			do // first loop to search a boundary component.
			{
				ep.NextE();
				if (ep.IsBorder())
					break;
			} while (startep != ep);
			if (!ep.IsBorder())
			{
				// it's a loop
				startVertex = ep.V();
			}
			else
			{
				// to keep the uniform resampling order-independent
				// start from the border with 'lowest' point
				edge::Pos<EdgeType> altEp = ep;
				do {
					altEp.NextE();
				} while (!altEp.IsBorder());

				if (altEp.V()->cP() < ep.V()->cP())
				{
					ep = altEp;
				}
			}

			ScalarType totalLen = 0;
			ep.FlipV();
			// second loop to compute length of the chain.
			do
			{
				ep.E()->SetV();
				totalLen += Distance(ep.V()->cP(), ep.VFlip()->cP());
				ep.NextE();
			} while (!ep.IsBorder() && ep.V() != startVertex);
			ep.E()->SetV();
			totalLen += Distance(ep.V()->cP(), ep.VFlip()->cP());

			// Third loop actually perform the sampling.
			int sampleNum = conservative ? floor(totalLen / radius) : ceil(totalLen / radius);

			ScalarType sampleDist = totalLen / sampleNum;
//			printf("Found a chain of %f with %i samples every %f (%f)\n", totalLen, sampleNum, sampleDist, radius);

			ScalarType curLen = 0;
			int sampleCnt     = 1;
			ps.AddEdge(*(ep.E()), ep.VInd() == 0 ? 0.0 : 1.0);
			do
			{
				ep.NextE();
				assert(ep.E()->IsV());
				ScalarType edgeLen = Distance(ep.V()->cP(), ep.VFlip()->cP());
				ScalarType d0      = curLen;
				ScalarType d1      = d0 + edgeLen;

				while (d1 > sampleCnt * sampleDist && sampleCnt < sampleNum)
				{
					ScalarType off = (sampleCnt * sampleDist - d0) / edgeLen;
//					printf("edgeLen %f off %f samplecnt %i\n", edgeLen, off, sampleCnt);
					ps.AddEdge(*(ep.E()), ep.VInd() == 0 ? 1.0 - off : off);
					sampleCnt++;
				}
				curLen += edgeLen;
			} while(!ep.IsBorder() && ep.V() != startVertex);

			if(ep.V() != startVertex)
				ps.AddEdge(*(ep.E()), ep.VInd() == 0 ? 0.0 : 1.0);
		}
	}
}


/// \brief Sample all the border corner vertices
///
/// It assumes that the border flag have been set over the mesh both for vertex and for faces.
/// All the vertices on the border where the edges of the boundary of the surface forms an angle smaller than the given threshold are sampled.
///
static void VertexBorderCorner(MeshType & m, VertexSampler &ps, float angleRad)
{
  typename MeshType::template PerVertexAttributeHandle<float> angleSumH = tri::Allocator<MeshType>:: template GetPerVertexAttribute<float> (m);

  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
    angleSumH[vi]=0;

  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    for(int i=0;i<3;++i)
    {
      angleSumH[fi->V(i)] += vcg::Angle(fi->P2(i) - fi->P0(i),fi->P1(i) - fi->P0(i));
    }
  }

  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    if((angleSumH[vi]<angleRad && vi->IsB())||
       (angleSumH[vi]>(360-angleRad) && vi->IsB()))
        ps.AddVert(*vi);
  }

  tri::Allocator<MeshType>:: template DeletePerVertexAttribute<float> (m,angleSumH);
}

/// \brief Sample all the border vertices
///
/// It assumes that the border flag have been set over the mesh.
/// All the vertices on the border are sampled.
///
static void VertexBorder(MeshType & m, VertexSampler &ps)
{
  VertexBorderCorner(m,ps,std::numeric_limits<ScalarType>::max());
}

/// Sample all the crease vertices.
/// It assumes that the crease edges had been marked as non-faux edges
/// for example by using
/// tri::UpdateFlags<MeshType>::FaceFauxCrease(mesh,creaseAngleRad);
/// Then it chooses all the vertices where there are at least three non faux edges.
///
static void VertexCrease(MeshType & m, VertexSampler &ps)
{
  typedef typename UpdateTopology<MeshType>::PEdge SimpleEdge;
  std::vector< SimpleEdge > Edges;
  typename std::vector< SimpleEdge >::iterator ei;
  UpdateTopology<MeshType>::FillUniqueEdgeVector(m,Edges,false);

  typename MeshType::template PerVertexAttributeHandle  <int> hv = tri::Allocator<MeshType>:: template GetPerVertexAttribute<int> (m);

  for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
  {
    hv[ei->v[0]]++;
    hv[ei->v[1]]++;
  }

  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    if(hv[vi]>2)
        ps.AddVert(*vi);
  }
}


static void FaceUniform(MeshType & m, VertexSampler &ps, int sampleNum)
{
    if(sampleNum>=m.fn) {
      AllFace(m,ps);
        return;
    }

    std::vector<FacePointer> faceVec;
    FillAndShuffleFacePointerVector(m,faceVec);

    for(int i =0; i< sampleNum; ++i)
        ps.AddFace(*faceVec[i],Barycenter(*faceVec[i]));
}

static void AllFace(MeshType & m, VertexSampler &ps)
{
    FaceIterator fi;
    for(fi=m.face.begin();fi!=m.face.end();++fi)
                if(!(*fi).IsD())
                {
                    ps.AddFace(*fi,Barycenter(*fi));
                }
}


static void AllEdge(MeshType & m, VertexSampler &ps)
{
  // Edge sampling.
  typedef typename UpdateTopology<MeshType>::PEdge SimpleEdge;
  std::vector< SimpleEdge > Edges;
  typename std::vector< SimpleEdge >::iterator ei;
  UpdateTopology<MeshType>::FillUniqueEdgeVector(m,Edges);

  for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    ps.AddFace(*(*ei).f,ei->EdgeBarycentricToFaceBarycentric(0.5));
}

// Regular Uniform Edge sampling
// Each edge is subdivided in a number of pieces proprtional to its length
// Sample are choosen without touching the vertices.

static void EdgeUniform(MeshType & m, VertexSampler &ps,int sampleNum, bool sampleFauxEdge=true)
{
  typedef typename UpdateTopology<MeshType>::PEdge SimpleEdge;

  std::vector< SimpleEdge > Edges;
  UpdateTopology<MeshType>::FillUniqueEdgeVector(m,Edges,sampleFauxEdge);
  // First loop compute total edge length;
  float edgeSum=0;
  typename std::vector< SimpleEdge >::iterator ei;
  for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    edgeSum+=Distance((*ei).v[0]->P(),(*ei).v[1]->P());

  float sampleLen = edgeSum/sampleNum;
  float rest=0;
  for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
  {
    float len = Distance((*ei).v[0]->P(),(*ei).v[1]->P());
    float samplePerEdge = floor((len+rest)/sampleLen);
    rest = (len+rest) - samplePerEdge * sampleLen;
    float step = 1.0/(samplePerEdge+1);
    for(int i=0;i<samplePerEdge;++i)
    {
      CoordType interp(0,0,0);
      interp[ (*ei).z     ]=step*(i+1);
      interp[((*ei).z+1)%3]=1.0-step*(i+1);
      ps.AddFace(*(*ei).f,interp);
    }
  }
}

// Generate the barycentric coords of a random point over a single face,
// with a uniform distribution over the triangle.
// It uses the parallelogram folding trick.
static CoordType RandomBarycentric()
{
  return math::GenerateBarycentricUniform<ScalarType>(SamplingRandomGenerator());
}

// Given a triangle return a random point over it
static CoordType RandomPointInTriangle(const FaceType &f)
{
    CoordType u = RandomBarycentric();
    return f.cP(0)*u[0] + f.cP(1)*u[1] + f.cP(2)*u[2];
}

static void StratifiedMontecarlo(MeshType & m, VertexSampler &ps,int sampleNum)
{
    ScalarType area = Stat<MeshType>::ComputeMeshArea(m);
    ScalarType samplePerAreaUnit = sampleNum/area;
    // Montecarlo sampling.
    double  floatSampleNum = 0.0;

    FaceIterator fi;
    for(fi=m.face.begin(); fi != m.face.end(); fi++)
        if(!(*fi).IsD())
        {
            // compute # samples in the current face (taking into account of the remainders)
            floatSampleNum += 0.5*DoubleArea(*fi) * samplePerAreaUnit;
            int faceSampleNum   = (int) floatSampleNum;

            // for every sample p_i in T...
            for(int i=0; i < faceSampleNum; i++)
                    ps.AddFace(*fi,RandomBarycentric());
            floatSampleNum -= (double) faceSampleNum;
        }
}

/**
  This function compute montecarlo distribution with an approximate number of
  samples exploiting the poisson distribution approximation of the binomial distribution.

  For a given triangle t of area a_t, in a Mesh of area A,
  if we take n_s sample over the mesh, the number of samples that falls in t
  follows the poisson distribution of P(lambda ) with lambda = n_s * (a_t/A).

  To approximate the Binomial we use a Poisson distribution with parameter
  \lambda = np can be used as an approximation to B(n,p)
  (it works if n is sufficiently large and p is sufficiently small).

  */

static void MontecarloPoisson(MeshType & m, VertexSampler &ps,int sampleNum)
{
  ScalarType area = Stat<MeshType>::ComputeMeshArea(m);
  ScalarType samplePerAreaUnit = sampleNum/area;

  FaceIterator fi;
  for(fi=m.face.begin(); fi != m.face.end(); fi++)
    if(!(*fi).IsD())
    {
      float areaT=DoubleArea(*fi) * 0.5f;
      int faceSampleNum = Poisson(areaT*samplePerAreaUnit);

      // for every sample p_i in T...
      for(int i=0; i < faceSampleNum; i++)
          ps.AddFace(*fi,RandomBarycentric());
//      SampleNum -= (double) faceSampleNum;
    }
}


/**
  This function computes a montecarlo distribution with an EXACT number of samples.
  it works by generating a sequence of consecutive segments proportional to the triangle areas
  and actually shooting sample over this line
  */

static void EdgeMontecarlo(MeshType & m, VertexSampler &ps, int sampleNum, bool sampleAllEdges)
{
  typedef typename UpdateTopology<MeshType>::PEdge SimpleEdge;
  std::vector< SimpleEdge > Edges;
  UpdateTopology<MeshType>::FillUniqueEdgeVector(m,Edges,sampleAllEdges);

  assert(!Edges.empty());

  typedef  std::pair<ScalarType, SimpleEdge*> IntervalType;
  std::vector< IntervalType > intervals (Edges.size()+1);
  int i=0;
  intervals[i]=std::make_pair(0,(SimpleEdge*)(0));
  // First loop: build a sequence of consecutive segments proportional to the edge lenghts.
  typename std::vector< SimpleEdge >::iterator ei;
  for(ei=Edges.begin(); ei != Edges.end(); ei++)
  {
    intervals[i+1]=std::make_pair(intervals[i].first+Distance((*ei).v[0]->P(),(*ei).v[1]->P()), &*ei);
    ++i;
  }

  // Second Loop get a point on the line 0...Sum(edgeLen) to pick a point;
  ScalarType edgeSum = intervals.back().first;
  for(i=0;i<sampleNum;++i)
  {
    ScalarType val = edgeSum * RandomDouble01();
    // lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value.
    // E.g. An iterator pointing to the first element "not less than" val, or end() if every element is less than val.
    typename std::vector<IntervalType>::iterator it = lower_bound(intervals.begin(),intervals.end(),std::make_pair(val,(SimpleEdge*)(0)) );
    assert(it != intervals.end() && it != intervals.begin());
    assert( ( (*(it-1)).first < val ) && ((*(it)).first >= val) );
    SimpleEdge * ep=(*it).second;
    ps.AddFace( *(ep->f), ep->EdgeBarycentricToFaceBarycentric(RandomDouble01()) );
  }
}

/**
  This function computes a montecarlo distribution with an EXACT number of samples.
  it works by generating a sequence of consecutive segments proportional to the triangle areas
  and actually shooting sample over this line
  */

static void Montecarlo(MeshType & m, VertexSampler &ps,int sampleNum)
{
    typedef  std::pair<ScalarType, FacePointer> IntervalType;
    std::vector< IntervalType > intervals (m.fn+1);
    FaceIterator fi;
    int i=0;
    intervals[i]=std::make_pair(0,FacePointer(0));
    // First loop: build a sequence of consecutive segments proportional to the triangle areas.
    for(fi=m.face.begin(); fi != m.face.end(); fi++)
        if(!(*fi).IsD())
        {
            intervals[i+1]=std::make_pair(intervals[i].first+0.5*DoubleArea(*fi), &*fi);
            ++i;
        }
    ScalarType meshArea = intervals.back().first;
    for(i=0;i<sampleNum;++i)
        {
            ScalarType val = meshArea * RandomDouble01();
            // lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value.
            // E.g. An iterator pointing to the first element "not less than" val, or end() if every element is less than val.
            typename std::vector<IntervalType>::iterator it = lower_bound(intervals.begin(),intervals.end(),std::make_pair(val,FacePointer(0)) );
            assert(it != intervals.end());
            assert(it != intervals.begin());
            assert( (*(it-1)).first <val );
            assert( (*(it)).first >= val);
            ps.AddFace( *(*it).second, RandomBarycentric() );
        }
}

static ScalarType WeightedArea(FaceType &f, PerVertexFloatAttribute &wH)
{
    ScalarType averageQ = ( wH[f.V(0)] + wH[f.V(1)] + wH[f.V(2)] )/3.0;
    return averageQ*averageQ*DoubleArea(f)/2.0;
}

/// Compute a sampling of the surface that is weighted by the quality and a variance
///
/// We use the quality as linear distortion of density.
/// We consider each triangle as scaled between 1 and 1/variance linearly according quality.
///
/// In practice with variance 2 the average distance between sample will double where the quality is maxima.
/// If you have two same area region A with q==-1 and B with q==1, if variance==2 the A will have 4 times more samples than B
///
static void WeightedMontecarlo(MeshType & m, VertexSampler &ps,int sampleNum, float variance)
{
  tri::RequirePerVertexQuality(m);
  tri::RequireCompactness(m);
  PerVertexFloatAttribute rH = tri::Allocator<MeshType>:: template GetPerVertexAttribute<float> (m,"radius");
  InitRadiusHandleFromQuality(m, rH, 1.0, variance, true);

  ScalarType weightedArea = 0;
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      weightedArea += WeightedArea(*fi,rH);

  ScalarType samplePerAreaUnit = sampleNum/weightedArea;
  // Montecarlo sampling.
  double  floatSampleNum = 0.0;
  for(FaceIterator fi=m.face.begin(); fi != m.face.end(); fi++)
    {
      // compute # samples in the current face (taking into account of the remainders)
      floatSampleNum += WeightedArea(*fi,rH) * samplePerAreaUnit;
      int faceSampleNum   = (int) floatSampleNum;

      // for every sample p_i in T...
      for(int i=0; i < faceSampleNum; i++)
        ps.AddFace(*fi,RandomBarycentric());

      floatSampleNum -= (double) faceSampleNum;
    }
}


// Subdivision sampling of a single face.
// return number of added samples

static int SingleFaceSubdivision(int sampleNum, const CoordType & v0, const CoordType & v1, const CoordType & v2, VertexSampler &ps, FacePointer fp, bool randSample)
{
    // recursive face subdivision.
    if(sampleNum == 1)
    {
        // ground case.
        CoordType SamplePoint;
        if(randSample)
        {
            CoordType rb=RandomBarycentric();
            SamplePoint=v0*rb[0]+v1*rb[1]+v2*rb[2];
        }
        else SamplePoint=((v0+v1+v2)*(1.0f/3.0f));

        ps.AddFace(*fp,SamplePoint);
        return 1;
    }

    int s0 = sampleNum /2;
    int s1 = sampleNum-s0;
    assert(s0>0);
    assert(s1>0);

    ScalarType w0 = ScalarType(s1)/ScalarType(sampleNum);
    ScalarType w1 = 1.0-w0;
    // compute the longest edge.
    ScalarType  maxd01 = SquaredDistance(v0,v1);
    ScalarType  maxd12 = SquaredDistance(v1,v2);
    ScalarType  maxd20 = SquaredDistance(v2,v0);
    int     res;
    if(maxd01 > maxd12)
        if(maxd01 > maxd20)     res = 0;
    else                    res = 2;
    else
        if(maxd12 > maxd20)     res = 1;
    else                    res = 2;

    int faceSampleNum=0;
    // break the input triangle along the midpoint of the longest edge.
    CoordType  pp;
    switch(res)
    {
    case 0 :    pp = v0*w0 + v1*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,pp,v2,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
        break;
    case 1 :    pp =  v1*w0 + v2*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,v0,pp,v2,ps,fp,randSample);
        break;
    case 2 :    pp = v0*w0 + v2*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
        break;
    }
    return faceSampleNum;
}


/// Compute a sampling of the surface where the points are regularly scattered over the face surface using a recursive longest-edge subdivision rule.
static void FaceSubdivision(MeshType & m, VertexSampler &ps,int sampleNum, bool randSample)
{

    ScalarType area = Stat<MeshType>::ComputeMeshArea(m);
    ScalarType samplePerAreaUnit = sampleNum/area;
    std::vector<FacePointer> faceVec;
    FillAndShuffleFacePointerVector(m,faceVec);
    vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
    double  floatSampleNum = 0.0;
    int faceSampleNum;
    // Subdivision sampling.
    typename std::vector<FacePointer>::iterator fi;
    for(fi=faceVec.begin(); fi!=faceVec.end(); fi++)
    {
        const CoordType b0(1.0, 0.0, 0.0);
        const CoordType b1(0.0, 1.0, 0.0);
        const CoordType b2(0.0, 0.0, 1.0);
        // compute # samples in the current face.
        floatSampleNum += 0.5*DoubleArea(**fi) * samplePerAreaUnit;
        faceSampleNum          = (int) floatSampleNum;
        if(faceSampleNum>0)
            faceSampleNum = SingleFaceSubdivision(faceSampleNum,b0,b1,b2,ps,*fi,randSample);
        floatSampleNum -= (double) faceSampleNum;
    }
}
//---------
// Subdivision sampling of a single face.
// return number of added samples

static int SingleFaceSubdivisionOld(int sampleNum, const CoordType & v0, const CoordType & v1, const CoordType & v2, VertexSampler &ps, FacePointer fp, bool randSample)
{
    // recursive face subdivision.
    if(sampleNum == 1)
    {
        // ground case.
        CoordType SamplePoint;
        if(randSample)
        {
            CoordType rb=RandomBarycentric();
            SamplePoint=v0*rb[0]+v1*rb[1]+v2*rb[2];
        }
        else SamplePoint=((v0+v1+v2)*(1.0f/3.0f));

        CoordType SampleBary;
        InterpolationParameters(*fp,SamplePoint,SampleBary);
        ps.AddFace(*fp,SampleBary);
        return 1;
    }

    int s0 = sampleNum /2;
    int s1 = sampleNum-s0;
    assert(s0>0);
    assert(s1>0);

    ScalarType w0 = ScalarType(s1)/ScalarType(sampleNum);
    ScalarType w1 = 1.0-w0;
    // compute the longest edge.
    ScalarType  maxd01 = SquaredDistance(v0,v1);
    ScalarType  maxd12 = SquaredDistance(v1,v2);
    ScalarType  maxd20 = SquaredDistance(v2,v0);
    int     res;
    if(maxd01 > maxd12)
        if(maxd01 > maxd20)     res = 0;
    else                    res = 2;
    else
        if(maxd12 > maxd20)     res = 1;
    else                    res = 2;

    int faceSampleNum=0;
    // break the input triangle along the midpoint of the longest edge.
    CoordType  pp;
    switch(res)
    {
    case 0 :    pp = v0*w0 + v1*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,pp,v2,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
        break;
    case 1 :    pp =  v1*w0 + v2*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,v0,pp,v2,ps,fp,randSample);
        break;
    case 2 :    pp = v0*w0 + v2*w1;
        faceSampleNum+=SingleFaceSubdivision(s0,v0,v1,pp,ps,fp,randSample);
        faceSampleNum+=SingleFaceSubdivision(s1,pp,v1,v2,ps,fp,randSample);
        break;
    }
    return faceSampleNum;
}


/// Compute a sampling of the surface where the points are regularly scattered over the face surface using a recursive longest-edge subdivision rule.
static void FaceSubdivisionOld(MeshType & m, VertexSampler &ps,int sampleNum, bool randSample)
{

    ScalarType area = Stat<MeshType>::ComputeMeshArea(m);
    ScalarType samplePerAreaUnit = sampleNum/area;
    std::vector<FacePointer> faceVec;
    FillAndShuffleFacePointerVector(m,faceVec);
    tri::UpdateNormal<MeshType>::PerFaceNormalized(m);
    double  floatSampleNum = 0.0;
    int faceSampleNum;
    // Subdivision sampling.
    typename std::vector<FacePointer>::iterator fi;
    for(fi=faceVec.begin(); fi!=faceVec.end(); fi++)
    {
        // compute # samples in the current face.
        floatSampleNum += 0.5*DoubleArea(**fi) * samplePerAreaUnit;
        faceSampleNum          = (int) floatSampleNum;
        if(faceSampleNum>0)
            faceSampleNum = SingleFaceSubdivision(faceSampleNum,(**fi).V(0)->cP(), (**fi).V(1)->cP(), (**fi).V(2)->cP(),ps,*fi,randSample);
        floatSampleNum -= (double) faceSampleNum;
    }
}


//---------

// Similar Triangles sampling.
// Skip vertex and edges
// Sample per edges includes vertexes, so here we should expect  n_samples_per_edge >=4

static int SingleFaceSimilar(FacePointer fp, VertexSampler &ps, int n_samples_per_edge)
{
        int n_samples=0;
    int         i, j;
    float segmentNum=n_samples_per_edge -1 ;
        float segmentLen = 1.0/segmentNum;
        // face sampling.
    for(i=1; i < n_samples_per_edge-1; i++)
        for(j=1; j < n_samples_per_edge-1-i; j++)
        {
            //AddSample( v0 + (V1*(double)i + V2*(double)j) );
                        CoordType sampleBary(i*segmentLen,j*segmentLen, 1.0 - (i*segmentLen+j*segmentLen) ) ;
            n_samples++;
                        ps.AddFace(*fp,sampleBary);
        }
    return n_samples;
}
static int SingleFaceSimilarDual(FacePointer fp, VertexSampler &ps, int n_samples_per_edge, bool randomFlag)
{
        int n_samples=0;
    float         i, j;
    float segmentNum=n_samples_per_edge -1 ;
        float segmentLen = 1.0/segmentNum;
        // face sampling.
    for(i=0; i < n_samples_per_edge-1; i++)
        for(j=0; j < n_samples_per_edge-1-i; j++)
        {
            //AddSample( v0 + (V1*(double)i + V2*(double)j) );
                        CoordType V0((i+0)*segmentLen,(j+0)*segmentLen, 1.0 - ((i+0)*segmentLen+(j+0)*segmentLen) ) ;
                        CoordType V1((i+1)*segmentLen,(j+0)*segmentLen, 1.0 - ((i+1)*segmentLen+(j+0)*segmentLen) ) ;
                        CoordType V2((i+0)*segmentLen,(j+1)*segmentLen, 1.0 - ((i+0)*segmentLen+(j+1)*segmentLen) ) ;
                        n_samples++;
                        if(randomFlag) 	{
                                            CoordType rb=RandomBarycentric();
                                            ps.AddFace(*fp, V0*rb[0]+V1*rb[1]+V2*rb[2]);
                            } else  ps.AddFace(*fp,(V0+V1+V2)/3.0);

                if( j < n_samples_per_edge-i-2 )
                            {
                                    CoordType V3((i+1)*segmentLen,(j+1)*segmentLen, 1.0 - ((i+1)*segmentLen+(j+1)*segmentLen) ) ;
                                    n_samples++;
                                    if(randomFlag) 	{
                                                        CoordType rb=RandomBarycentric();
                                                        ps.AddFace(*fp, V3*rb[0]+V1*rb[1]+V2*rb[2]);
                                        } else  ps.AddFace(*fp,(V3+V1+V2)/3.0);
                            }
        }
    return n_samples;
}

// Similar sampling
// Each triangle is subdivided into similar triangles following a generalization of the classical 1-to-4 splitting rule of triangles.
// According to the level of subdivision <k> you get 1, 4 , 9, 16 , <k^2> triangles.
// Depending on the kind of the sampling strategies we can have two different approach to choosing the sample points.
// 1) you have already sampled both edges and vertices
// 2) you are not going to take samples on edges and vertices.
//
// In the first case you have to consider only internal vertices of the subdivided triangles (to avoid multiple sampling of edges and vertices).
// Therefore the number of internal points is ((k-3)*(k-2))/2. where k is the number of points on an edge (vertex included)
// E.g. for k=4 you get 3 segments on each edges and the original triangle is subdivided
// into 9 smaller triangles and you get (1*2)/2 == 1 only a single internal point.
// So if you want N samples in a triangle you have to solve  k^2 -5k +6 - 2N = 0
// from which you get:
//
//      5 + sqrt( 1 + 8N )
// k = -------------------
//             2
//
// In the second case if you are not interested to skip the sampling on edges and vertices you have to consider as sample number the number of triangles.
// So if you want N samples in a triangle, the number <k> of points on  an edge (vertex included) should be simply:
//      k = 1 + sqrt(N)
// examples:
// N = 4 -> k = 3
// N = 9 -> k = 4



//template <class MeshType>
//void Sampling<MeshType>::SimilarFaceSampling()
static void FaceSimilar(MeshType & m, VertexSampler &ps,int sampleNum, bool dualFlag, bool randomFlag)
{
        ScalarType area = Stat<MeshType>::ComputeMeshArea(m);
        ScalarType samplePerAreaUnit = sampleNum/area;

        // Similar Triangles sampling.
    int n_samples_per_edge;
    double  n_samples_decimal = 0.0;
    FaceIterator fi;

    for(fi=m.face.begin(); fi != m.face.end(); fi++)
    {
        // compute # samples in the current face.
        n_samples_decimal += 0.5*DoubleArea(*fi) * samplePerAreaUnit;
        int n_samples          = (int) n_samples_decimal;
        if(n_samples>0)
        {
            // face sampling.
            if(dualFlag)
                            {
                                    n_samples_per_edge = (int)((sqrt(1.0+8.0*(double)n_samples) +5.0)/2.0); // original for non dual case
                                    n_samples = SingleFaceSimilar(&*fi,ps, n_samples_per_edge);
                            } else {
                                    n_samples_per_edge = (int)(sqrt((double)n_samples) +1.0);
                                    n_samples = SingleFaceSimilarDual(&*fi,ps, n_samples_per_edge,randomFlag);
                        }
        }
        n_samples_decimal -= (double) n_samples;
    }
}


    // Rasterization fuction
    // Take a triangle
    // T deve essere una classe funzionale che ha l'operatore ()
    // con due parametri x,y di tipo S esempio:
    // class Foo { public void operator()(int x, int y ) { ??? } };

// This function does rasterization with a safety buffer area, thus accounting some points actually outside triangle area
// The safety area samples are generated according to face flag BORDER which should be true for texture space border edges
// Use correctSafePointsBaryCoords = true to map safety texels to closest point barycentric coords (on edge).
    static void SingleFaceRaster(typename MeshType::FaceType &f,  VertexSampler &ps,
                            const Point2<typename MeshType::ScalarType> & v0,
                            const Point2<typename MeshType::ScalarType> & v1,
                            const Point2<typename MeshType::ScalarType> & v2,
                            bool correctSafePointsBaryCoords=true)
    {
    typedef typename MeshType::ScalarType S;
    // Calcolo bounding box
    Box2i bbox;
    Box2<S> bboxf;
    bboxf.Add(v0);
    bboxf.Add(v1);
    bboxf.Add(v2);

    bbox.min[0] = floor(bboxf.min[0]);
    bbox.min[1] = floor(bboxf.min[1]);
    bbox.max[0] = ceil(bboxf.max[0]);
    bbox.max[1] = ceil(bboxf.max[1]);

    // Calcolo versori degli spigoli
    Point2<S> d10 = v1 - v0;
    Point2<S> d21 = v2 - v1;
    Point2<S> d02 = v0 - v2;

    // Preparazione prodotti scalari
    S b0  = (bbox.min[0]-v0[0])*d10[1] - (bbox.min[1]-v0[1])*d10[0];
    S b1  = (bbox.min[0]-v1[0])*d21[1] - (bbox.min[1]-v1[1])*d21[0];
    S b2  = (bbox.min[0]-v2[0])*d02[1] - (bbox.min[1]-v2[1])*d02[0];
    // Preparazione degli steps
    S db0 = d10[1];
    S db1 = d21[1];
    S db2 = d02[1];
    // Preparazione segni
    S dn0 = -d10[0];
    S dn1 = -d21[0];
    S dn2 = -d02[0];

    //Calculating orientation
    bool flipped = !(d02 * vcg::Point2<S>(-d10[1], d10[0]) >= 0);

    // Calculating border edges
    Segment2<S> borderEdges[3];
    S edgeLength[3];
    unsigned char edgeMask = 0;

    if (f.IsB(0)) {
        borderEdges[0] = Segment2<S>(v0, v1);
        edgeLength[0] = borderEdges[0].Length();
        edgeMask |= 1;
    }
    if (f.IsB(1)) {
        borderEdges[1] = Segment2<S>(v1, v2);
        edgeLength[1] = borderEdges[1].Length();
        edgeMask |= 2;
    }
    if (f.IsB(2)) {
        borderEdges[2] = Segment2<S>(v2, v0);
        edgeLength[2] = borderEdges[2].Length();
        edgeMask |= 4;
    }

    // Rasterizzazione
    double de = v0[0]*v1[1]-v0[0]*v2[1]-v1[0]*v0[1]+v1[0]*v2[1]-v2[0]*v1[1]+v2[0]*v0[1];

    for(int x=bbox.min[0]-1;x<=bbox.max[0]+1;++x)
    {
        bool in = false;
        S n[3]  = { b0-db0-dn0, b1-db1-dn1, b2-db2-dn2};
        for(int y=bbox.min[1]-1;y<=bbox.max[1]+1;++y)
        {
            if( ((n[0]>=0 && n[1]>=0 && n[2]>=0) || (n[0]<=0 && n[1]<=0 && n[2]<=0))  && (de != 0))
            {
                typename MeshType::CoordType baryCoord;
                baryCoord[0] =  double(-y*v1[0]+v2[0]*y+v1[1]*x-v2[0]*v1[1]+v1[0]*v2[1]-x*v2[1])/de;
                baryCoord[1] = -double( x*v0[1]-x*v2[1]-v0[0]*y+v0[0]*v2[1]-v2[0]*v0[1]+v2[0]*y)/de;
                baryCoord[2] = 1-baryCoord[0]-baryCoord[1];

                ps.AddTextureSample(f, baryCoord, Point2i(x,y), 0);
                in = true;
            } else {
                // Check whether a pixel outside (on a border edge side) triangle affects color inside it
                Point2<S> px(x, y);
                Point2<S> closePoint;
                int closeEdge = -1;
                S minDst = FLT_MAX;

                // find the closest point (on some edge) that lies on the 2x2 squared neighborhood of the considered point
                for (int i=0; i<3; ++i)
                {
                    if (edgeMask & (1 << i))
                    {
                        Point2<S> close;
                        S dst;
                        if ( ((!flipped) && (n[i]<0)) ||
                             (  flipped  && (n[i]>0))   )
                        {
                            dst = ((close = ClosestPoint(borderEdges[i], px)) - px).Norm();
                            if(dst < minDst &&
                               close.X() > px.X()-1 && close.X() < px.X()+1 &&
                               close.Y() > px.Y()-1 && close.Y() < px.Y()+1)
                            {
                                minDst = dst;
                                closePoint = close;
                                closeEdge = i;
                            }
                        }
                    }
                }

                if (closeEdge >= 0)
                {
                    typename MeshType::CoordType baryCoord;
                    if (correctSafePointsBaryCoords)
                    {
                        // Add x,y sample with closePoint barycentric coords (on edge)
                        baryCoord[closeEdge] = (closePoint - borderEdges[closeEdge].P1()).Norm()/edgeLength[closeEdge];
                        baryCoord[(closeEdge+1)%3] = 1 - baryCoord[closeEdge];
                        baryCoord[(closeEdge+2)%3] = 0;
                    } else {
                        // Add x,y sample with his own barycentric coords (off edge)
                        baryCoord[0] =  double(-y*v1[0]+v2[0]*y+v1[1]*x-v2[0]*v1[1]+v1[0]*v2[1]-x*v2[1])/de;
                        baryCoord[1] = -double( x*v0[1]-x*v2[1]-v0[0]*y+v0[0]*v2[1]-v2[0]*v0[1]+v2[0]*y)/de;
                        baryCoord[2] = 1-baryCoord[0]-baryCoord[1];
                    }
                    ps.AddTextureSample(f, baryCoord, Point2i(x,y), minDst);
                    in = true;
                }
            }
            n[0] += dn0;
            n[1] += dn1;
            n[2] += dn2;
        }
        b0 += db0;
        b1 += db1;
        b2 += db2;
    }
}

// check the radius constrain
static bool checkPoissonDisk(SampleSHT & sht, const Point3<ScalarType> & p, ScalarType radius)
{
    // get the samples closest to the given one
    std::vector<VertexType*> closests;
  typedef EmptyTMark<MeshType> MarkerVert;
  static MarkerVert mv;

    Box3f bb(p-Point3f(radius,radius,radius),p+Point3f(radius,radius,radius));
    GridGetInBox(sht, mv, bb, closests);

  ScalarType r2 = radius*radius;
    for(int i=0; i<closests.size(); ++i)
        if(SquaredDistance(p,closests[i]->cP()) < r2)
            return false;

    return true;
}

struct PoissonDiskParam
{
  PoissonDiskParam()
  {
    adaptiveRadiusFlag = false;
    bestSampleChoiceFlag = true;
    bestSamplePoolSize = 10;
    radiusVariance =1;
    MAXLEVELS = 5;
    invertQuality = false;
    preGenFlag = false;
    preGenMesh = NULL;
    geodesicDistanceFlag = false;
    randomSeed = 0;
  }

  struct Stat
  {
    int montecarloTime;
    int gridTime;
    int pruneTime;
    int totalTime;
    Point3i gridSize;
    int gridCellNum;
    size_t sampleNum;
    int montecarloSampleNum;
  };

  bool geodesicDistanceFlag;
  bool bestSampleChoiceFlag; // In poisson disk pruning when we choose a sample in a cell, we choose the sample that remove the minimal number of other samples. This previlege the "on boundary" samples.
  int bestSamplePoolSize;
  bool adaptiveRadiusFlag;
  float radiusVariance;
  bool invertQuality;
  bool preGenFlag;            // when generating a poisson distribution, you can initialize the set of computed points with
                              // ALL the vertices of another mesh. Useful for building progressive//prioritize refinements.
  MeshType *preGenMesh;      // There are two ways of passing the pregen vertexes to the pruning, 1) is with a mesh pointer
                              // 2) with a per vertex attribute.
  int MAXLEVELS;
  int randomSeed;

  Stat pds;
};


// generate Poisson-disk sample using a set of pre-generated samples (with the Montecarlo algorithm)
// It always return a point.
static VertexPointer getSampleFromCell(Point3i &cell, MontecarloSHT & samplepool)
{
    MontecarloSHTIterator cellBegin, cellEnd;
    samplepool.Grid(cell, cellBegin, cellEnd);
    return *cellBegin;
}

// Given a cell of the grid it search the point that remove the minimum number of other samples
// it linearly scan all the points of a cell.

static VertexPointer getBestPrecomputedMontecarloSample(Point3i &cell, MontecarloSHT & samplepool, ScalarType diskRadius, const PoissonDiskParam &pp)
{
  MontecarloSHTIterator cellBegin,cellEnd;
  samplepool.Grid(cell, cellBegin, cellEnd);
  VertexPointer bestSample=0;
  int minRemoveCnt = std::numeric_limits<int>::max();
  std::vector<typename MontecarloSHT::HashIterator> inSphVec;
  int i=0;
  for(MontecarloSHTIterator ci=cellBegin; ci!=cellEnd && i<pp.bestSamplePoolSize; ++ci,i++)
  {
    VertexPointer sp = *ci;
    if(pp.adaptiveRadiusFlag)  diskRadius = sp->Q();
    int curRemoveCnt = samplepool.CountInSphere(sp->cP(),diskRadius,inSphVec);
    if(curRemoveCnt < minRemoveCnt)
    {
      bestSample = sp;
      minRemoveCnt = curRemoveCnt;
    }
  }
  return bestSample;
}

/// \brief Estimate the radius r that you should give to get a certain number of samples in a Poissson Disk Distribution of radius r.
///
static ScalarType ComputePoissonDiskRadius(MeshType &origMesh, int sampleNum)
{
    ScalarType meshArea = Stat<MeshType>::ComputeMeshArea(origMesh);
    // Manage approximately the PointCloud Case, use the half a area of the bbox.
    // TODO: If you had the radius a much better approximation could be done.
    if(meshArea ==0)
        {
                    meshArea = (origMesh.bbox.DimX()*origMesh.bbox.DimY() +
                                            origMesh.bbox.DimX()*origMesh.bbox.DimZ() +
                                            origMesh.bbox.DimY()*origMesh.bbox.DimZ());
        }
    ScalarType diskRadius = sqrt(meshArea / (0.7 * M_PI * sampleNum)); // 0.7 is a density factor
    return diskRadius;
}

static int ComputePoissonSampleNum(MeshType &origMesh, ScalarType diskRadius)
{
    ScalarType meshArea = Stat<MeshType>::ComputeMeshArea(origMesh);
    int sampleNum = meshArea /  (diskRadius*diskRadius *M_PI *0.7)  ; // 0.7 is a density factor
    return sampleNum;
}

/// When performing an adptive pruning for each sample we expect a varying radius to be removed.
/// The radius is a PerVertex attribute that we compute from the current quality
///
/// the expected radius of the sample is computed so that
/// it linearly maps the quality between diskradius and diskradius*variance
/// in other words the radius

static void InitRadiusHandleFromQuality(MeshType &sampleMesh, PerVertexFloatAttribute &rH, ScalarType diskRadius, ScalarType radiusVariance, bool invert)
{
  std::pair<float,float> minmax = tri::Stat<MeshType>::ComputePerVertexQualityMinMax( sampleMesh);
  float minRad = diskRadius ;
  float maxRad = diskRadius * radiusVariance;
  float deltaQ = minmax.second-minmax.first;
  float deltaRad = maxRad-minRad;
  for (VertexIterator vi = sampleMesh.vert.begin(); vi != sampleMesh.vert.end(); vi++)
    {
      rH[*vi] = minRad + deltaRad*((invert ? minmax.second - (*vi).Q() : (*vi).Q() - minmax.first )/deltaQ);
    }
}

// initialize spatial hash table for searching
// radius is the radius of empty disk centered over the samples (e.g. twice of the empty space disk)
// This radius implies that when we pick a sample in a cell all that cell probably will not be touched again.
// Howvever we must ensure that we do not put too many vertices inside each hash cell

static void InitSpatialHashTable(MeshType &montecarloMesh, MontecarloSHT &montecarloSHT, ScalarType diskRadius,
                                 struct PoissonDiskParam pp=PoissonDiskParam())
{
  ScalarType cellsize = 2.0f* diskRadius / sqrt(3.0);
  float occupancyRatio=0;
  do
  {
    // inflating
    BoxType bb=montecarloMesh.bbox;
    assert(!bb.IsNull());
    bb.Offset(cellsize);

    int sizeX = std::max(1,int(bb.DimX() / cellsize));
    int sizeY = std::max(1,int(bb.DimY() / cellsize));
    int sizeZ = std::max(1,int(bb.DimZ() / cellsize));
    Point3i gridsize(sizeX, sizeY, sizeZ);

    montecarloSHT.InitEmpty(bb, gridsize);

    for (VertexIterator vi = montecarloMesh.vert.begin(); vi != montecarloMesh.vert.end(); vi++)
      if(!(*vi).IsD())
      {
        montecarloSHT.Add(&(*vi));
      }

    montecarloSHT.UpdateAllocatedCells();
    pp.pds.gridSize = gridsize;
    pp.pds.gridCellNum = (int)montecarloSHT.AllocatedCells.size();
    cellsize/=2.0f;
    occupancyRatio = float(montecarloMesh.vn) / float(montecarloSHT.AllocatedCells.size());
    //    qDebug(" %i / %i = %6.3f", montecarloMesh.vn , montecarloSHT.AllocatedCells.size(),occupancyRatio);
  }
  while( occupancyRatio> 100);
}

static void PoissonDiskPruningByNumber(VertexSampler &ps, MeshType &m,
                                       size_t sampleNum, ScalarType &diskRadius,
                                       PoissonDiskParam &pp,
                                       float tolerance=0.04,
                                       int maxIter=20)

{
  size_t sampleNumMin = int(float(sampleNum)*(1.0f-tolerance));
  size_t sampleNumMax = int(float(sampleNum)*(1.0f+tolerance));
  float RangeMinRad = m.bbox.Diag()/50.0;
  float RangeMaxRad = m.bbox.Diag()/50.0;
  size_t RangeMinRadNum;
  size_t RangeMaxRadNum;
   // Note     RangeMinRad          <       RangeMaxRad
  //  but      RangeMinRadNum > sampleNum > RangeMaxRadNum
  do {
    ps.reset();
    RangeMinRad/=2.0f;
    PoissonDiskPruning(ps, m ,RangeMinRad,pp);
    RangeMinRadNum = pp.pds.sampleNum;
//    qDebug("PoissonDiskPruning Iteratin Min (%6.3f:%5i) instead of %i",RangeMinRad,RangeMinRadNum,sampleNum);
  } while(RangeMinRadNum < sampleNum); // if the number of sample is still smaller you have to make radius larger.

  do {
    ps.reset();
    RangeMaxRad*=2.0f;
    PoissonDiskPruning(ps, m ,RangeMaxRad,pp);
    RangeMaxRadNum =  pp.pds.sampleNum;
//    qDebug("PoissonDiskPruning Iteratin Max (%6.3f:%5i) instead of %i",RangeMaxRad,RangeMaxRadNum,sampleNum);
  } while(RangeMaxRadNum > sampleNum);


  float curRadius=RangeMaxRad;
  int iterCnt=0;
  while(iterCnt<maxIter &&
        (pp.pds.sampleNum < sampleNumMin || pp.pds.sampleNum  > sampleNumMax) )
  {
    iterCnt++;
    ps.reset();
    curRadius=(RangeMaxRad+RangeMinRad)/2.0f;
    PoissonDiskPruning(ps, m ,curRadius,pp);
//    qDebug("PoissonDiskPruning Iteratin (%6.3f:%5lu %6.3f:%5lu) Cur Radius %f -> %lu sample instead of %lu",RangeMinRad,RangeMinRadNum,RangeMaxRad,RangeMaxRadNum,curRadius,pp.pds.sampleNum,sampleNum);
    if(pp.pds.sampleNum > sampleNum){
      RangeMinRad = curRadius;
      RangeMinRadNum = pp.pds.sampleNum;
    }
    if(pp.pds.sampleNum < sampleNum){
      RangeMaxRad = curRadius;
      RangeMaxRadNum =  pp.pds.sampleNum;
    }
  }
  diskRadius = curRadius;
}


/// This is the main function that is used to build a poisson distribuition
/// starting from a dense sample cloud (the montecarloMesh) by 'pruning' it.
/// it puts all the samples in a hashed UG and randomly choose a sample
/// and remove all the points in the sphere centered on the chosen sample
/// 
/// You can impose some constraint: all the vertices in the montecarloMesh 
/// that are marked with a bool attribute called "fixed" are surely chosen 
/// (if you also set the  preGenFlag option)
/// 
static void PoissonDiskPruning(VertexSampler &ps, MeshType &montecarloMesh,
                               ScalarType diskRadius, PoissonDiskParam &pp)
{
  tri::RequireCompactness(montecarloMesh);
  if(pp.randomSeed) SamplingRandomGenerator().initialize(pp.randomSeed);
  if(pp.adaptiveRadiusFlag)
    tri::RequirePerVertexQuality(montecarloMesh);
  int t0 = clock();
    // spatial index of montecarlo samples - used to choose a new sample to insert
    MontecarloSHT montecarloSHT;
    InitSpatialHashTable(montecarloMesh,montecarloSHT,diskRadius,pp);

    // if we are doing variable density sampling we have to prepare the handle that keeps the the random samples expected radii.
    // At this point we just assume that there is the quality values as sampled from the base mesh
    PerVertexFloatAttribute rH = tri::Allocator<MeshType>:: template GetPerVertexAttribute<float> (montecarloMesh,"radius");
    if(pp.adaptiveRadiusFlag)
        InitRadiusHandleFromQuality(montecarloMesh, rH, diskRadius, pp.radiusVariance, pp.invertQuality);

    unsigned int (*p_myrandom)(unsigned int) = RandomInt;
    std::random_shuffle(montecarloSHT.AllocatedCells.begin(),montecarloSHT.AllocatedCells.end(), p_myrandom);
    int t1 = clock();
    pp.pds.montecarloSampleNum = montecarloMesh.vn;
    pp.pds.sampleNum =0;
    int removedCnt=0;
    // Initial pass for pruning the Hashed grid with the an eventual pre initialized set of samples
    if(pp.preGenFlag)
    {
      if(pp.preGenMesh==0)
      {
        typename MeshType::template PerVertexAttributeHandle<bool> fixed;
        fixed = tri::Allocator<MeshType>:: template GetPerVertexAttribute<bool> (montecarloMesh,"fixed");
        for(VertexIterator vi=montecarloMesh.vert.begin();vi!=montecarloMesh.vert.end();++vi)
          if(fixed[*vi]) {
            pp.pds.sampleNum++;
            ps.AddVert(*vi);
            removedCnt += montecarloSHT.RemoveInSphere(vi->cP(),diskRadius);
          }
      }
      else
      {
        for(VertexIterator vi =pp.preGenMesh->vert.begin(); vi!=pp.preGenMesh->vert.end();++vi)
        {
          ps.AddVert(*vi);
          pp.pds.sampleNum++;
          removedCnt += montecarloSHT.RemoveInSphere(vi->cP(),diskRadius);
        }
      }
      montecarloSHT.UpdateAllocatedCells();
    }
    vertex::ApproximateGeodesicDistanceFunctor<VertexType> GDF;
    while(!montecarloSHT.AllocatedCells.empty())
    {
        removedCnt=0;
        for (size_t i = 0; i < montecarloSHT.AllocatedCells.size(); i++)
        {
            if( montecarloSHT.EmptyCell(montecarloSHT.AllocatedCells[i])  ) continue;
            ScalarType currentRadius =diskRadius;
            VertexPointer sp;
            if(pp.bestSampleChoiceFlag)
              sp = getBestPrecomputedMontecarloSample(montecarloSHT.AllocatedCells[i], montecarloSHT, diskRadius, pp);
            else
              sp = getSampleFromCell(montecarloSHT.AllocatedCells[i], montecarloSHT);

            if(pp.adaptiveRadiusFlag)
              currentRadius = rH[sp];

            ps.AddVert(*sp);
            pp.pds.sampleNum++;
            if(pp.geodesicDistanceFlag) removedCnt += montecarloSHT.RemoveInSphereNormal(sp->cP(),sp->cN(),GDF,currentRadius);
                            else        removedCnt += montecarloSHT.RemoveInSphere(sp->cP(),currentRadius);
        }
        montecarloSHT.UpdateAllocatedCells();
    }
    int t2 = clock();
    pp.pds.gridTime = t1-t0;
    pp.pds.pruneTime = t2-t1;
}

/** Compute a Poisson-disk sampling of the surface.
 *  The radius of the disk is computed according to the estimated sampling density.
 *
 * This algorithm is an adaptation of the algorithm of White et al. :
 *
 * "Poisson Disk Point Set by Hierarchical Dart Throwing"
 * K. B. White, D. Cline, P. K. Egbert,
 * IEEE Symposium on Interactive Ray Tracing, 2007,
 * 10-12 Sept. 2007, pp. 129-132.
 */
static void HierarchicalPoissonDisk(MeshType &origMesh, VertexSampler &ps, MeshType &montecarloMesh, ScalarType diskRadius, const struct PoissonDiskParam pp=PoissonDiskParam())
{
//  int t0=clock();
    // spatial index of montecarlo samples - used to choose a new sample to insert
    MontecarloSHT montecarloSHTVec[5];



    // initialize spatial hash table for searching
    // radius is the radius of empty disk centered over the samples (e.g. twice of the empty space disk)
    // This radius implies that when we pick a sample in a cell all that cell will not be touched again.
    ScalarType cellsize = 2.0f* diskRadius / sqrt(3.0);

    // inflating
    BoxType bb=origMesh.bbox;
    bb.Offset(cellsize);

  int sizeX = std::max(1.0f,bb.DimX() / cellsize);
  int sizeY = std::max(1.0f,bb.DimY() / cellsize);
  int sizeZ = std::max(1.0f,bb.DimZ() / cellsize);
    Point3i gridsize(sizeX, sizeY, sizeZ);

    // spatial hash table of the generated samples - used to check the radius constrain
    SampleSHT checkSHT;
    checkSHT.InitEmpty(bb, gridsize);


    // sampling algorithm
    // ------------------
    //
    // - generate millions of samples using montecarlo algorithm
    // - extract a cell (C) from the active cell list (with probability proportional to cell's volume)
    // - generate a sample inside C by choosing one of the contained pre-generated samples
    //   - if the sample violates the radius constrain discard it, and add the cell to the cells-to-subdivide list
    // - iterate until the active cell list is empty or a pre-defined number of subdivisions is reached
    //

    int level = 0;

    // initialize spatial hash to index pre-generated samples
    montecarloSHTVec[0].InitEmpty(bb, gridsize);
    // create active cell list
    for (VertexIterator vi = montecarloMesh.vert.begin(); vi != montecarloMesh.vert.end(); vi++)
        montecarloSHTVec[0].Add(&(*vi));
    montecarloSHTVec[0].UpdateAllocatedCells();

  // if we are doing variable density sampling we have to prepare the random samples quality with the correct expected radii.
    PerVertexFloatAttribute rH = tri::Allocator<MeshType>:: template GetPerVertexAttribute<float> (montecarloMesh,"radius");
    if(pp.adaptiveRadiusFlag)
            InitRadiusHandleFromQuality(montecarloMesh, rH, diskRadius, pp.radiusVariance, pp.invertQuality);

    do
    {
        MontecarloSHT &montecarloSHT = montecarloSHTVec[level];

        if(level>0)
        {// initialize spatial hash with the remaining points
            montecarloSHT.InitEmpty(bb, gridsize);
            // create active cell list
            for (typename MontecarloSHT::HashIterator hi = montecarloSHTVec[level-1].hash_table.begin(); hi != montecarloSHTVec[level-1].hash_table.end(); hi++)
                montecarloSHT.Add((*hi).second);
            montecarloSHT.UpdateAllocatedCells();
        }
        // shuffle active cells
        unsigned int (*p_myrandom)(unsigned int) = RandomInt;
        std::random_shuffle(montecarloSHT.AllocatedCells.begin(),montecarloSHT.AllocatedCells.end(), p_myrandom);

        // generate a sample inside C by choosing one of the contained pre-generated samples
        //////////////////////////////////////////////////////////////////////////////////////////
    int removedCnt=montecarloSHT.hash_table.size();
    int addedCnt=checkSHT.hash_table.size();
        for (int i = 0; i < montecarloSHT.AllocatedCells.size(); i++)
        {
            for(int j=0;j<4;j++)
            {
                if( montecarloSHT.EmptyCell(montecarloSHT.AllocatedCells[i])  ) continue;

            // generate a sample chosen from the pre-generated one
            typename MontecarloSHT::HashIterator hi = montecarloSHT.hash_table.find(montecarloSHT.AllocatedCells[i]);

            if(hi==montecarloSHT.hash_table.end()) {break;}
            VertexPointer sp = (*hi).second;
            // vr spans between 3.0*r and r / 4.0 according to vertex quality
            ScalarType sampleRadius = diskRadius;
            if(pp.adaptiveRadiusFlag)  sampleRadius = rH[sp];
            if (checkPoissonDisk(checkSHT, sp->cP(), sampleRadius))
            {
               ps.AddVert(*sp);
               montecarloSHT.RemoveCell(sp);
               checkSHT.Add(sp);
               break;
            }
            else
                montecarloSHT.RemovePunctual(sp);
        }
        }
        addedCnt = checkSHT.hash_table.size()-addedCnt;
        removedCnt = removedCnt-montecarloSHT.hash_table.size();

        // proceed to the next level of subdivision
        // increase grid resolution
        gridsize *= 2;

        //
        level++;
    } while(level < 5);
}

//template <class MeshType>
//void Sampling<MeshType>::SimilarFaceSampling()

// This function also generates samples outside faces if those affects faces in texture space.
// Use correctSafePointsBaryCoords = true to map safety texels to closest point barycentric coords (on edge)
// otherwise obtained samples will map to barycentric coord actually outside face
//
// If you don't need to get those extra points clear faces Border Flags
// vcg::tri::UpdateFlags<Mesh>::FaceClearB(m);
//
// Else make sure to update border flags from texture space FFadj
// vcg::tri::UpdateTopology<Mesh>::FaceFaceFromTexCoord(m);
// vcg::tri::UpdateFlags<Mesh>::FaceBorderFromFF(m);
static void Texture(MeshType & m, VertexSampler &ps, int textureWidth, int textureHeight, bool correctSafePointsBaryCoords=true)
{
typedef Point2<ScalarType> Point2x;
        printf("Similar Triangles face sampling\n");
        for(FaceIterator fi=m.face.begin(); fi != m.face.end(); fi++)
            if (!fi->IsD())
            {
                Point2x ti[3];
                for(int i=0;i<3;++i)
                    ti[i]=Point2x((*fi).WT(i).U() * textureWidth - 0.5, (*fi).WT(i).V() * textureHeight - 0.5);
                    // - 0.5 constants are used to obtain correct texture mapping

                SingleFaceRaster(*fi,  ps, ti[0],ti[1],ti[2], correctSafePointsBaryCoords);
            }
}

typedef GridStaticPtr<FaceType, ScalarType > TriMeshGrid;

class RRParam
{
public:
float offset;
float minDiag;
tri::FaceTmark<MeshType> markerFunctor;
TriMeshGrid gM;
};

static void RegularRecursiveOffset(MeshType & m, std::vector<CoordType> &pvec, ScalarType offset, float minDiag)
{
    Box3<ScalarType> bb=m.bbox;
    bb.Offset(offset*2.0);

    RRParam rrp;

    rrp.markerFunctor.SetMesh(&m);

    rrp.gM.Set(m.face.begin(),m.face.end(),bb);


    rrp.offset=offset;
    rrp.minDiag=minDiag;
    SubdivideAndSample(m, pvec, bb, rrp, bb.Diag());
}

static void SubdivideAndSample(MeshType & m, std::vector<CoordType> &pvec, const Box3<ScalarType> bb, RRParam &rrp, float curDiag)
{
  CoordType startPt = bb.Center();

  ScalarType dist;
  // Compute mesh point nearest to bb center
  FaceType   *nearestF=0;
  ScalarType dist_upper_bound = curDiag+rrp.offset;
  CoordType closestPt;
  vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
  dist=dist_upper_bound;
  nearestF =  rrp.gM.GetClosest(PDistFunct,rrp.markerFunctor,startPt,dist_upper_bound,dist,closestPt);
  curDiag /=2;
  if(dist < dist_upper_bound)
  {
    if(curDiag/3 < rrp.minDiag) //store points only for the last level of recursion (?)
    {
      if(rrp.offset==0)
        pvec.push_back(closestPt);
      else
      {
        if(dist>rrp.offset) // points below the offset threshold cannot be displaced at the right offset distance, we can only make points nearer.
        {
          CoordType delta = startPt-closestPt;
          pvec.push_back(closestPt+delta*(rrp.offset/dist));
        }
      }
    }
    if(curDiag < rrp.minDiag) return;
    CoordType hs = (bb.max-bb.min)/2;
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
        for(int k=0;k<2;k++)
          SubdivideAndSample(m, pvec,
                             BoxType(CoordType( bb.min[0]+i*hs[0],  bb.min[1]+j*hs[1],  bb.min[2]+k*hs[2]),
                                     CoordType(startPt[0]+i*hs[0], startPt[1]+j*hs[1], startPt[2]+k*hs[2]) ),
                             rrp,curDiag
              );

  }
}
}; // end sampling class

template <class MeshType>
typename MeshType::ScalarType ComputePoissonDiskRadius(MeshType &origMesh, int sampleNum)
{
  typedef typename MeshType::ScalarType ScalarType;
  ScalarType meshArea = Stat<MeshType>::ComputeMeshArea(origMesh);
  // Manage approximately the PointCloud Case, use the half a area of the bbox.
  // TODO: If you had the radius a much better approximation could be done.
  if(meshArea ==0)
  {
    meshArea = (origMesh.bbox.DimX()*origMesh.bbox.DimY() +
                origMesh.bbox.DimX()*origMesh.bbox.DimZ() +
                origMesh.bbox.DimY()*origMesh.bbox.DimZ());
  }
  ScalarType diskRadius = sqrt(meshArea / (0.7 * M_PI * sampleNum)); // 0.7 is a density factor
  return diskRadius;
}



template <class MeshType>
void MontecarloSampling(MeshType &m, // the mesh that has to be sampled
                     MeshType &mm, // the mesh that will contain the samples
                     int sampleNum) // the desired number sample, if zero you must set the radius to the wanted value
{
  typedef tri::MeshSampler<MeshType> BaseSampler;
  MeshSampler<MeshType> mcSampler(&mm);
  tri::SurfaceSampling<MeshType,BaseSampler>::Montecarlo(m, mcSampler, sampleNum);
}


template <class MeshType>
void MontecarloSampling(MeshType &m, // the mesh that has to be sampled
                     std::vector<Point3f> &montercarloSamples, // the vector that will contain the set of points
                     int sampleNum) // the desired number sample, if zero you must set the radius to the wanted value
{
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  BaseSampler mcSampler(montercarloSamples);
  tri::SurfaceSampling<MeshType,BaseSampler>::Montecarlo(m, mcSampler, sampleNum);
}

// Yet another simpler wrapper for the generation of a poisson disk distribution over a mesh.
//
template <class MeshType>
void PoissonSampling(MeshType &m, // the mesh that has to be sampled
                     std::vector<typename MeshType::CoordType> &poissonSamples, // the vector that will contain the set of points
                     int sampleNum, // the desired number sample, if zero you must set the radius to the wanted value
                     typename MeshType::ScalarType &radius,  // the Poisson Disk Radius (used if sampleNum==0, setted if sampleNum!=0)
                     typename MeshType::ScalarType radiusVariance=1,
                     typename MeshType::ScalarType PruningByNumberTolerance=0.04f,
                     unsigned int randSeed=0)

{
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  typedef tri::MeshSampler<MeshType> MontecarloSampler;

  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;
  int t0=clock();

//  if(sampleNum>0) radius = tri::SurfaceSampling<MeshType,BaseSampler>::ComputePoissonDiskRadius(m,sampleNum);
  if(radius>0 && sampleNum==0) sampleNum = tri::SurfaceSampling<MeshType,BaseSampler>::ComputePoissonSampleNum(m,radius);

  pp.pds.sampleNum = sampleNum;
  pp.randomSeed = randSeed;
  poissonSamples.clear();
//  std::vector<Point3f> MontecarloSamples;
  MeshType MontecarloMesh;

  // First step build the sampling
  MontecarloSampler mcSampler(MontecarloMesh);
  BaseSampler pdSampler(poissonSamples);

  if(randSeed) tri::SurfaceSampling<MeshType,MontecarloSampler>::SamplingRandomGenerator().initialize(randSeed);
  tri::SurfaceSampling<MeshType,MontecarloSampler>::Montecarlo(m, mcSampler, std::max(10000,sampleNum*40));
  tri::UpdateBounding<MeshType>::Box(MontecarloMesh);
//  tri::Build(MontecarloMesh, MontecarloSamples);
  int t1=clock();
  pp.pds.montecarloTime = t1-t0;

  if(radiusVariance !=1)
  {
    pp.adaptiveRadiusFlag=true;
    pp.radiusVariance=radiusVariance;
  }
  if(sampleNum==0) tri::SurfaceSampling<MeshType,BaseSampler>::PoissonDiskPruning(pdSampler, MontecarloMesh, radius,pp);
  else tri::SurfaceSampling<MeshType,BaseSampler>::PoissonDiskPruningByNumber(pdSampler, MontecarloMesh, sampleNum, radius,pp,PruningByNumberTolerance);
  int t2=clock();
  pp.pds.totalTime = t2-t0;
}

/// \brief Low level wrapper for Poisson Disk Pruning
///
/// This function simply takes a mesh and a radius and returns a vector of vertex pointers listing the "surviving" points.
//
template <class MeshType>
void PoissonPruning(MeshType &m, // the mesh that has to be pruned
                    std::vector<typename MeshType::VertexPointer> &poissonSamples, // the vector that will contain the chosen set of points
                    float radius, unsigned int randSeed=0)
{
  typedef tri::TrivialPointerSampler<MeshType> BaseSampler;
  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;
  pp.randomSeed = randSeed;

  tri::UpdateBounding<MeshType>::Box(m);
  BaseSampler pdSampler;
  tri::SurfaceSampling<MeshType,BaseSampler>::PoissonDiskPruning(pdSampler, m, radius,pp);
  std::swap(pdSampler.sampleVec,poissonSamples);
}


/// \brief Low level wrapper for Poisson Disk Pruning
///
/// This function simply takes a mesh containing a point cloud to be pruned and a radius 
/// It returns a vector of CoordType listing the "surviving" points.
///
template <class MeshType>
void PoissonPruning(MeshType &m, // the mesh that has to be pruned
                    std::vector<typename MeshType::CoordType> &poissonSamples, // the vector that will contain the chosen set of points
                    float radius, unsigned int randSeed=0)
{
  std::vector<typename MeshType::VertexPointer> poissonSamplesVP;
  PoissonPruning(m,poissonSamplesVP,radius,randSeed);
  poissonSamples.resize(poissonSamplesVP.size());
  for(size_t i=0;i<poissonSamplesVP.size();++i)
    poissonSamples[i]=poissonSamplesVP[i]->P();
}



/// \brief Very simple wrapping for the Exact Poisson Disk Pruning
///
/// This function simply takes a mesh and an expected number of points and returns
/// vector of points. It performs multiple attempts with varius radii to correctly get the expected number of samples.
/// It is obviously much slower than the other versions...
template <class MeshType>
void PoissonPruningExact(MeshType &m, /// the mesh that has to be pruned
                         std::vector<typename MeshType::VertexPointer> &poissonSamples, /// the vector that will contain the chosen set of points
                         typename MeshType::ScalarType & radius,
                         int sampleNum,
                         float tolerance=0.04,
                         int maxIter=20,
                         unsigned int randSeed=0)
{
  size_t sampleNumMin = int(float(sampleNum)*(1.0f-tolerance));  // the expected values range.
  size_t sampleNumMax = int(float(sampleNum)*(1.0f+tolerance));  // e.g. any sampling in [sampleNumMin, sampleNumMax] is OK
  float RangeMinRad = m.bbox.Diag()/10.0f;
  float RangeMaxRad = m.bbox.Diag()/10.0f;
  size_t RangeMinSampleNum;
  size_t RangeMaxSampleNum;
  std::vector<typename MeshType::VertexPointer> poissonSamplesTmp;

  do
  {
    RangeMinRad/=2.0f;
    PoissonPruning(m,poissonSamplesTmp,RangeMinRad,randSeed);
    RangeMinSampleNum = poissonSamplesTmp.size();
  } while(RangeMinSampleNum < sampleNumMin);

  do
  {
    RangeMaxRad*=2.0f;
    PoissonPruning(m,poissonSamplesTmp,RangeMaxRad,randSeed);
    RangeMaxSampleNum = poissonSamplesTmp.size();
  } while(RangeMaxSampleNum > sampleNumMax);

  float curRadius;
  int iterCnt=0;
  while(iterCnt<maxIter &&
        (poissonSamplesTmp.size() < sampleNumMin || poissonSamplesTmp.size() > sampleNumMax) )
  {
    curRadius=(RangeMaxRad+RangeMinRad)/2.0f;
    PoissonPruning(m,poissonSamplesTmp,curRadius,randSeed);
    //qDebug("(%6.3f:%5i %6.3f:%5i) Cur Radius %f -> %i sample instead of %i",RangeMinRad,RangeMinSampleNum,RangeMaxRad,RangeMaxSampleNum,curRadius,poissonSamplesTmp.size(),sampleNum);
    if(poissonSamplesTmp.size() > size_t(sampleNum))
      RangeMinRad = curRadius;
    if(poissonSamplesTmp.size() < size_t(sampleNum))
      RangeMaxRad = curRadius;
  }

  swap(poissonSamples,poissonSamplesTmp);
  radius = curRadius;
}
} // end namespace tri
} // end namespace vcg

#endif

