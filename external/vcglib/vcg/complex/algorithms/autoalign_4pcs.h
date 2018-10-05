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
#ifndef _AUTOALIGN_4PCS_H_
#define _AUTOALIGN_4PCS_H_

/**
implementation of the 4PCS method from the paper:
"4-Points Congruent Sets for Robust Pairwise Surface Registration"
D.Aiger, N.Mitra D.Cohen-Or, SIGGRAPH 2008
ps: the name of the variables are out of vcg standard but like the one
used in the paper pseudocode.
*/

#include <vcg/complex/complex.h>
#include <vcg/space/point_matching.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/math/random_generator.h>
#include <ctime>
namespace vcg{
namespace tri{

template <class MeshType>
class FourPCS {
public:
    /* mesh only for using spatial indexing functions (to remove) */
  class PVertex;    // dummy prototype never used
  class PFace;

  class PUsedTypes: public vcg::UsedTypes < vcg::Use<PVertex>::template AsVertexType,
                                            vcg::Use<PFace  >::template AsFaceType >{};

  class PVertex : public vcg::Vertex< PUsedTypes,vcg::vertex::BitFlags,vcg::vertex::Coord3f,vcg::vertex::Mark>{};
  class PFace   : public vcg::Face<   PUsedTypes> {};
  class PMesh   : public vcg::tri::TriMesh< std::vector<PVertex>, std::vector<PFace> > {};

  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename vcg::Matrix44<ScalarType> Matrix44x;
  typedef typename vcg::Box3<ScalarType> Box3x;

  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::VertexType VertexType;
  typedef vcg::Point4< vcg::Point3<ScalarType> > FourPoints;
  typedef vcg::GridStaticPtr<typename PMesh::VertexType, ScalarType > GridType;

  /* class for Parameters */
  struct Param
  {
    ScalarType overlap;    // overlap estimation as a percentage of overlapping points.

    int sampleNumP;        // number of samples on moving mesh P (it determines the sampling radius to be used to sample Q too)
    float samplingRadius;

    ScalarType deltaPerc;  // Approximation Level (expressed as a percentage of the avg distance between samples)
    ScalarType deltaAbs;   // Approximation Level
    int feetSize;          // how many points in the neighborhood of each of the 4 points
    int scoreFeet;         // how many of the feetsize points must match (max feetsize*4) to try an early interrupt
    ScalarType cosAngle;   // max admittable angle that can be admitted between matching points in alignments (expressed as cos(ang) )
    int seed;              // random seed used. Need for repeatability.

    void Default(){
      overlap = 0.5;
      sampleNumP=500;
      samplingRadius=0;

      deltaPerc = 0.5;
      deltaAbs = 0;
      feetSize = 25;
      scoreFeet = 50;
      seed =0;
      cosAngle = 0; // normals must differ more than 90 degree to be considered bad.
    }
  };

  struct Stat
  {
    Stat() : initTime(0),selectCoplanarBaseTime(0),findCongruentTime(0),testAlignmentTime(0)
    {}
    clock_t initTime;
    clock_t selectCoplanarBaseTime;
    clock_t findCongruentTime;
    clock_t testAlignmentTime;
    float init()   {return 1000.0f*float(initTime)/float(CLOCKS_PER_SEC);}
    float select() {return 1000.0f*float(selectCoplanarBaseTime)/float(CLOCKS_PER_SEC);}
    float findCongruent() {return 1000.0f*float(findCongruentTime)/float(CLOCKS_PER_SEC);}
    float testAlignment() {return 1000.0f*float(testAlignmentTime)/float(CLOCKS_PER_SEC);}
  };

  class Couple
  {
  public:
    VertexPointer p0,p1;
    Couple(VertexPointer i, VertexPointer j, float d) : p0(i),p1(j),dist(d){}
    float dist;
    bool operator < (const   Couple & o) const {return dist < o.dist;}
    VertexPointer operator[](const int &i) const {return (i==0)? this->p0 : this->p1;}
  };

  struct Candidate
  {
    Candidate():score(0){}
    Candidate(FourPoints _p, vcg::Matrix44<ScalarType>_T):p(_p),T(_T){}
    FourPoints  p;
    vcg::Matrix44<ScalarType> T;
    int score;
    inline bool operator <(const Candidate & o) const {return score > o.score;}
  };


  // class for the point  'ei'
  struct EPoint{
    EPoint(vcg::Point3<ScalarType> _p, int _i):pos(_p),pi(_i){}
    vcg::Point3<ScalarType> pos;
    int pi;        //index to R[1|2]
    void GetBBox(vcg::Box3<ScalarType> & b){b.Add(pos);}
  };


  Param par;    /// parameters
  Stat stat;

  MeshType  *P;    // Moving Mesh (from which the coplanar base is selected)
  MeshType  *Q;    // Fixed Mesh  (mesh where to find the correspondences)

  math::MarsenneTwisterRNG rnd;

  std::vector<VertexPointer> subsetQ;  // subset of the vertices in Q
  std::vector<VertexPointer> subsetP; // random selection on P

  ScalarType side;               // side

  PMesh     Invr;                // invariants
  std::vector< Candidate > U;    // the
  int iwinner;                   // winner == U[iwinner]
  std::vector<FourPoints> bases; // used bases
  std::vector<VertexType*> ExtB[4]; // selection of vertices "close" to the four point

  vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType > ugridQ;
  vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType > ugridP;

  /* returns the closest point between to segments x1-x2 and x3-x4.  */
    void IntersectionLineLine(const CoordType & x1,const CoordType & x2,const CoordType & x3,const CoordType & x4, CoordType&x)
    {
      CoordType a = x2-x1, b = x4-x3, c = x3-x1;
      x = x1 + a * ((c^b).dot(a^b)) / (a^b).SquaredNorm();
    }


void Init(MeshType &_movP,MeshType &_fixQ)
{
  clock_t t0= clock();
  P = &_movP;
  Q = &_fixQ;
  tri::UpdateBounding<MeshType>::Box(*P);
  if(par.seed==0) rnd.initialize(time(0));
  else rnd.initialize(par.seed);

  ugridQ.Set(Q->vert.begin(),Q->vert.end());
  ugridP.Set(P->vert.begin(),P->vert.end());

  if(par.samplingRadius==0)
    par.samplingRadius = tri::ComputePoissonDiskRadius(*P,par.sampleNumP);
  tri::PoissonPruning(*P, subsetP, par.samplingRadius, par.seed);
  tri::PoissonPruning(*Q, subsetQ, par.samplingRadius, par.seed);
  par.deltaAbs = par.samplingRadius * par.deltaPerc;
  side = P->bbox.Dim()[P->bbox.MaxDim()]*par.overlap; //rough implementation
  stat.initTime+=clock()-t0;
}


// Try to select four coplanar points such that they are at least side distance and
//
bool SelectCoplanarBase(FourPoints &B, ScalarType &r1, ScalarType &r2)
{
  clock_t t0= clock();

  // choose the inter point distance
  ScalarType dtol = side*0.1; //rough implementation


  // **** first point: random
  B[0] = P->vert[ rnd.generate(P->vert.size())].P();

  // **** second point: a random point at distance side +-dtol
  size_t i;
  for(i = 0; i < P->vert.size(); ++i){
    int id = rnd.generate(P->vert.size());
    ScalarType dd = (P->vert[id].P() - B[0]).Norm();
    if(  ( dd < side + dtol) && (dd > side - dtol)){
      B[1] = P->vert[id].P();
      break;
    }
  }
  if(i ==  P->vert.size()) return false;

  // **** third point: at distance less than side*0.8 from middle way between B[0] and B[1]
  const CoordType middle = (B[0]+B[1])/2.0;
  for(i = 0; i < P->vert.size(); ++i){
    int id = rnd.generate(P->vert.size());
    if( Distance(P->vert[id].P(),middle) < side*0.8 ){
      B[2] = P->vert[id].P();
      break;
    }
  }
  if(i ==  P->vert.size()) return false;

  // **** fourth point:
  ScalarType cpr = rnd.generate01();
  CoordType crossP = B[0] *(1-cpr)+B[1]*cpr;
  CoordType B4 = B[2]+(crossP-B[2]).Normalize()*side;
  CoordType n = ((B[0]-B[1]).normalized() ^ (B[2]-B[1]).normalized()).normalized();
  ScalarType radius = dtol;

  std::vector<typename MeshType::VertexType*> closests;
  std::vector<ScalarType> distances;
  std::vector<CoordType> points;

  vcg::tri::GetInSphereVertex<
      MeshType,
      vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType >,
      std::vector<typename MeshType::VertexType*>,
      std::vector<ScalarType>,
      std::vector<CoordType>
      >(*P,ugridP,B4,radius,closests,distances,points);

  if(closests.empty())
    return false;
  int bestInd = -1;   ScalarType bestv=std::numeric_limits<float>::max();
  for(i = 0; i <closests.size(); ++i){
    ScalarType dist_from_plane = fabs((closests[i]->P() - B[1]).normalized().dot(n));
    if( dist_from_plane < bestv){
      bestv = dist_from_plane;
      bestInd = i;
    }
  }
  if(bestv >dtol)
    return false;
  B[3] =  closests[bestInd]->P();

  //printf("B[3] %d\n", (typename MeshType::VertexType*)closests[best] - &(*P->vert.begin()));

  // compute r1 and r2
  CoordType x;
  //        std::swap(B[1],B[2]);
  IntersectionLineLine(B[0],B[1],B[2],B[3],x);

  r1 = (x - B[0]).dot(B[1]-B[0]) / (B[1]-B[0]).SquaredNorm();
  r2 = (x - B[2]).dot(B[3]-B[2]) / (B[3]-B[2]).SquaredNorm();

  if( ((B[0]+(B[1]-B[0])*r1)-(B[2]+(B[3]-B[2])*r2)).Norm() > par.deltaAbs )
    return false;

  radius  = side*0.5;
  std::vector< CoordType > samples;
  std::vector<ScalarType > dists;

  for(int i  = 0 ; i< 4; ++i){
    vcg::tri::GetKClosestVertex<
        MeshType,
        vcg::GridStaticPtr<typename MeshType::VertexType, ScalarType >,
        std::vector<VertexType*>,
        std::vector<ScalarType>,
        std::vector< CoordType > >(*P,ugridP, par.feetSize ,B[i],radius, ExtB[i], dists, samples);
  }

  //qDebug("ExtB %i",ExtB[0].size()+ExtB[1].size()+ExtB[2].size()+ExtB[3].size());
  stat.selectCoplanarBaseTime+=clock()-t0;
  return true;
}

bool IsTransfCongruent(const FourPoints &B, const FourPoints &fp, vcg::Matrix44<ScalarType> & mat)
{
  std::vector<vcg::Point3<ScalarType> > fix(4);
  std::vector<vcg::Point3<ScalarType> > mov(4);
  for(int i = 0 ; i < 4; ++i) {
    mov[i]=B[i];
    fix[i]=fp[i];
  }

  if(fabs( Distance(fix[0],fix[1]) - Distance(mov[0],mov[1]) ) > par.deltaAbs) return false;
  if(fabs( Distance(fix[0],fix[2]) - Distance(mov[0],mov[2]) ) > par.deltaAbs) return false;
  if(fabs( Distance(fix[0],fix[3]) - Distance(mov[0],mov[3]) ) > par.deltaAbs) return false;
  if(fabs( Distance(fix[1],fix[2]) - Distance(mov[1],mov[2]) ) > par.deltaAbs) return false;
  if(fabs( Distance(fix[1],fix[3]) - Distance(mov[1],mov[3]) ) > par.deltaAbs) return false;
  if(fabs( Distance(fix[2],fix[3]) - Distance(mov[2],mov[3]) ) > par.deltaAbs) return false;

  vcg::ComputeRigidMatchMatrix(fix,mov,mat);

  ScalarType maxSquaredDistance = 0.0;
  for(int i = 0; i < 4; ++i)
    maxSquaredDistance =std::max(maxSquaredDistance, SquaredDistance(mat * mov[i] ,fix[i]));
  return  sqrt(maxSquaredDistance)  < par.deltaAbs;
}

/// Compute the vector R1 of couple of points on FixQ at a given distance.
/// Used by FindCongruent
void ComputeR1(std::vector<Couple > &R1)
{
  R1.clear();
  for(size_t vi = 0; vi  < subsetQ.size(); ++vi)
    for(size_t vj = vi; vj < subsetQ.size(); ++vj){
      ScalarType d = Distance(subsetQ[vi]->P(),subsetQ[vj]->P());
      if( (d < side+par.deltaAbs))
      {
        R1.push_back(Couple(subsetQ[vi],subsetQ[vj], d));
        R1.push_back(Couple(subsetQ[vj],subsetQ[vi], d));
      }
    }

  std::sort(R1.begin(),R1.end());
}

// Find congruent elements of a base B, on Q, with approximation delta
// and put them in the U vector.
bool FindCongruent(const std::vector<Couple > &R1, const FourPoints &B, const ScalarType r1, const ScalarType r2)
{
  clock_t t0=clock();
  int n_base=0;
  bool done = false;
  int n_closests = 0, n_congr = 0;
  int ac =0 ,acf = 0,tr = 0,trf =0;
  ScalarType d1,d2;
  d1 = (B[1]-B[0]).Norm();
  d2 = (B[3]-B[2]).Norm();

  typename std::vector<Couple>::const_iterator bR1,eR1,bR2,eR2,ite;
  bR1 = std::lower_bound(R1.begin(),R1.end(),Couple(0,0,d1-par.deltaAbs));
  eR1 = std::lower_bound(R1.begin(),R1.end(),Couple(0,0,d1+par.deltaAbs));
  bR2 = std::lower_bound(R1.begin(),R1.end(),Couple(0,0,d2-par.deltaAbs));
  eR2 = std::lower_bound(R1.begin(),R1.end(),Couple(0,0,d2+par.deltaAbs));

  // in  [bR1,eR1) there are all the pairs at a distance d1 +- par.delta
  // in  [bR1,eR1) there are all the pairs at a distance d2 +- par.delta

  if(bR1 == R1.end()) return false;// if there are no such pairs return
  if(bR2 == R1.end()) return false; // if there are no such pairs return

  // put [bR1,eR1) in a mesh to have the search operator for free (lazy me)
  Invr.Clear();
  typename PMesh::VertexIterator vii;
  int i = &(*bR1)-&(*R1.begin());
  for(ite = bR1; ite != eR1;++ite){
    vii = vcg::tri::Allocator<PMesh>::AddVertices(Invr,1);
    //      (*vii).P() = Q->vert[R1[i][0]].P() + (Q->vert[R1[i][1]].P()-Q->vert[R1[i][0]].P()) * r1;
    (*vii).P() .Import(         ite->p0->P() + (        ite->p1->P() -       ite->p0->P()) * r1);
    ++i;
  }
  if(Invr.vert.empty() ) return false;

  // per vertex attribute 'index' remaps a vertex of Invr to its corresponding point in R1
  typename PMesh::template PerVertexAttributeHandle<int> id = vcg::tri::Allocator<PMesh>::template AddPerVertexAttribute<int>(Invr,std::string("index"));
  i = &(*bR1)-&(*R1.begin());
  for(vii = Invr.vert.begin(); vii != Invr.vert.end();++vii,++i)  id[vii] = i;

  vcg::tri::UpdateBounding<PMesh>::Box(Invr);


  std::vector<EPoint> R2inv;
  i = &(*bR2)-&(*R1.begin());
  // R2inv contains all the points generated by the couples in R2 (with the reference to remap into R2)
  for(ite = bR2; ite != eR2;++ite){
    //        R2inv.push_back( EPoint( Q->vert[R1[i][0]].P() + (Q->vert[R1[i][1]].P()-Q->vert[R1[i][0]].P()) * r2,i));
    R2inv.push_back( EPoint( R1[i].p0->P() + (R1[i].p1->P() - R1[i].p0->P()) * r2,i));
    ++i;
  }

  GridType ugrid; // griglia
  ugrid.Set(Invr.vert.begin(),Invr.vert.end());
  n_closests = 0; n_congr = 0; ac =0 ; acf = 0; tr = 0; trf = 0;
  printf("R2Inv.size  = %d \n",R2inv.size());
  for(unsigned int i = 0 ; i < R2inv.size() ; ++i)
  {
    std::vector<typename PMesh::VertexType*> closests;

    // for each point in R2inv get all the points in R1 closer than par.delta
    vcg::Matrix44<ScalarType> mat;
    Box3x bb;
    bb.Add(R2inv[i].pos+CoordType(par.deltaAbs,par.deltaAbs, par.deltaAbs));
    bb.Add(R2inv[i].pos-CoordType(par.deltaAbs,par.deltaAbs, par.deltaAbs));

    vcg::tri::GetInBoxVertex<PMesh,GridType,std::vector<typename PMesh::VertexType*> >
        (Invr,ugrid,bb,closests);

    if(closests.size() > 5)
      closests.resize(5);

    n_closests+=closests.size();
    for(unsigned int ip = 0; ip < closests.size(); ++ip)
    {
      FourPoints p;
      p[0] = R1[id[closests[ip]]][0]->cP();
      p[1] = R1[id[closests[ip]]][1]->cP();
      p[2] = R1[ R2inv[i].pi][0]->cP();
      p[3] = R1[ R2inv[i].pi][1]->cP();

      n_base++;
      if(!IsTransfCongruent(B,p,mat)) {
        trf++;
      }
      else{
        tr++;
        n_congr++;
        Candidate c(p,mat);
        EvaluateAlignment(c);

        if( c.score > par.scoreFeet)
          U.push_back(c);
      }
    }
  }

  vcg::tri::Allocator<PMesh>::DeletePerVertexAttribute(Invr,id);
  printf("n_closests %5d = (An %5d ) + ( Tr %5d ) + (OK) %5d\n",n_closests,acf,trf,n_congr);

  stat.findCongruentTime += clock()-t0;
  return done;
}


int EvaluateSample(Candidate & fp, const CoordType & tp, const CoordType & np)
{
  CoordType ttp = fp.T * tp;
  vcg::Point4<ScalarType> np4 = fp.T * vcg::Point4<ScalarType>(np[0],np[1],np[2],0.0);
  CoordType tnp(np4[0],np4[1],np4[2]);

  ScalarType   dist ;
  VertexType* v = vcg::tri::GetClosestVertex(*Q, ugridQ, ttp, par.deltaAbs*2.0,  dist  );

  if(v!=0)
  {
    if( v->N().dot(tnp) > par.cosAngle )  return 1;
    else return -1;
  }
  else return 0;
}

// Check a candidate against the small subset of points ExtB
void EvaluateAlignment(Candidate  & fp){
        int n_delta_close = 0;
        for(int i  = 0 ; i< 4; ++i) {
            for(unsigned int j = 0; j < ExtB[i].size();++j){
                n_delta_close+=EvaluateSample(fp, ExtB[i][j]->P(), ExtB[i][j]->cN());
            }
        }
        fp.score = n_delta_close;
}

void TestAlignment(Candidate  & fp)
{
  clock_t t0 = clock();
  int n_delta_close = 0;
  for(unsigned int j = 0; j < subsetP.size();++j){
    CoordType np = subsetP[j]->N();
    CoordType tp = subsetP[j]->P();
    n_delta_close+=EvaluateSample(fp,tp,np);
  }
  fp.score =  n_delta_close;
  stat.testAlignmentTime += clock()-t0;
}


bool Align(Matrix44x & result, vcg::CallBackPos * cb )
{
  int maxAttempt =100;
  int scoreThr = par.sampleNumP*0.8;

  Candidate bestC;

  std::vector<Couple > R1;
  ComputeR1(R1);
  for(int i  = 0; i  < maxAttempt && bestC.score<scoreThr ; ++i )
  {
    FourPoints B;
    ScalarType r1,r2;
    if(SelectCoplanarBase(B,r1,r2))
    {
      U.clear();
      FindCongruent(R1,B,r1,r2);
      qDebug("Attempt %i found %i candidate best score %i",i,U.size(),bestC.score);
      for(size_t i = 0 ; i <  U.size() ;++i)
      {
        TestAlignment(U[i]);
        if(U[i].score > bestC.score)
          bestC = U[i];
      }
    }
  }
  result = bestC.T;
  return bestC.score >0;
}

bool Align(int L, Matrix44x & result, vcg::CallBackPos * cb )
{
    int bestv = 0;
    bool found;
    int n_tries = 0;
    U.clear();

    if(L==0)
    {
        // overlap is expressed as the probability that a point in P(mov) can be found in Q (fix)
        L = (log(1.0-0.9) / log(1.0-pow((float)par.overlap,3.f)))+1;
        printf("using %d bases\n",L);
    }
    std::vector<Couple > R1;
    ComputeR1(R1);

    for(int t  = 0; t  < L; ++t )
    {
      FourPoints B;
      ScalarType r1,r2;
      do
      {
        n_tries = 0;
        do
        {
          n_tries++;
          found = SelectCoplanarBase(B,r1,r2);
        }
        while(!found && (n_tries < 50));
        if(!found) {
          par.overlap*=0.9;
          side = P->bbox.Dim()[P->bbox.MaxDim()]*par.overlap; //rough implementation
          ComputeR1(R1);
        }
      } while (!found && (par.overlap >0.1));

      if(par.overlap < 0.1) {
        printf("FAILED");
        return false;
      }
      bases.push_back(B);
      if(cb) cb(t*100/L,"Trying bases");
      if(FindCongruent(R1,B,r1,r2))
        break;
    }

    if(U.empty()) return false;

//    std::sort(U.begin(),U.end());
    if(cb) cb(90,"TestAlignment");
    bestv  = -std::numeric_limits<float>::max();
    iwinner = 0;

    for(int i = 0 ; i <  U.size() ;++i)
     {
        TestAlignment(U[i]);
        if(U[i].score > bestv){
            bestv = U[i].score;
            iwinner = i;
            }
    }

    result = U[iwinner].T;
    Invr.Clear();
    return true;
}

}; // end class

} // namespace tri
} // namespace vcg
#endif
