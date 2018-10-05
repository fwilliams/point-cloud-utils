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
#ifndef __VCG_HISTOGRAM
#define __VCG_HISTOGRAM

#include <vector>
#include <algorithm>
#include <assert.h>
#include <string>
#include <limits>
#include <vcg/math/base.h>
#include <stdio.h>

namespace vcg {

template <class ScalarType>
class Distribution
{
private:
  std::vector<ScalarType> vec;
  bool dirty;
  double valSum;
  double sqrdValSum;
  double avg;
  double sqrdAvg;
  double rms;
  double min_v;
  double max_v;


public:

  Distribution() { Clear(); }

  void Clear()
  {
    vec.clear();
    dirty=true;
    min_v =  std::numeric_limits<float>::max();
    max_v = -std::numeric_limits<float>::max();
  }

  void Add(const ScalarType v)
  {
    vec.push_back(v);
    dirty=true;
    if(v<min_v) min_v=v;
    if(v>max_v) max_v=v;
  }

  ScalarType Min() const { return min_v; }
  ScalarType Max() const { return max_v; }
  ScalarType Cnt() const { return ScalarType(vec.size()); }

  ScalarType Sum(){ DirtyCheck(); return valSum; }
  ScalarType Avg(){ DirtyCheck(); return avg;}
  //! Returns the Root Mean Square of the data.
  ScalarType RMS(){ DirtyCheck(); return rms;}

  //! \brief Returns the variance of the data.
  /// the average of the squares less the square of the average.
  ScalarType Variance(){ DirtyCheck(); return sqrdAvg - avg*avg ;}

  //! Returns the standard deviation of the data.
  ScalarType StandardDeviation(){ DirtyCheck(); return sqrt( Variance() );}

  void DirtyCheck()
  {
    if(!dirty) return;
    std::sort(vec.begin(),vec.end());
    valSum=0;
    sqrdValSum=0;
    typename std::vector<ScalarType>::iterator vi;
    for(vi=vec.begin();vi!=vec.end();++vi)
    {
      valSum += double(*vi);
      sqrdValSum += double(*vi)*double(*vi);
    }
    avg = valSum/double(vec.size());
    sqrdAvg = sqrdValSum/double(vec.size());
    rms = math::Sqrt(sqrdAvg);
    dirty=false;
  }

  ScalarType Percentile(ScalarType perc)
  {
    assert(!vec.empty());
    assert(perc>=0 && perc<=1);
    DirtyCheck();
    int index = vec.size() *perc -1;

    if(index< 0 ) index = 0;

    return vec[index];
  }
};



/**
 * Histogram.
 *
 * This class implements a single-value histogram.
 */
template <class ScalarType>
class Histogram
{

  // public data members
protected:

  std::vector <ScalarType> H; 	//! Counters for bins.
  std::vector <ScalarType> R; 	//! Range for bins.
  ScalarType minv; 	//! Minimum value.
  ScalarType maxv;	//! Maximum value.
  ScalarType minElem; 	//! Minimum value.
  ScalarType maxElem;	//! Maximum value.
  int n;	//! Number of vaild intervals stored between minv and maxv.


  /// incrementally updated values
  ScalarType cnt;	//! Number of accumulated samples.
  ScalarType sum;	//! Average.
  ScalarType rms; 	//! Root mean square.

  /**
        * Returns the index of the bin which contains a given value.
        */
  int BinIndex(ScalarType val) ;

  // public methods
public:

  /**
     * Set the histogram values.
     *
     * This method is used to correctly initialize the bins of the histogram.
     * n is the number of valid intervals between minv and maxv.
     * for a more robust working, the Histogram class stores also the two out of range intervals (-inf, minv] and [maxv, +inf)
     * Each bin is left closed (eg it contains the value
     * The \a gamma parameter is applied to modify the distribution of the ranges of the bins. Default uniform distibution.
     *
     */
  void SetRange(ScalarType _minv, ScalarType _maxv, int _n,ScalarType gamma=1.0 );

  ScalarType MinV() {return minv;} 	//! Minimum value of the range where the histogram is defined.
  ScalarType MaxV() {return maxv;} 	//! Maximum value of the range where the histogram is defined.
  ScalarType Sum()  {return sum;} 	//! Total sum of inserted values.
  ScalarType Cnt()  {return cnt;}

  ScalarType MinElem() {return minElem;} 	//! Minimum element that has been added to the histogram. It could be < or > than MinV;.
  ScalarType MaxElem() {return maxElem;} 	//! Maximum element that has been added to the histogram. It could be < or > than MinV;..

  /**
     * Add a new value to the histogram.
     *
     * The statistics related to the histogram data (average, RMS, etc.) are
     * also updated.
     */
  void Add(ScalarType v, ScalarType increment=ScalarType(1.0));

  ScalarType MaxCount() const;        //! Max number of elements among all buckets (including the two infinity bounded buckets)
  ScalarType MaxCountInRange() const; //! Max number of elements among all buckets between MinV and MaxV.
  int BinNum() const {return n;}
  ScalarType BinCount(ScalarType v);
  ScalarType BinCountInd(int index) {return H[index];}
  ScalarType BinCount(ScalarType v, ScalarType width);
  ScalarType BinLowerBound(int index) {return R[index];}
  ScalarType BinUpperBound(int index) {return R[index+1];}
  ScalarType RangeCount(ScalarType rangeMin, ScalarType rangeMax);
  ScalarType BinWidth(ScalarType v);

  /**
     * Returns the value corresponding to a given percentile of the data.
     *
     * The percentile range between 0 and 1.
     */
  ScalarType Percentile(ScalarType frac) const;

  //! Returns the average of the data.
  ScalarType Avg(){ return sum/cnt;}

  //! Returns the Root Mean Square of the data.
  ScalarType RMS(){ return sqrt(rms/double(cnt));}

  //! Returns the variance of the data.
  ScalarType Variance(){ return fabs(rms/cnt-Avg()*Avg());}

  //! Returns the standard deviation of the data.
  ScalarType StandardDeviation(){ return sqrt(Variance());}

  //! Dump the histogram to a file.
  void FileWrite(const std::string &filename);

  //! Reset histogram data.
  void Clear();
};

template <class ScalarType>
void Histogram<ScalarType>::Clear()
{
  H.clear();
  R.clear();
  cnt=0;
  sum=0;
  rms=0;
  n=0;
  minv=0;
  maxv=1;
  minElem = std::numeric_limits<ScalarType>::max();
  maxElem = -std::numeric_limits<ScalarType>::max();
}

/*
Note that the histogram holds <n> valid bins plus two semi-infinite bins.

R[0]   = -inf
R[1]   = minv
R[n+1] = maxv
R[n+2] = +inf


Eg. SetRange(0, 10, 5) asks for 5 intervals covering the 0..10 range

    H[0]  H[1]   H[2]   H[3]   H[4]   H[5]    H[6]
-inf    0      2      4      6      8      10    +inf
R[0]   R[1]   R[2]   R[3]   R[4]   R[5]   R[6]   R[7]

*/

template <class ScalarType>
void Histogram<ScalarType>::SetRange(ScalarType _minv, ScalarType _maxv, int _n, ScalarType gamma)
{
  // reset data
  Clear();

  minv=_minv;maxv=_maxv;n=_n;
  H.resize(n+2);
  fill(H.begin(),H.end(),0);
  R.resize(n+3);

  R[0]   = - std::numeric_limits< ScalarType >::max();
  R[n+2] =   std::numeric_limits< ScalarType >::max();

  double delta=(maxv-minv);
  if(gamma==1)
  {
    for(int i=0; i<=n; ++i)
      R[i+1] = minv + delta*ScalarType(i)/n;
  }
  else
  {
    for(int i=0; i<=n; ++i)
      R[i+1] = minv + delta*pow(ScalarType(i)/n,gamma);
  }
}



template <class ScalarType>
int Histogram<ScalarType>::BinIndex(ScalarType val)
{
  // lower_bound returns the furthermost iterator i in [first, last) such that, for every iterator j in [first, i), *j < value.
  // E.g. An iterator pointing to the first element "not less than" val, or end() if every element is less than val.
  typename std::vector<ScalarType>::iterator it = lower_bound(R.begin(),R.end(),val);

  assert(it!=R.begin());
  assert(it!=R.end());
  assert((*it)>=val);

  int pos = it-R.begin();
  assert(pos >=1);
  pos -= 1;
  assert (R[pos] < val);
  assert (         val <= R[pos+1] );
  return pos;
}

/*
    H[0]  H[1]   H[2]   H[3]   H[4]   H[5]    H[6]
-inf    0      2      4      6      8      10    +inf
R[0]   R[1]   R[2]   R[3]   R[4]   R[5]   R[6]   R[7]

asking for 3.14 lower bound will return an iterator pointing to R[3]==4; and will increase H[2]
asking for 4    lower bound will return an iterator pointing to R[3]==4; and will increase H[2]

*/
template <class ScalarType>
void Histogram<ScalarType>::Add(ScalarType v, ScalarType increment)
{
  int pos=BinIndex(v);
  if(v<minElem) minElem=v;
  if(v>maxElem) maxElem=v;
  assert((pos>=0)&&(pos<=n+1));
  H[pos]+=increment;
  cnt+=increment;
  sum+=v*increment;
  rms += (v*v)*increment;
}

template <class ScalarType>
ScalarType Histogram<ScalarType>::BinCount(ScalarType v)
{
  return H[BinIndex(v)];
}

template <class ScalarType>
ScalarType Histogram<ScalarType>::BinCount(ScalarType v, ScalarType width)
{
  return RangeCount(v-width/2.0,v+width/2.0);
}

template <class ScalarType>
ScalarType Histogram<ScalarType>::RangeCount(ScalarType rangeMin, ScalarType rangeMax)
{
  int firstBin=BinIndex(rangeMin);
  int lastBin=BinIndex (rangeMax);
  ScalarType sum=0;
  for(int i=firstBin; i<=lastBin;++i)
    sum+=H[i];
  return sum;
}

template <class ScalarType>
ScalarType Histogram<ScalarType>::BinWidth(ScalarType v)
{
  int pos=BinIndex(v);
  return R[pos+1]-R[pos];
}

template <class ScalarType>
void Histogram<ScalarType>::FileWrite(const std::string &filename)
{
  FILE *fp;
  fp=fopen(filename.c_str(),"w");

  for(unsigned int i=0; i<H.size(); i++)
    fprintf (fp,"%12.8lf , %12.8lf \n",R[i],double(H[i])/cnt);

  fclose(fp);
}


template <class ScalarType>
ScalarType Histogram<ScalarType>::MaxCount() const
{
  return *(std::max_element(H.begin(),H.end()));
}

template <class ScalarType>
ScalarType Histogram<ScalarType>::MaxCountInRange() const
{
  return *(std::max_element(H.begin()+1,H.end()-1));
}

// Return the scalar value <r> such that there are <frac> samples <= <r>.
// E.g. Percentile(0.0) will return  R[1]   e.g. min value
// E.g. Percentile(1.0) will return  R[n+1] e.g max value

template <class ScalarType>
ScalarType Histogram<ScalarType>::Percentile(ScalarType frac) const
{
  if(H.size()==0 && R.size()==0)
    return 0;

  // check percentile range
  assert(frac >= 0 && frac <= 1);

  ScalarType sum=0,partsum=0;
  size_t i;

  // useless summation just to be sure
  for(i=0;i<H.size();i++) sum+=H[i];
  assert(sum==cnt);

  sum*=frac;
  for(i=0; i<H.size(); i++)
  {
    partsum+=H[i];
    if(partsum>=sum) break;
  }

  assert(i<H.size());

  return R[i+1];
}

typedef Histogram<double> Histogramd ;
typedef Histogram<float> Histogramf ;

} // end namespace (vcg)

#endif  /* __VCG_HISTOGRAM */
