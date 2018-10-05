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
#ifndef VCG__SKELETON_H
#define VCG__SKELETON_H
#include<vcg/complex/algorithms/voronoi_volume_sampling.h>

namespace vcg
{
namespace tri
{

template <class MeshType> 
class SampledSkeleton
{
public:
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::BoxType BoxType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertexPointer VertexPointer;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::FacePointer FacePointer;
  typedef typename MeshType::FaceType FaceType;
  typedef VoronoiVolumeSampling<MeshType> VoronoiVolumeSamplingType;
  SampledSkeleton(VoronoiVolumeSamplingType &_vvs):vvs(_vvs){}
  
  VoronoiVolumeSamplingType &vvs;
  

/**
 * @brief Compute an evaulation of the thickness as distance from the medial axis.
 * It starts from a montecarlo volume sampling and try to search for the samples that can be part of the medial axis.
 * It use a sampled representation of the surface. A volume sample is considered part 
 * of the medial axis if there are at least two points that are (almost) the same minimal distance to that point.
 *
 * 
 */
void ThicknessEvaluator(float distThr, int smoothSize, int smoothIter, MeshType *skelM=0)
{
  tri::UpdateQuality<MeshType>::VertexConstant(vvs.psd.poissonSurfaceMesh,0);
  std::vector<VertexPointer> medialSrc(vvs.psd.poissonSurfaceMesh.vert.size(),0);
  for(VertexIterator vi=vvs.montecarloVolumeMesh.vert.begin(); vi!=vvs.montecarloVolumeMesh.vert.end(); ++vi)
   {
    unsigned int ind;
    ScalarType sqdist;
    this->vvs.psd.surfTree->doQueryClosest(vi->P(),ind,sqdist);
    VertexPointer vp = &vvs.psd.poissonSurfaceMesh.vert[ind];
    ScalarType minDist = math::Sqrt(sqdist);
    if(vp->Q() < minDist) 
    {
      std::vector<unsigned int> indVec;
      std::vector<ScalarType> sqDistVec;
      
      this->vvs.psd.surfTree->doQueryDist( vi->P(), minDist*distThr,indVec,sqDistVec);
      if(indVec.size()>1)
      {
        for(size_t i=0;i<indVec.size();++i)
        {
          VertexPointer vp = &vvs.psd.poissonSurfaceMesh.vert[indVec[i]];
          //ScalarType dist = math::Sqrt(sqDistVec[i]);
          if(vp->Q() < minDist) {
            vp->Q()=minDist;
            medialSrc[indVec[i]]=&*vi;
          }             
        }
      }       
    }
  }
  // Now collect the vertexes of the volume mesh that are on the medial surface 
  if(skelM)
  {
    tri::UpdateFlags<MeshType>::VertexClearV(vvs.montecarloVolumeMesh);
    for(size_t i=0;i<medialSrc.size();++i)
      medialSrc[i]->SetV();
    for(VertexIterator vi=vvs.montecarloVolumeMesh.vert.begin(); vi!=vvs.montecarloVolumeMesh.vert.end(); ++vi)
      if(vi->IsV()) tri::Allocator<MeshType>::AddVertex(*skelM,vi->P());
    printf("Generated a medial surf of %i vertexes\n",skelM->vn);
  }
  
  
  tri::Smooth<MeshType>::PointCloudQualityMedian(vvs.psd.poissonSurfaceMesh);
  tri::Smooth<MeshType>::PointCloudQualityAverage(vvs.psd.poissonSurfaceMesh,smoothSize,smoothIter);
  tri::UpdateColor<MeshType>::PerVertexQualityRamp(vvs.psd.poissonSurfaceMesh);
  tri::RedetailSampler<MeshType> rs;
  rs.init(&vvs.psd.poissonSurfaceMesh);
  rs.dist_upper_bound = vvs.psd.poissonSurfaceMesh.bbox.Diag()*0.05 ;
  rs.qualityFlag = true;
  tri::SurfaceSampling<MeshType, RedetailSampler<MeshType> >::VertexUniform(vvs.baseMesh, rs, vvs.baseMesh.vn, false);
}

void RefineSkeletonVolume(MeshType &skelMesh)
{
  CoordType closestP;
  int trialNum=0;
  for(int i=0;i<skelMesh.vn;++i)
   {
       CoordType point = math::GeneratePointInBox3Uniform(vvs.rng,vvs.baseMesh.bbox);
       trialNum++;
       ScalarType d = this->DistanceFromSurface(point, closestP);
       if(d<0){
         vcg::tri::Allocator<MeshType>::AddVertex(vvs.montecarloVolumeMesh,point);
         vvs.montecarloVolumeMesh.vert.back().Q() = fabs(d);
       }
   }
}


}; // end class


} // end namespace vcg
} // end namespace vcg
#endif // VCG__SKELETON_H
