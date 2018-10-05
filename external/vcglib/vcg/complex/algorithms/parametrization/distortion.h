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

#ifndef VCG_PARAM_DISTORTION
#define VCG_PARAM_DISTORTION
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

namespace vcg {
namespace tri{
template <class MeshType, bool PerWedgeFlag>


/*
 *  Energy types:
 *
 *    AreaDist : 0 for equiareal (equipotent) mappings
 *    EdgeDist (hack): 0 for isometric mappings (computed on edges only)
 *    AngleDist (hack): 0 for conformal mappings
 *    CrossDist : as above, but computed on tangent directions (not UVs)
 *    L2Stretch : 1 for isometric mappings (averaged case on the mesh),
 *                +inf on degenerate / folded cases
 *                Described in [1]
 *    LInfStretch : as above, but WORST case
 *                  (returns the worst stretch on any position and direction)
 *                  Described in [1]
 *
 * [1] Sander, P. V., Snyder, J., Gortler, S. J., & Hoppe, H.
 *     "Texture mapping progressive meshes."
 *      In Proc. ACM SIGGRAPH (pp. 409-416). 2001
 */

class Distortion
{
public:
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType::CurVecType CurVecType;
    typedef typename MeshType::FaceType::TexCoordType::ScalarType TexScalarType;
    typedef Point2<TexScalarType> TexCoordType;

    static ScalarType Area3D(const FaceType *f)
    {
        return DoubleArea(*f)*(0.5);
    }

    static ScalarType AreaUV(const FaceType *f)
    {
        TexCoordType uv0,uv1,uv2;
        if(PerWedgeFlag) {
            uv0=f->cWT(0).P();
            uv1=f->cWT(1).P();
            uv2=f->cWT(2).P();
        } else {
            uv0=f->cV(0)->T().P();
            uv1=f->cV(1)->T().P();
            uv2=f->cV(2)->T().P();
        }
        ScalarType AreaUV=((uv1-uv0)^(uv2-uv0))/2.0;
        return AreaUV;
    }

    static ScalarType EdgeLenght3D(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
        ScalarType length=(f->cP0(e)-f->cP1(e)).Norm();
        return (length);
    }

    static ScalarType EdgeLenghtUV(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
        Point2<TexScalarType> uv0,uv1;
        if(PerWedgeFlag) {
            uv0=f->cWT(e+0).P();
            uv1=f->cWT((e+1)%3).P();
        } else {
            uv0=f->cV0(e)->T().P();
            uv1=f->cV1(e)->T().P();
        }
        ScalarType UVlength=Distance(uv0,uv1);
        return UVlength;
    }

    static ScalarType AngleCos3D(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
        CoordType p0=f->P((e+2)%3);
        CoordType p1=f->P(e);
        CoordType p2=f->P((e+1)%3);
        typedef typename CoordType::ScalarType ScalarType;
        CoordType dir0=p2-p1;
        CoordType dir1=p0-p1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType angle=dir0*dir1;
        return angle;
    }

    static ScalarType AngleCosUV(const FaceType *f,int e)
    {
        Point2<ScalarType> uv0,uv1,uv2;
        if(PerWedgeFlag) {
            uv0=f->cWT((e+2)%3).P();
            uv1=f->cWT((e+0)%3).P();
            uv2=f->cWT((e+1)%3).P();
        } else {
            uv0=f->V2(e)->T().P();
            uv1=f->V0(e)->T().P();
            uv2=f->V1(e)->T().P();
        }
        vcg::Point2<ScalarType> dir0=uv2-uv1;
        vcg::Point2<ScalarType> dir1=uv0-uv1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType angle=dir0*dir1;
        return angle;
    }

    static ScalarType AngleRad3D(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
        CoordType p0=f->cP((e+2)%3);
        CoordType p1=f->cP(e);
        CoordType p2=f->cP((e+1)%3);
        typedef typename CoordType::ScalarType ScalarType;
        CoordType dir0=p2-p1;
        CoordType dir1=p0-p1;
        return Angle(dir0,dir1);
    }

    static ScalarType AngleRadUV(const FaceType *f,int e)
    {
        Point2<TexScalarType> uv0,uv1,uv2;
        if(PerWedgeFlag) {
            uv0=f->cWT((e+2)%3).P();
            uv1=f->cWT((e+0)%3).P();
            uv2=f->cWT((e+1)%3).P();
        } else {
            uv0=f->cV2(e)->T().P();
            uv1=f->cV0(e)->T().P();
            uv2=f->cV1(e)->T().P();
        }
        vcg::Point2<TexScalarType> dir0=uv2-uv1;
        vcg::Point2<TexScalarType> dir1=uv0-uv1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType t=dir0*dir1;
        if(t>1) t = 1;
        else if(t<-1) t = -1;
        return acos(t);
    }


public:
    enum DistType{AreaDist,EdgeDist,AngleDist,CrossDist,L2Stretch,LInfStretch};

    ///return the absolute difference between angle in 3D space and texture space
    ///Actually the difference in cos space
    static ScalarType AngleCosDistortion(const FaceType *f,int e)
    {
        ScalarType Angle_3D=AngleCos3D(f,e);
        ScalarType Angle_UV=AngleCosUV(f,e);
        ScalarType diff=fabs(Angle_3D-Angle_UV);///Angle_3D;
        return diff;
    }
    ///return the absolute difference between angle in 3D space and texture space
    ///Actually the difference in cos space
    static ScalarType AngleRadDistortion(const FaceType *f,int e)
    {
        ScalarType Angle_3D=AngleRad3D(f,e);
        ScalarType Angle_UV=AngleRadUV(f,e);
        ScalarType diff=fabs(Angle_3D-Angle_UV)/Angle_3D;///Angle_3D;
        return diff;
    }

    ///return the variance of angle, normalized
    ///in absolute value
    static ScalarType AngleDistortion(const FaceType *f)
    {
        return  (AngleRadDistortion(f,0) +
                AngleRadDistortion(f,1) +
                AngleRadDistortion(f,2))/3.0;
    }

    ///return the global scaling factors  from 3D to UV
    static void MeshScalingFactor(const MeshType &m,
                                  ScalarType &AreaScale,
                                  ScalarType &EdgeScale)
    {
        ScalarType SumArea3D=0;
        ScalarType SumArea2D=0;
        ScalarType SumEdge3D=0;
        ScalarType SumEdge2D=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            SumArea3D+=Area3D(&m.face[i]);
            SumArea2D+=AreaUV(&m.face[i]);
            for (int j=0;j<3;j++)
            {
                SumEdge3D+=EdgeLenght3D(&m.face[i],j);
                SumEdge2D+=EdgeLenghtUV(&m.face[i],j);
            }
        }
        AreaScale=SumArea3D/SumArea2D;
        EdgeScale=SumEdge3D/SumEdge2D;
    }

    ///return the variance of edge length, normalized in absolute value,
    ///the needed scaling factor EdgeScaleVal may be calculated
    ///by using the ScalingFactor function
    static ScalarType EdgeDistortion(const FaceType *f,int e,
                                     ScalarType EdgeScaleVal)
    {
        ScalarType edgeUV=EdgeLenghtUV(f,e)*EdgeScaleVal;
        ScalarType edge3D=EdgeLenght3D(f,e);
        assert(edge3D > 0);
        ScalarType diff=fabs(edge3D-edgeUV)/edge3D;
        assert(!math::IsNAN(diff));
        return diff;
    }

    ///return the variance of area, normalized
    ///in absolute value, the scalar AreaScaleVal may be calculated
    ///by using the ScalingFactor function
    static ScalarType AreaDistortion(const FaceType *f,
                                     ScalarType AreaScaleVal)
    {
        ScalarType areaUV=AreaUV(f)*AreaScaleVal;
        ScalarType area3D=Area3D(f);
        assert(area3D > 0);
        ScalarType diff=fabs(areaUV-area3D)/area3D;
        assert(!math::IsNAN(diff));
        return diff;
    }

    static ScalarType L2StretchEnergySquared(const FaceType *f,
                                             ScalarType AreaScaleVal)
    {
        TexCoordType p0 = (PerWedgeFlag)? f->cWT(0).P() :  f->cV(0)->T().P() ;
        TexCoordType p1 = (PerWedgeFlag)? f->cWT(1).P() :  f->cV(1)->T().P() ;
        TexCoordType p2 = (PerWedgeFlag)? f->cWT(2).P() :  f->cV(2)->T().P() ;

        CoordType q0 = f->cP(0);
        CoordType q1 = f->cP(1);
        CoordType q2 = f->cP(2);

        TexScalarType A2 = ((p1-p0)^(p2-p0));

        if (A2<0) A2 = 0; // will be NAN, +infinity

        CoordType Ss = ( q0 * ( p1[1]-p2[1] ) + q1 * (p2[1]-p0[1]) + q2 * (p0[1]-p1[1]) ) / A2;
        CoordType St = ( q0 * ( p2[0]-p1[0] ) + q1 * (p0[0]-p2[0]) + q2 * (p1[0]-p0[0]) ) / A2;

        ScalarType a = Ss.SquaredNorm() / AreaScaleVal;
        ScalarType c = St.SquaredNorm() / AreaScaleVal;

        return ((a+c)/2);
    }



    static ScalarType LInfStretchEnergy(const FaceType *f,  ScalarType AreaScaleVal)
    {
        TexCoordType p0 = (PerWedgeFlag)? f->cWT(0).P() :  f->cV(0)->T().P() ;
        TexCoordType p1 = (PerWedgeFlag)? f->cWT(1).P() :  f->cV(1)->T().P() ;
        TexCoordType p2 = (PerWedgeFlag)? f->cWT(2).P() :  f->cV(2)->T().P() ;

        CoordType q0 = f->cP(0);
        CoordType q1 = f->cP(1);
        CoordType q2 = f->cP(2);

        TexScalarType A2 = ((p1-p0)^(p2-p0));

        if (A2<0) A2 = 0; // will be NAN, +infinity

        CoordType Ss = ( q0 * ( p1[1]-p2[1] ) + q1 * (p2[1]-p0[1]) + q2 * (p0[1]-p1[1]) ) / A2;
        CoordType St = ( q0 * ( p2[0]-p1[0] ) + q1 * (p0[0]-p2[0]) + q2 * (p1[0]-p0[0]) ) / A2;

        ScalarType a = Ss.SquaredNorm() / AreaScaleVal;
        ScalarType b = Ss*St / AreaScaleVal;
        ScalarType c = St.SquaredNorm() / AreaScaleVal;

        ScalarType delta = sqrt((a-c)*(a-c)+4*b*b);
        ScalarType G =  sqrt( (a+c+delta)/2 );
        //ScalarType g = sqrt( (a+c-delta)/2 ); // not needed
        return G;
    }


    ///return the number of folded faces
    static bool Folded(const FaceType *f)
    {
        ScalarType areaUV=AreaUV(f);
        /*if (areaUV<0)
                    printf("area %5.5f \n",areaUV);*/
        return (areaUV<0);
    }

    static int Folded(const MeshType &m)
    {
        int folded=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            if(Folded(&m.face[i]))folded++;
        }
        return folded;
    }

    static bool GloballyUnFolded(const MeshType &m)
    {
        int num=Folded(m);
        return (num>(m.fn)/2);
    }

    static ScalarType MeshAngleDistortion(const MeshType &m)
    {
        ScalarType UDdist=0;
        for (int i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            const FaceType *f=&(m.face[i]);
            UDdist+=AngleDistortion(f)*Area3D(f);
        }
        return UDdist;
    }

    static ScalarType SetFQAsCrossDirDistortion(MeshType &m)
    {
        //first save the old UV dir
        std::vector<CurVecType> Dir1,Dir2;
        for (size_t i=0;i<m.face.size();i++)
        {
            Dir1.push_back(m.face[i].PD1());
            Dir2.push_back(m.face[i].PD2());
        }
        vcg::tri::CrossField<MeshType>::InitDirFromWEdgeUV(m);

        ScalarType tot = 0, totA = 0;

        //then compute angle deficit
        for (size_t i=0;i<m.face.size();i++)
        {

            FaceType &f( m.face[i] );
            CoordType transfPD1=vcg::tri::CrossField<MeshType>::K_PI(CoordType::Construct( Dir1[i] ),
                                                                     CoordType::Construct( f.PD1() ),
                                                                     f.N());
            transfPD1.Normalize();
            ScalarType AngleDeficit=vcg::Angle(transfPD1,CoordType::Construct( f.PD1() ));
            AngleDeficit=math::ToDeg(AngleDeficit);
            if ((AngleDeficit>45)||(AngleDeficit<0))
            {
                std::cout<<"Warnign A Deficit "<<AngleDeficit<<std::endl;
            }
//            assert(AngleDeficit<45);
//            assert(AngleDeficit>=0);

            ScalarType doubleArea = vcg::DoubleArea( f );
            ScalarType distortion = (AngleDeficit)/ 45 ;

            m.face[i].Q()= distortion;
            tot += distortion * doubleArea;
            totA += doubleArea;
        }

        //finally restore the original directions
        for (size_t i=0;i<m.face.size();i++)
        {
            m.face[i].PD1()=Dir1[i];
            m.face[i].PD2()=Dir2[i];
        }

        return tot / totA;
    }

    static ScalarType SetQasDistorsion(MeshType &m, DistType DType=AreaDist)
    {
        if (DType==CrossDist)
        {
            ScalarType res = SetFQAsCrossDirDistortion(m);

            vcg::tri::UpdateQuality<MeshType>::VertexFromFace(m,true);
            return res;
        }

        ScalarType edge_scale,area_scale;
        MeshScalingFactor(m,area_scale,edge_scale);

        float tot = 0;
        float totA = 0;

        for (int i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            ScalarType q;
            switch (DType) {
            case CrossDist:
                // make compiler happy
                q = 0;
                break;
            case AreaDist:
                q = AreaDistortion(&m.face[i],area_scale);
                break;
            case AngleDist:
                q = AngleDistortion(&m.face[i]);
                break;
            case EdgeDist:
                q =( EdgeDistortion(&m.face[i],0,edge_scale)+
                     EdgeDistortion(&m.face[i],1,edge_scale)+
                     EdgeDistortion(&m.face[i],2,edge_scale) )/3;
                break;
            case L2Stretch:
                q = L2StretchEnergySquared( &m.face[i],area_scale );
                break;
            case LInfStretch:
                q = LInfStretchEnergy( &m.face[i],area_scale );
                break;
            }

            m.face[i].Q() = q; // note: for L2Stretch, we are puttning E^2 on Q

            // aggregate:
            if (DType==LInfStretch) {
                tot = std::max( tot, q );
            } else {
                ScalarType a = Area3D(&m.face[i]);
                tot += q*a;
                totA += a;
            }

        }

        vcg::tri::UpdateQuality<MeshType>::VertexFromFace(m,true);

        switch (DType) {
            case L2Stretch: return sqrt(tot/totA);
            case LInfStretch: return tot;
            default:  return tot/totA;
        }
    }
};

}} // namespace end

#endif
